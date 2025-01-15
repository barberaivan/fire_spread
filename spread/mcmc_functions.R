# Packages ---------------------------------------------------------------

library(tidyverse)
library(rstan)         # sample with HMC to check
library(LaplacesDemon) # rinvwishart
library(invgamma)
library(MASS)          # mvrnorm, not so sensitive to positive-definiteness
library(trialr)        # rlkjcorr
library(truncnorm)
library(truncreg)      # tnorm regression

# Notes ------------------------------------------------------------------

# Functions to update parameters in the mcmc run.
# Spread hyperparameters are updated with gibbs for linear regression, including
#   each sigma2;
# The correlation matrix among random effects is updated using an inverse-wishart
#   prior, but only the correlation is used. This is done to avoid the correlation
#   between scale and correlations in the inverse-wishart distribution.
# The random effects are updated through MH, but using as proposal the samples
#   from stage1. Here the stage1 prior must be subtracted, and the jacobian for the
#   steps transform must be considered in the posterior.
# The upper limit for steps (logit-normal distribution) is updated through
#   MH.
# For the area ~ steps regression, with a truncated-normal likelihood, MH updates
#   are performed.

# In this model, we assume a scaled-logit-normal distribution for the
# steps, with a logit-linear effect of the FWI (linear at logit scale).
# The lower limit for the steps is fixed at (stepsL = 2), so the spread model 
# can always spread at least one cell from the ignition cell (ignition is 
# considered as step 1). The upper limit for the logit transform is estimated
# (stepsU), but we set a flat prior between Umin and Umax. Umax is fixed at
# 2000, which implies that a fire can spread at most, 2000 * 30 = 60000 m 
# (60 km) from the ignition point. That is approximately the width of the 
# Andean-Patagonian forest. Umin is fixed at the maximum of step samples 
# obtained in stage 1. 
# stepsU will be sampled at the scaled_logit(x, Umin, Umax) scale, with a 
# logistic(0, 1) prior.

# Functions ---------------------------------------------------------------

# Transformations ---------------------------------------------------------

# function to transfrom steps (constrained) to logit scale.
logit_scaled_raw <- function(x, L, U) {
  x_unit <- (x - L) / (U - L)
  return(qlogis(x_unit))
}

# alternative form, from wolfram alpha:
logit_scaled <- function(x, L, U) {
  log((L - x) / (x - U))
}

# its derivative
logit_scaled_d <- function(x, L, U) {
  (U - L) / ((L - x) * (x - U)) 
}

# its log-derivative. Do not take abs because the transform is monotonic positive.
logit_scaled_d_log <- function(x, L, U) {
  log(U - L) - log((L - x) * (x - U)) 
}

# function to transfrom logit-steps to constrained and scaled scale.
invlogit_scaled <- function(x, L, U) {
  plogis(x) * (U - L) + L
}

# Updates -----------------------------------------------------------------

# ChatGPT helped understanding where to use jacobians, see chat here:
# https://chatgpt.com/c/6733fe6c-4018-800b-97de-d93aec6e0563

# In summary, spread-steps must be sampled at their original scale (2, stepsU),
# because their logit is a transformed parameter, not a parameter itself. 
# Hence, the Jacobian has to be included in all computations where we place a 
# prior for the logit-steps. That happens when we update stepsU, but also when
# updating the steps itself. There, to remove the prior from stage1, we must
# compute it considering its Jacobian too, as the placeholder prior was
# defined at another logit scale.

# The conditional unnormalized density for steps is

# Pr(steps | mu, sigma, eta) = 
#   Pr(eta | mu, sigma) * 
#   abs(Jacobian(logit from steps to eta)) / 
#   Pr1(steps)
# where eta = logit((steps - stepsL) / (stepsU - stepsL)) 

# Pr1 is the placeholder prior, defined as
# Pr1(steps) ~ 
#   Pr(theta | mu0, sigma0) * 
#   abs(Jacobian(logit from steps to theta))
# where theta = logit((steps - L_j) / (U_j - L_j))

# (The logit_translate will not be used.)
# log(abs(Jacobian(logit(...)))) is computed with logit_scaled_d_log()
 
# ______________________________

# Perform the gibbs update for a linear regression.
# Borrowing code from
# https://www.r-bloggers.com/2021/05/bayesian-linear-regression-with-gibbs-sampling-using-r-code/

# y: response.
# X: design matrix, including the intercept.
# tXX: t(X) %*% X, which is computed ahead of sampling.
# s2: current variance, which will be updated.
# b0: prior mean for b.
# S0_inv: the inverse of the prior vcov (usually diagonal).
# t0 and d0: prior parameters for the inverse-gamma.

# Returns a vector of updated parameters, as c(b, s2).
update_lm <- function(y, X, tXX, s2, b0, S0_inv, t0, d0) {
  # Draw b conditional on s2
  V <- solve(S0_inv + (1 / s2) * tXX)
  M <- V %*% (S0_inv %*% b0 + (1 / s2) * t(X) %*% y)
  b <- mvrnorm(1, M, V) |> as.numeric()

  # Sample sigma2 conditional on b
  resids <- y - X %*% b # residuals
  t1 <- t0 + length(y)
  d1 <- d0 + t(resids) %*% resids
  s2 <- rinvgamma(1, t1, d1)

  return(c(b, s2))
}

# Update correlation matrix for random effects errors.
# Note that raneffs must be previously substracted their expected values to
# be converted to error.
update_corr <- function(Y) {
  A <- t(Y) %*% Y

  d <- ncol(Y)
  nu_prior <- d + 1 # uniform on correlations
  S_prior <- diag(rep(1, d))

  nu_posterior <- nu_prior + nrow(Y)
  S_posterior <- S_prior + A

  R <- rinvwishart(nu_posterior, S_posterior) |> cov2cor()
  return(R)
}

# Perform the gibbs update for a linear regression with truncated normal
# likelihood.
# y: response.
# x: predictor.
# area_coef: current coefficients, to be updated (b0, b1, s2).
# L: lower truncation value.
# b0: prior mean for b.
# S0: prior vcov (usually diagonal).
# t0 and d0: prior parameters for the inverse-gamma.
# sd_jump: proposal sd for b[1], b[2] and sqrt(s2).

# Returns a vector of updated parameters, as c(b, s2).
update_truncnorm <- function(y, x, coef, L, b0, S0, t0, d0,
                             sd_jump) {
  s0 <- sqrt(diag(S0)) # prior scales for b
  b <- coef[1:2]
  s <- sqrt(coef[3])
  # _________________
  # update sigma
  # _________________
  mu <- b[1] + b[2] * x
  s_try <- rtruncnorm(1, a = 0, mean = s, sd = sd_jump[3])

  lp_new <-
    sum(dtruncnorm(y, a = L, mean = mu, sd = s_try) |> log()) +  # likelihood
    dinvgamma(s_try, t0, d0, log = T) +                          # prior
    dtruncnorm(s, a = 0, mean = s_try, sd = sd_jump[3]) |> log() # jump
  # this terms take into account the asymmetric proposal (truncnorm)

  lp_old <-
    sum(dtruncnorm(y, a = L, mean = mu, sd = s) |> log()) +      # likelihood
    dinvgamma(s, t0, d0, log = T) +                              # prior
    dtruncnorm(s_try, a = 0, mean = s, sd = sd_jump[3]) |> log() # jump

  s <- ifelse(exp(lp_new - lp_old) > runif(1), s_try, s)

  # _________________
  # update b[1]
  # _________________
  b1_try <- rnorm(1, b[1], sd_jump[1])
  mu_try <- b1_try + b[2] * x

  lp_new <-
    sum(dtruncnorm(y, a = L, mean = mu_try, sd = s) |> log()) +
    dnorm(b1_try, b0[1], s0[1], log = T)
  lp_old <-
    sum(dtruncnorm(y, a = L, mean = mu, sd = s) |> log()) +
    dnorm(b[1], b0[1], s0[1], log = T)
  b[1] <- ifelse(exp(lp_new - lp_old) > runif(1), b1_try, b[1])

  # update mu
  mu <- b[1] + b[2] * x

  # _________________
  # update b[2]
  # _________________
  b2_try <- rnorm(1, b[2], sd_jump[2])
  mu_try <- b[1] + b2_try * x

  lp_new <-
    sum(dtruncnorm(y, a = L, mean = mu_try, sd = s) |> log()) +
    dnorm(b2_try, b0[2], s0[2], log = T)
  lp_old <-
    sum(dtruncnorm(y, a = L, mean = mu, sd = s) |> log()) +
    dnorm(b[2], b0[2], s0[2], log = T)
  b[2] <- ifelse(exp(lp_new - lp_old) > runif(1), b2_try, b[2])

  # merge parameters
  return(c(b, s^2))
}

# update spread parameters (random effects).

# Y and Ytry carry the steps
# at the constrained (simulator) scale, between stepsL and stepsU. Both the 
# hierarchical and stage1 prior for steps are defined at logits, but with
# different bounds. So their jacobian is similar, but not the same.

# Steps have two logit scales: from the stage 1, each fire had steps bounds 
#   based on the landscape size, stored in steps_bounds. Then, all fires have 
#   a shared logit scale, bounded between stepsL and stepsU.

# Y: matrix of random effects, with replicates in columns. 
#   Steps are in the constrained (simulator) scale.
# mu: matrix of expected values for random effects, with the same dimension
#   as Y, but the steps-mu is at the logit (shared) scale.
# Sigma: vcov matrix for random effects, computed from scales and correlations.
# mu0: matrix with the same dimension as Y, with the prior mean from stage1 for
#   each random effect. All at logit scale, but remember that each fire here has
#   their steps at their own logit-scale, defined by the landscape size.
# Sigma0: array with the prior vcov for each random effect in each slice.
# Ytry: array with samples from stage1 for all random effects, resampled 
#   using weights, with dimension [n_coef, n_fires, n_resamples_stage1].
#   Steps are in the constrained (simulator) scale.
# steps_bounds: matrix with steps bounds to transform from constrained to
#   fire-wise logit. Fires in rows. 
# s: position of steps parameter in the parameter vector
# area: log-area for spread fires.
# areaL: lower bound for log-area.
# area_coef: b1, b2 and s2 for area, to compute likelihood conditional on steps.
update_ranef <- function(Y, mu, Sigma, mu0, Sigma0, Ytry, 
                         steps_bounds, stepsL, stepsU,
                         area, areaL, area_coef, s = n_coef) {

  area_sd <- sqrt(area_coef[3])
  # matrix with updated random effects
  Yout <- Y
  for(j in 1:J1) {
    Yold <- Y[, j]
    id <- sample(1:dim(Ytry)[3], size = 1)
    Ynew <- Ytry[, j, id]
    
    # extract unconstrained steps
    steps_old <- Yold[s]
    steps_new <- Ynew[s]
    
    # Transform steps in Yold and Ynew to the shared logit scale, to compute
    # hierarchical priors
    Yold_h <- Yold # duplicate to transform steps
    Ynew_h <- Ynew  
    Yold_h[s] <- logit_scaled(steps_old, stepsL, stepsU)
    Ynew_h[s] <- logit_scaled(steps_new, stepsL, stepsU)
    
    # (hierarchical priors, hpri)
    hpri_new <- 
      LaplacesDemon::dmvn(Ynew_h, mu[, j], Sigma, log = T) +
      logit_scaled_d_log(steps_new, stepsL, stepsU)
    hpri_old <- 
      LaplacesDemon::dmvn(Yold_h, mu[, j], Sigma, log = T) +
      logit_scaled_d_log(steps_old, stepsL, stepsU)
    
    # Transform steps in Yold and Ynew to the fire-wise logit scale, to compute
    # stage1 prior (vectors denoted with _pri1)
    Yold_pri1 <- Yold
    Ynew_pri1 <- Ynew
    
    Yold_pri1[s] <- logit_scaled(steps_old, steps_bounds[j, 1], steps_bounds[j, 2])
    Ynew_pri1[s] <- logit_scaled(steps_new, steps_bounds[j, 1], steps_bounds[j, 2])
    
    # stage1 priors (placeholders, ph)
    phpri_new <- 
      LaplacesDemon::dmvn(Ynew_pri1, mu0[, j], Sigma0[, , j], log = T) +
      logit_scaled_d_log(steps_new, steps_bounds[j, 1], steps_bounds[j, 2])
    
    phpri_old <- 
      LaplacesDemon::dmvn(Yold_pri1, mu0[, j], Sigma0[, , j], log = T) +
      logit_scaled_d_log(steps_old, steps_bounds[j, 1], steps_bounds[j, 2])
    
    # compute area mu from log-steps and area likelihood
    area_mu_new <- area_coef[1] + area_coef[2] * log(steps_new)
    area_mu_old <- area_coef[1] + area_coef[2] * log(steps_old)
    
    area_like_new <- dtruncnorm(area[j], a = areaL, mean = area_mu_new, 
                                sd = area_sd) |> log()
    area_like_old <- dtruncnorm(area[j], a = areaL, mean = area_mu_old, 
                                sd = area_sd) |> log()
    
    # log-posterior for raneff vector is
    
    #   hierarchical prior +  
    #   area likelihood    +
    #   stage1 prior * (-1)   
    
    # Log-posterior
    lp_new <- hpri_new + area_like_new - phpri_new
    lp_old <- hpri_old + area_like_old - phpri_old
      
    if(exp(lp_new - lp_old) > runif(1)) {
      Yout[, j] <- Ynew
    }
  }
  
  return(Yout)
}

# Update the logit-step parameter of fires for which spread was not simulated. They
# depend on the hierarchical prior and the log-area likelihood.

# steps: logit-steps current value.
# mu: mean of logit-steps, defined by the FWI.
# s: sd of logit-steps (not s^2).
# stepsL and stepsU: lower and upper limit for steps.
# area: log(area)
# areaL: lower limit for log(area) ~= log(10) = log(9.999999).
# area_coef: vector with intercept, slope and s2 for the area ~ steps regression.
# sd_jump: vector with proposal sd for each step parameter.

# returns a vector of updated (logit-)steps.
update_steps <- function(steps, mu, s, stepsL, stepsU, 
                         area, areaL, area_coef, sd_jump) {

  area_s <- sqrt(area_coef[3])
  steps_try <- rnorm(J2, steps, sd_jump)

  # log_posterior = area_likelihood + hierarchical_prior
  
  # get log-steps to compute area likelihood
  steps_log <- log(invlogit_scaled(steps, stepsL, stepsU))
  steps_log_try <- log(invlogit_scaled(steps_try, stepsL, stepsU))
  area_mu <- area_coef[1] + area_coef[2] * steps_log
  area_mu_try <- area_coef[1] + area_coef[2] * steps_log_try

  lp_new <-
    dtruncnorm(area, a = areaL, mean = area_mu_try, sd = area_s) |> log() +
    dnorm(steps_try, mu, s, log = TRUE)

  lp_old <-
    dtruncnorm(area, a = areaL, mean = area_mu, sd = area_s) |> log() +
    dnorm(steps, mu, s, log = TRUE)

  new_ids <- exp(lp_new - lp_old) > runif(J2)
  steps[new_ids] <- steps_try[new_ids]

  return(steps)
}

# Upper steps parameter (stepsU).

# The logit-steps for non-spread fires is a parameter independent of stepsU. 
# Hence, the logit-steps non-spread are not updated with U, but their log-step 
# depends on U, and that affects the area likelihood.
# Conversely, the logit-steps for spread fires is a derived quantity of U and 
# spread-steps.

# The log-posterior for stepsU is 
#   stepsU prior + 
#   area likelihood (only for non-spread fires) +
#   spread steps likelihood 

# As spread steps are in the original scale but their likelihood is at the logit,
# their density is dnorm(steps_logit, ...) + log(abs(jacobian{logit transform}))

# stepsU_logit: logit of the upper steps parameter at the current iteration, to be 
#   updated.
# Umin, Umax: lower and upper bounds for stepsU.
# steps1: steps of spread fires, at the constrained (simulator) scale. Used
#   to compute the logit-steps1, based on stepsU.
# mu_steps_logit1: logit-mean of the steps parameter for spread fires, to
#     compute the logit-steps likelihood.
# s: standard deviation of logit-steps.
# steps_logit2: steps for non-spread fires at the logit scale. stepsU changes 
#   their log(steps), and so the area likelihood.
# area: log(area) for non-spread fires.
# areaL: lower limit for log(area) ~= log(10) = log(9.999999).
# area_coef: vector with intercept, slope and s2 for the area ~ steps regression.
# sd_jump: proposal sd for stepsU_logit.
# stepsL: lower bound for steps, defaults to 2.

update_stepsU <- function(stepsU_logit, Umin, Umax, 
                          steps1, mu_steps_logit1, s,
                          steps_logit2, 
                          area, areaL, area_coef, 
                          sd_jump,
                          stepsL = 2) {
  
  stepsU_logit_try <- rnorm(1, stepsU_logit, sd_jump)
  stepsU <- invlogit_scaled(stepsU_logit, Umin, Umax)
  stepsU_try <- invlogit_scaled(stepsU_logit_try, Umin, Umax)
  
  # Compute steps_logit1
  steps_logit1 <- logit_scaled(steps1, stepsL, stepsU)
  steps_logit1_try <- logit_scaled(steps1, stepsL, stepsU_try)
  
  # Constrained steps2 
  steps2 <- invlogit_scaled(steps_logit2, stepsL, stepsU)
  steps2_try <- invlogit_scaled(steps_logit2, stepsL, stepsU_try)
  
  # area mu for non-spread fires
  area_s <- sqrt(area_coef[3])
  area_mu <- area_coef[1] + area_coef[2] * log(steps2)
  area_mu_try <- area_coef[1] + area_coef[2] * log(steps2_try)
  
  lp_new <-
    # area likelihood (non-spread)
    dtruncnorm(area, a = areaL, mean = area_mu_try, sd = area_s) |> log() |> sum() +
    # steps1 likelihood with jacobian adjustment
    dnorm(steps_logit1_try, mu_steps_logit1, s, log = TRUE) |> sum() +
    logit_scaled_d_log(steps1, stepsL, stepsU_try) |> sum() + 
    # stepsU_logit prior
    dlogis(stepsU_logit_try, log = T)
  
  lp_old <-
    # area likelihood (non-spread)
    dtruncnorm(area, a = areaL, mean = area_mu, sd = area_s) |> log() |> sum() +
    # steps1 likelihood with jacobian adjustment
    dnorm(steps_logit1, mu_steps_logit1, s, log = TRUE) |> sum() +
    logit_scaled_d_log(steps1, stepsL, stepsU) |> sum() + 
    # stepsU_logit prior
    dlogis(stepsU_logit, log = T)
  
  out <- ifelse(exp(lp_new - lp_old) > runif(1), 
                stepsU_logit_try, stepsU_logit)
  return(out)
}