# Packages ---------------------------------------------------------------

library(tidyverse)
library(rstan)         # sample with HMC to check
library(LaplacesDemon) # rinvwishart
library(invgamma)
library(MASS)          # mvrnorm, not so sensitive to positive-definiteness
library(trialr)        # rlkjcorr
library(truncnorm)
library(truncreg)      # tnorm regression

# Functions ---------------------------------------------------------------

# Functions to update parameters in the mcmc run.
# Spread hyperparameters are updated with gibbs for linear regression, including
#   each sigma2;
# The correlation matrix among random effects is updated using an inverse-wishart
#   prior, but only the correlation is used. This is done to avoide the correlation
#   between scale and correlations in the inverse-wishart distribution.
# The random effects are updated through MH, but using as proposal the samples
#   from stage1. Here the stage1 prior must be removed, and the jacobian for the
#   steps transform must be considered in the posterior.
# For the area ~ steps regression, with a trunctad-normal likelihood, MH updates
#   are performed.

# Perform the gibbs update for a linear regresion.
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

# function to transfrom log-steps to logit scale
exp_logit <- function(x, L, U) {
  exp_x <- exp(x)
  exp_unit <- (exp_x - L) / (U - L)
  return(qlogis(exp_unit))
}

# derivative of exp_logit, used to add the log-jacobian of exp_logit
# to the posterior density.
exp_logit_d <- function(x, L, U) {
  exp_x <- exp(x)
  d <- (exp_x * (L - U)) / ((exp_x - L) * (exp_x - U))
  return(d)
}

# update spread parameters.
# Y: matrix of random effects, with replicates in columns
# mu: matrix of expected values for random effects, with the same dimension
#   as Y.
# Sigma: vcov matrix for random effects, computed from scales and correlations.
# mu0: matrix with the same dimension as Y, with the prior mean from stage1 for
#   each random effect.
# Sigma0: array with the prior vcov for each random effect in each slice.
# Ytry: list of matrices with samples from stage1 for all random effects,
#   with an additional column indicating the weight of each sample.
#   with dimension [[n_fires]][n_samples, n_coef + 1]
# steps_bounds: matrix with steps bounds to convert from log to logit.
#   Fires in rows.
# s: position of steps parameter in the parameter vector
# area: log-area.
# areaL: lower bound for log-area.
# area_coef: b1, b2 and s2 for area, to compute likelihood conditional on steps.
update_ranef <- function(Y, mu, Sigma, mu0, Sigma0, Ytry, steps_bounds,
                         area, areaL, area_coef, s = n_coef) {

  area_sd <- sqrt(area_coef[3])
  # matrix with updated random effects
  Yout <- Y
  for(j in 1:J1) {
    Yold <- Y[, j]
    id <- sample(1:dim(Ytry)[3], size = 1)
    Ynew <- Ytry[, j, id]

    # transform the steps to logit, to compute the prior
    Yold_logit <- Yold
    Ynew_logit <- Ynew

    Yold_logit[s] <- exp_logit(Yold[s], steps_bounds[j, 1], steps_bounds[j, 2])
    Ynew_logit[s] <- exp_logit(Ynew[s], steps_bounds[j, 1], steps_bounds[j, 2])

    # log-posterior for raneff is
    #   hierarchical prior +
    #   stage1 prior * (-1) +
    #   log-jacobian for exp_logit(steps) transform +
    #   area likelihood

    area_mu_new <- area_coef[1] + area_coef[2] * Ynew[s]
    area_mu_old <- area_coef[1] + area_coef[2] * Yold[s]

    lp_new <-
      LaplacesDemon::dmvn(Ynew, mu[, j], Sigma, log = T) -
      LaplacesDemon::dmvn(Ynew_logit, mu0[, j], Sigma0[, , j], log = T) +
      exp_logit_d(Ynew[s], steps_bounds[j, 1], steps_bounds[j, 2]) |> abs() |> log() +
      dtruncnorm(area[j], a = areaL, mean = area_mu_new, sd = area_sd) |> log()

    lp_old <-
      LaplacesDemon::dmvn(Yold, mu[, j], Sigma, log = T) -
      LaplacesDemon::dmvn(Yold_logit, mu0[, j], Sigma0[, , j], log = T) +
      exp_logit_d(Yold[s], steps_bounds[j, 1], steps_bounds[j, 2]) |> abs() |> log() +
      dtruncnorm(area[j], a = areaL, mean = area_mu_old, sd = area_sd) |> log()

    if(exp(lp_new - lp_old) > runif(1)) {
      Yout[, j] <- Ynew
    }
  }

  return(Yout)
}

# Update the step parameter of fires for which spread was not simulated. They
# depend on the hierarchical prior and the log-area likelihood.

# steps: log-steps current value.
# mu: mean of log-steps, defined by the FWI.
# s: sd of log-steps (not s^2!).
# area: log(area)
# areaL: lower limit for log(area) = log(10).
# area_coef: vector with intercept, slope and s2 for the area ~ steps regression.
# sd_jump: vector with proposal sd for each step parameter.

# returns a vector of updated steps.
update_steps <- function(steps, mu, s, area, areaL, area_coef, sd_jump) {

  area_s <- sqrt(area_coef[3])
  steps_try <- rnorm(J2, steps, sd_jump)

  # log_posterior = area_likelihood + hierarchical_prior
  area_mu <- area_coef[1] + area_coef[2] * steps
  area_mu_try <- area_coef[1] + area_coef[2] * steps_try

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
