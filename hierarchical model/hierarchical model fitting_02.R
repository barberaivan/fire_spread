# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(bayesplot)
library(bayestestR) # hdi
library(posterior)
library(deeptime)
library(terra)
library(kdevine)
library(truncnorm)
library(mgcv)          # models to tune proposals
library(foreach)       # parallelization
library(doMC)          # doRNG had problems with cpp functions
theme_set(theme_bw())

# Functions ---------------------------------------------------------------

summarise <- function(x) {
  q <- quantile(x, probs = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
                method = 8)
  names(q) <- c("eti_lower_95", "eti_lower_90", "eti_lower_80", "median",
                "eti_upper_80", "eti_upper_90", "eti_upper_95")

  hdi_95 <- hdi(x, ci = 0.95)
  hdi_90 <- hdi(x, ci = 0.90)
  hdi_80 <- hdi(x, ci = 0.80)

  hdis <- c(
    "hdi_lower_95" = hdi_95$CI_low, "hdi_lower_90" = hdi_90$CI_low,
    "hdi_lower_80" = hdi_80$CI_low,
    "hdi_upper_80" = hdi_80$CI_high, "hdi_upper_90" = hdi_90$CI_high,
    "hdi_upper_95" = hdi_95$CI_high
  )

  res <- c("mean" = mean(x), hdis, q)

  return(res)
}

# Function to transform from original (constrained) support to the
# unconstrained one, using scaled logit and log (for steps)
unconstrain <- function(x, support) {
  xun <- x

  names_logit <- colnames(x)[colnames(x) != "steps"]
  for(j in names_logit) {
    xun[, j] <- qlogis((x[, j] - support[1, j]) / (support[2, j] - support[1, j]))
  }

  xun[, "steps"] <- log(x[, "steps"])

  return(xun)
}

# The inverse of unconstrain
constrain <- function(xun, support) {
  xc <- xun

  names_logit <- colnames(xun)[colnames(xun) != "steps"]
  for(j in names_logit) {
    xc[, j] <- plogis(xun[, j]) * (support[2, j] - support[1, j]) + support[1, j]
  }

  xc[, "steps"] <- exp(xun[, "steps"])

  return(xc)
}

# In the first stage, steps were assigned a flat prior, between 5 and S, with S
# varying among fires. As we aim to estimate the steps at the log scale, we
# need to know which was the prior distribution at the log scale, which is not
# a uniform, but an exponential-uniform distribution with pdf:

dexpunif <- function(x, l = 1, u = 10) { # for x \in [log(l), log(u)]
  d <- ifelse(x >= log(l) & x <= log(u),
              exp(x) / (u-l), 0)
  return(d)
}

expunif_lpdf <- function(x, l = 1, u = 10) { # for x \in [log(l), log(u)]
  ld <- ifelse(x >= log(l) & x <= log(u),
               x, -Inf)
  return(ld)
}

expunif_lpdf2 <- function(x, l = 1, u = 10) { # for x \in [log(l), log(u)]
  ld <- ifelse(x >= log(l-1e-12) & x <= log(u+1e-12),
               x, -Inf)
  return(ld)
}

# simulator of random effects, based on the ranef_kde list (a list with a
# kdevine density for each fire.)
ranef_rng <- function(ids) {
  x <- sapply(ids, function(i) {
    rkdevine(1, ranef_kde[[i]])
  })
  return(x)
}

# Inverse-Gamma densities
# parameterizations from Hooten and Hefley 2019:
digamma <- function(x, q = 1, r = 1000) {
  x ^ (-(q + 1)) * exp(-1 / r / x) / (r ^ q) / gamma(q)
}

digamma2 <- function(x, q = 1, r = 1000) {
  x ^ (-(q + 1)) * exp(-1 / (r * x)) / ((r ^ q) * gamma(q))
}

dinvgamma <- function(x, q = 1, r = 1000, log = F) {
  d <- 1 / (r ^ q * gamma(q)) * x ^ (-(q + 1)) * exp(-(1 / (r * x)))
  if(!log) {
    return(d)
  } else {
    return(log(d))
  }
}

# from wikipedia, parameterization based on alpha and beta, as the gamma
digamma4 <- function(x, a = q, b = 1 / r, q = 1, r = 1000) {
  # if(is.null(a)) a <- q
  # if(is.null(b)) b <- 1 / r
  b ^ a / gamma(a) * x ^ (-a-1) * exp(-(b / x))
}

# Simulate parameters from the empirical kde of each fire
rkde <- function(fnum, nsim = 1e5) {
  rkdevine(n = nsim, obj = ranef_kde[[fnum]])
}

# MCMC: sampler for the joint posterior

# nsim: number of iterations to run.
# sd_jump: list containing the proposal sd for each parameter that is updated
#   through m-h, using either normal or truncated normal in the case of sigma.
# start: list with initial values.
# samples: list with the output of a previous mcmc run, used to get starting
#   values.
# jump_factor: in the case of missing(sd_jump), the proposal sigmas are computed
#   from the MLE distribution, but multiplying them by the jump_factor.
# sd_jump_out: return the used sd_jump list?

# returns a list with 3 arrays of samples: fixef, ranef and steps_extra.
mcmc <- function(nsim = 50, sd_jump, start, samples, jump_factor = 4,
                 sd_jump_out = F) {

  ###### TESTO
  # nsim <- 50; jump_factor <- 4
  ###### ENDO TESTO

  #### Allocate memory to save samples,

  # fixef:       [parname, partype (a, b, s), iter]
  # ranef:       [parname, fire_id,           iter]
  # steps_extra: [1,       fire_id,           iter]

  # b param is NA for intercepts, and fixef include the area ~ steps
  # regression.

  fixef <- array(NA, dim = c(n_coef + 1, 3, nsim),
                 dimnames = list(
                   par_names = c(par_names, "area"),
                   par_class = c("a", "b", "s"),
                   iter = 1:nsim
                 ))

  ranef <- array(NA, dim = c(n_coef, nfires_spread, nsim),
                 dimnames = list(
                   par_names = par_names,
                   fire_id = fire_ids,
                   iter = 1:nsim
                 ))

  steps_extra <- array(NA, dim = c(nfires_nonspread, nsim),
                       dimnames = list(
                         fire_id = fire_ids_nonspread,
                         iter = 1:nsim
                       ))

  #### Define proposals sd (if missing)
  if(missing(sd_jump)) {

    # intercepts' alpha and sigma are updated through gibbs
    fixef_jump <- array(NA, dim = c(n_coef + 1 - 5, 3), # no intercepts
                        dimnames = list(
                          par_names = c(par_names[-(1:5)], "area"),
                          par_class = c("a", "b", "s")
                        ))

    fixef_jump["slope", ] <- apply(par_start$slope, 2, sd) * jump_factor
    fixef_jump["wind", ] <- apply(par_start$wind, 2, sd) * jump_factor
    fixef_jump["steps", ] <- apply(par_start$steps_fwi, 2, sd) * jump_factor
    fixef_jump["area", ] <- apply(par_start$area_steps, 2, sd) * jump_factor

    steps_extra_jump <- apply(par_start$steps_area$ranef, 2, sd) * jump_factor

    if(sd_jump_out) {
      sd_jump <- list(fixef = fixef_jump,
                      steps_extra = steps_extra_jump)
    }

  } else {
    fixef_jump <- sd_jump$fixef
    steps_extra_jump <- sd_jump$steps_extra
  }

  #### Define initial values (if missing)

  if(!missing(start)) {
    fixef_start <- start$fixef
    ranef_start <- start$ranef
    steps_extra_start <- start$steps_extra
  }

  # If there is no start, but there is a previous sample, set as start
  # the last value of the previous sample
  if(missing(start) & !missing(samples)) {
    nn <- dim(samples$fixef)[3]
    fixef_start <- samples$fixef[, , nn]
    ranef_start <- samples$ranef[, , nn]
    steps_extra_start <- samples$steps_extra[, nn]
  }

  # if no start nor samples are provide, make random beginning
  if(missing(start) & missing(samples)) {

    is <- sample(1:5000, size = 1)

    fixef_start <- rbind(
      c(par_start$wet[is, 1], NA, par_start$wet[is, 2]),
      c(par_start$subalpine[is, 1], NA, par_start$subalpine[is, 2]),
      c(par_start$dry[is, 1], NA, par_start$dry[is, 2]),
      c(par_start$shrubland[is, 1], NA, par_start$shrubland[is, 2]),
      c(par_start$grassland[is, 1], NA, par_start$grassland[is, 2]),

      par_start$slope[is, ],
      par_start$wind[is, ],
      par_start$steps_fwi[is, ],
      par_start$area_steps[is, ]
    )
    colnames(fixef_start) <- c("a", "b", "s")
    rownames(fixef_start) <- c(par_names, "area")

    # spread random effects
    ranef_start <- ranef_proposal[, , is]
    steps_extra_start <- par_start$steps_area$ranef[is, ]
  }

  #### Initialize chain
  fixef[, , 1] <- fixef_start
  ranef[, , 1] <- ranef_start
  steps_extra[, 1] <- steps_extra_start

  #### Transient matrices to store intermediate results
  mu_fitted_keep <- array(NA, dim = c(nfires_spread, 2),
                          dimnames = list(
                            fire_id = fire_ids,
                            param = c("slope", "wind")
                          ))

  # extract a few constant priors
  q <- fprior[1, "s", "mu_q"]
  r <- fprior[1, "s", "s_r"]

  # MCMC loop -------------------------------------------------------------
  for(k in 2:nsim) {
    # print(k)

    ### Fixed effects:
    # intercepts (a, s) updated through gibbs;
    # slope, wind, steps, use mh for a and b, gibbs for s;
    # area updated with mh because it uses truncated-normal likelihood.


    # Intercepts ----------------------------------------------------------
    for(v in 1:5) {
      # sample s from the conditional distribution
      q_tmp <- J1 / 2 + q
      r_tmp <- 1 / (sum((ranef[v, , k-1] - fixef[v, "a", k-1]) ^ 2) / 2 + 1 / r)
      s2 <- 1 / rgamma(1, q_tmp, 1 / r_tmp) # r = 1/beta = scale
      fixef[v, "s", k] <- sqrt(s2)

      # sample mu (alpha) from the conditional distribution
      mu0 <- fprior[v, "a", "mu_q"] # prior mu
      s0 <- fprior[v, "a", "s_r"]  # prior sigma
      tmp_var <- 1 / (J1 / s2 + 1 / s0^2)
      tmp_mn <- tmp_var * (sum(ranef[v, , k-1]) / s2 + mu0 / s0^2)
      fixef[v, "a", k] <- rnorm(1, tmp_mn, sqrt(tmp_var))
    }


    # Slope and wind -----------------------------------------------------
    for(v in c("slope", "wind")) {
      # sample s from the conditional distribution
      # (compute mu_fitted, from the fwi, to get error terms)
      mu_fitted <- fixef[v, "a", k-1] + fixef[v, "b", k-1] * fwi_sub
      q_tmp <- J1 / 2 + q
      r_tmp <- 1 / (sum((ranef[v, , k-1] - mu_fitted) ^ 2) / 2 + 1 / r)
      fixef[v, "s", k] <- sqrt(1 / rgamma(1, q_tmp, 1 / r_tmp)) # r = 1/beta = scale

      # a and b using metropolis
      a_try <- rnorm(1, fixef[v, "a", k-1], fixef_jump[v, "a"])
      mu_try <- a_try + fixef[v, "b", k-1] * fwi_sub
      lp_new <-
        sum(dnorm(ranef[v, , k-1], mu_try, fixef[v, "s", k], log = T)) + # likelihood of random effects
        dnorm(a_try,
              fprior[v, "a", "mu_q"], fprior[v, "a", "s_r"], log = T) # prior
      lp_old <-
        sum(dnorm(ranef[v, , k-1], mu_fitted, fixef[v, "s", k], log = T)) + # likelihood of random effects
        dnorm(fixef[v, "a", k-1],
              fprior[v, "a", "mu_q"], fprior[v, "a", "s_r"], log = T) # prior
      fixef[v, "a", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                                 a_try, fixef[v, "a", k-1])
      # update mu with updated a
      mu_fitted <- fixef[v, "a", k] + fixef[v, "b", k-1] * fwi_sub

      b_try <- rnorm(1, fixef[v, "b", k-1], fixef_jump[v, "b"])
      mu_try <- fixef[v, "a", k] + b_try * fwi_sub
      lp_new <-
        sum(dnorm(ranef[v, , k-1], mu_try, fixef[v, "s", k], log = T)) + # likelihood of random effects
        dnorm(b_try,
              fprior[v, "b", "mu_q"], fprior[v, "b", "s_r"], log = T) # prior
      lp_old <-
        sum(dnorm(ranef[v, , k-1], mu_fitted, fixef[v, "s", k], log = T)) + # likelihood of random effects
        dnorm(fixef[v, "b", k-1],
              fprior[v, "b", "mu_q"], fprior[v, "b", "s_r"], log = T) # prior
      fixef[v, "b", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                                 b_try, fixef[v, "b", k-1])
      # update mu with updated b, and save (will be used to update random effects)
      mu_fitted_keep[, v] <- fixef[v, "a", k] + fixef[v, "b", k] * fwi_sub
    }


    # Area parameters (~ steps) -------------------------------------------
    # As the likelihood is truncated-normal, the conditional distribution for
    # area_s (sigma) is unknown. Use M-H, but keep the inv-gamma prior for
    # consistency.
    v <- "area" # ease subsetting

    steps_all <- c(ranef["steps", , k-1], steps_extra[, k-1])
    mu_fitted <- fixef[v, "a", k-1] + fixef[v, "b", k-1] * steps_all
    area_s_try <- rtruncnorm(1, a = 0, mean = fixef[v, "s", k-1],
                             sd = fixef_jump[v, "s"])

    lp_new <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_fitted, sd = area_s_try) |> log()) +     # likelihood
      dinvgamma(area_s_try, log = T,
                q = fprior[v, "s", "mu_q"],
                r = fprior[v, "s", "s_r"]) +                           # prior
      dtruncnorm(fixef[v, "s", k-1], a = 0,
                 mean = area_s_try, sd = fixef_jump[v, "s"]) |> log()   # jump
    # this terms take into account the asymmetric proposal (truncnorm)

    lp_old <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_fitted, sd = fixef[v, "s", k-1]) |> log()) +   # likelihood
      dinvgamma(fixef[v, "s", k-1], log = T,
                q = fprior[v, "s", "mu_q"],
                r = fprior[v, "s", "s_r"]) +                                 # prior
      dtruncnorm(area_s_try, a = 0,
                 mean = fixef[v, "s", k-1], sd = fixef_jump[v, "s"]) |> log() # jump

    fixef[v, "s", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               area_s_try, fixef[v, "s", k-1])

    # a and b using metropolis
    a_try <- rnorm(1, fixef[v, "a", k-1], fixef_jump[v, "a"])
    mu_try <- a_try + fixef[v, "b", k-1] * steps_all
    lp_new <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_try, sd = fixef[v, "s", k]) |> log()) +      # likelihood
      dnorm(a_try,
            fprior[v, "a", "mu_q"], fprior[v, "a", "s_r"], log = T) # prior
    lp_old <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_fitted, sd = fixef[v, "s", k]) |> log()) +   # likelihood
      dnorm(fixef[v, "a", k-1],
            fprior[v, "a", "mu_q"], fprior[v, "a", "s_r"], log = T) # prior
    fixef[v, "a", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               a_try, fixef[v, "a", k-1])
    # update mu with updated a
    mu_fitted <- fixef[v, "a", k] + fixef[v, "b", k-1] * steps_all

    b_try <- rnorm(1, fixef[v, "b", k-1], fixef_jump[v, "b"])
    mu_try <- fixef[v, "a", k] + b_try * steps_all
    lp_new <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_try, sd = fixef[v, "s", k]) |> log()) +      # likelihood
      dnorm(b_try,
            fprior[v, "b", "mu_q"], fprior[v, "b", "s_r"], log = T) # prior
    lp_old <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_fitted, sd = fixef[v, "s", k]) |> log()) +   # likelihood
      dnorm(fixef[v, "b", k-1],
            fprior[v, "b", "mu_q"], fprior[v, "b", "s_r"], log = T) # prior
    fixef[v, "b", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               b_try, fixef[v, "b", k-1])

    # compute mu to compute the likelihood when updating steps
    area_mu <- fixef[v, "a", k] + fixef[v, "b", k] * steps_all


    # Steps parameters (~fwi) --------------------------------------------
    v <- "steps" # ease subsetting
    mu_fitted <- fixef[v, "a", k-1] + fixef[v, "b", k-1] * fwi_all

    # Update sigma with gibbs
    q_tmp <- J / 2 + q
    r_tmp <- 1 / (sum((steps_all - mu_fitted) ^ 2) / 2 + 1 / r) # prior is inv-gamma
    fixef[v, "s", k] <- sqrt(1 / rgamma(1, q_tmp, 1 / r_tmp)) # r = 1/b = scale

    # a and b using metropolis
    a_try <- rnorm(1, fixef[v, "a", k-1], fixef_jump[v, "a"])
    mu_try <- a_try + fixef[v, "b", k-1] * fwi_all
    lp_new <-
      sum(dnorm(steps_all, mu_try, fixef[v, "s", k], log = T)) +    # likelihood for ranef
      dnorm(a_try,
            fprior[v, "a", "mu_q"], fprior[v, "a", "s_r"], log = T) # prior
    lp_old <-
      sum(dnorm(steps_all, mu_fitted, fixef[v, "s", k], log = T)) + # likelihood for ranef
      dnorm(fixef[v, "a", k-1],
            fprior[v, "a", "mu_q"], fprior[v, "a", "s_r"], log = T) # prior
    fixef[v, "a", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               a_try, fixef[v, "a", k-1])
    # update mu with updated a
    mu_fitted <- fixef[v, "a", k] + fixef[v, "b", k-1] * fwi_all

    b_try <- rnorm(1, fixef[v, "b", k-1], fixef_jump[v, "b"])
    mu_try <- fixef[v, "a", k] + b_try * fwi_all
    lp_new <-
      sum(dnorm(steps_all, mu_try, fixef[v, "s", k], log = T)) +      # likelihood
      dnorm(b_try,
            fprior[v, "b", "mu_q"], fprior[v, "b", "s_r"], log = T) # prior
    lp_old <-
      sum(dnorm(steps_all, mu_fitted, fixef[v, "s", k], log = T)) +   # likelihood
      dnorm(fixef[v, "b", k-1],
            fprior[v, "b", "mu_q"], fprior[v, "b", "s_r"], log = T) # prior
    fixef[v, "b", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               b_try, fixef[v, "b", k-1])

    # keep steps mu to update random effects
    steps_fitted <- fixef[v, "a", k] + fixef[v, "b", k] * fwi_all


    # Random effects (Lunn method) ---------------------------------------
    ii <- sample(1:ns1, 1)
    ranef_try <- ranef_proposal[, , ii]

    # Log-posterior of random effects is
    # area loglik + hierarchical prior - stage1 prior

    # area likelihood
    area_mu_try <- fixef["area", "a", k] +
      fixef["area", "b", k] * ranef_try["steps", ]
    area_like_try <- dtruncnorm(area_all[ids1], a = areaL,
                                mean = area_mu_try,
                                sd = fixef["area", "s", k]) |> log()

    # hierarchical prior
    priorh_try <- ranef_try # just initialize
    for(v in 1:5) {
      priorh_try[v, ] <- dnorm(
        ranef_try[v, ], fixef[v, "a", k], fixef[v, "s", k], log = T
      )
    }
    for(v in c("slope", "wind")) {
      priorh_try[v, ] <- dnorm(
        ranef_try[v, ], mu_fitted_keep[, v], fixef[v, "s", k], log = T
      )
    }
    priorh_try["steps", ] <- dnorm(
      ranef_try["steps", ], steps_fitted[ids1], fixef["steps", "s", k], log = T
    )
    priorh_try <- colSums(priorh_try)

    # stage1 prior
    prior1_try <- ranef_try # just initialize
    prior1_try["steps", ] <- expunif_lpdf(ranef_try["steps", ],
                                          l = steps_bounds[, 1],
                                          u = steps_bounds[, 2])
    prior1_try[1:7, ] <- dlogis(ranef_try[1:7, ], log = T)
    prior1_try <- colSums(prior1_try)

    # sum like and priors
    lp_new <- area_like_try + priorh_try - prior1_try

    # log posterior for previous ranef:

    # area likelihood
    area_like_old <- dtruncnorm(area_all[ids1], a = areaL,
                                mean = area_mu[ids1],
                                sd = fixef["area", "s", k]) |> log()

    # hierarchical prior
    priorh_old <- ranef[, , k-1] # just initialize
    for(v in 1:5) {
      priorh_old[v, ] <- dnorm(
        ranef[v, , k-1], fixef[v, "a", k], fixef[v, "s", k], log = T
      )
    }
    for(v in c("slope", "wind")) {
      priorh_old[v, ] <- dnorm(
        ranef[v, , k-1], mu_fitted_keep[, v], fixef[v, "s", k], log = T
      )
    }
    priorh_old["steps", ] <- dnorm(
      ranef["steps", , k-1], steps_fitted[ids1], fixef["steps", "s", k], log = T
    )
    priorh_old <- colSums(priorh_old)

    # stage1 prior
    prior1_old <- ranef[, , k-1] # just initialize
    prior1_old["steps", ] <- expunif_lpdf(ranef["steps", , k-1],
                                          l = steps_bounds[, 1],
                                          u = steps_bounds[, 2])
    prior1_old[1:7, ] <- dlogis(ranef[1:7, , k-1], log = T)
    prior1_old <- colSums(prior1_old)

    # sum like and priors
    lp_old <- area_like_old + priorh_old - prior1_old

    # MH random effects
    idx_keep <- exp(lp_new - lp_old) > runif(J1)
    ranef[, , k] <- ranef[, , k-1]                # duplicate old values
    ranef[, idx_keep, k] <- ranef_try[, idx_keep] # use accepted


    # Random effects (steps extra) ----------------------------------------
    steps_try <- rnorm(J2, steps_extra[, k-1], steps_extra_jump)

    # log_posterior = area_likelihood + hierarchical_prior
    area_mu_try2 <- fixef["area", "a", k] + fixef["area", "b", k] * steps_try

    lp_new <-
      dtruncnorm(area_all[ids2], a = areaL,
                 mean = area_mu_try2, sd = fixef["area", "s", k]) |> log() +
      dnorm(steps_try,
            steps_fitted[ids2], fixef["steps", "s", k], log = TRUE)

    lp_old <-
      dtruncnorm(area_all[ids2], a = areaL,
                 mean = area_mu[ids2], sd = fixef["area", "s", k]) |> log() +
      dnorm(steps_extra[, k-1],
            steps_fitted[ids2], fixef["steps", "s", k], log = TRUE)

    keep_idx <- exp(lp_new - lp_old) > runif(J2)
    steps_extra[keep_idx, k] <- steps_try[keep_idx]
    steps_extra[!keep_idx, k] <- steps_extra[!keep_idx, k-1]
  }

  #### Merge samples
  out <- list(
    fixef = fixef,
    ranef = ranef,
    steps_extra = steps_extra
  )

  if(sd_jump_out) out$sd_jump <- sd_jump

  return(out)
}

# function to count succesive changes in a vector
count_changes <- function(x) sum(abs(diff(x)) > 0)

# Compute the acceptance rate for the m-h updated parameters.
acceptance <- function(samples) {
  nn <- dim(samples$fixef)[3]
  nt <- nn - 1 # transitions
  pfixef <- apply(samples$fixef[6:9, , ], 1:2, count_changes) / nt
  psteps <- unname(apply(samples$steps_extra, 1, count_changes)) / nt
  out <- list(fixef = pfixef, steps_extra = psteps)
  return(out)
}

thin <- function(samples, rate = 100) {
  nn <- dim(samples$fixef)[3]
  if(nn < (rate * 10)) return(samples)

  ids <- rev(seq(nn, 1, by = -rate))

  out <- list(
    fixef = samples$fixef[, , ids],
    ranef = samples$ranef[, , ids],
    steps_extra = samples$steps_extra[, ids]
  )
  return(out)
}

# climatic data ------------------------------------------------------------

fwi_data <- read.csv(file.path("data", "climatic_data_by_fire_fwi-fortnight-cumulative.csv"))

# fire polygons (to get area) ---------------------------------------------

ff <- vect("data/patagonian_fires_spread.shp")
ff$area_ha <- expanse(ff) / 1e4 # turn m2 to ha

# A few constants ---------------------------------------------------------

par_names <- c("wet", "subalpine", "dry", "shrubland", "grassland",
               "slope", "wind", "steps")
par_names_all <- c(par_names, "area")

n_coef <- length(par_names)
# load file with constants to standardize
ndvi_params <- readRDS(file.path("data", "NDVI_regional_data",
                                 "ndvi_optim_and_proportion.rds"))
slope_sd <- ndvi_params$slope_term_sd

# support for parameters
n_veg <- 5
n_terrain <- 2

ext_alpha <- 50
ext_beta <- 30

params_lower <- c(rep(-ext_alpha, n_veg), rep(0, n_terrain))
params_upper <- c(rep(ext_alpha, n_veg), ext_beta / slope_sd, ext_beta)

support <- rbind(params_lower, params_upper)
colnames(support) <- names(params_lower) <- names(params_upper) <- par_names[-n_coef]
support_width <- apply(support, 2, diff)

# Load and prepare data ---------------------------------------------------

# dir to load files
target_dir <- file.path("files", "overlaps")

samples_files <- list.files(target_dir, pattern = "-posterior_samples.rds")

# import a list with matrix-samples
samples_list <- lapply(samples_files, function(s) {
  rrr <- readRDS(file.path(target_dir, s))
  return(rrr)
})

fire_ids <- lapply(samples_list, function(x) attr(x, "fire_id")) |> unlist()
names(samples_list) <- fire_ids

samples_df <- do.call("rbind", lapply(samples_list, function(x) {
  fire_id <- attr(x, "fire_id")
  r <- cbind(
    data.frame(fire_id = fire_id),
    as.data.frame(x)
  )
  return(r)
}))

# make list of samples in the unconstrained scale
samples_list_unc <- lapply(samples_list, function(x) {
  x2 <- unconstrain(x, support)
})

# merge with FWI data
# two fires were split, but they have the same FWI.
fires_data_0 <- fwi_data[!(fwi_data$fire_id %in% c("2011_19", "2015_47")), ]
fires_data_1 <- fwi_data[fwi_data$fire_id %in% c("2011_19", "2015_47"), ]
fires_data_2 <- fires_data_1[c(1, 1, 2, 2), ]

fires_data_2$fire_id[fires_data_2$fire_id == "2011_19"] <-
  fire_ids[grep("2011_19", fire_ids)]

fires_data_2$fire_id[fires_data_2$fire_id == "2015_47"] <-
  fire_ids[grep("2015_47", fire_ids)]

# bring area of separated fires
for(i in 1:nrow(fires_data_2)) {
  fires_data_2$area_ha[i] <- ff$area_ha[ff$fire_id == fires_data_2$fire_id[i]]
}

# put together all fires
fires_data <- rbind(fires_data_0, fires_data_2)
rownames(fires_data) <- NULL

# add log area and scaled fwi
fires_data$area_ha_log <- log(fires_data$area_ha)

fwi_mean <- mean(fires_data$fwi_expquad_fortnight)
fwi_sd <- sd(fires_data$fwi_expquad_fortnight)

fires_data$fwi <- (fires_data$fwi_expquad_fortnight - fwi_mean) / fwi_sd
rownames(fires_data) <- fires_data$fire_id

fires_data_spread <- fires_data[fires_data$fire_id %in% fire_ids, ]
fires_data_spread <- fires_data_spread[fire_ids, ]

fires_data_nonspread <- fires_data[!(fires_data$fire_id %in% fire_ids), ]

fire_ids_nonspread <- fires_data_nonspread$fire_id
nfires_nonspread <- length(fire_ids_nonspread)
nfires_spread <- length(fire_ids)

# Get bounds for steps parameters
steps_bounds <- do.call("rbind", lapply(samples_list, function(x) {
  bb <- attr(x, "support")
  return(bb[, "steps"])
}))

# lower bound for fire area (10 ha)
areaL <- log(10)

# N of random effects
J1 <- nfires_spread
J2 <- nfires_nonspread
J <- J1 + J2
ids1 <- 1:J1
ids2 <- (J1+1):J

fwi_sub <- fires_data_spread$fwi
fwi_all <- c(fires_data_spread$fwi, fires_data_nonspread$fwi)
area_all <- c(fires_data_spread$area_ha_log, fires_data_nonspread$area_ha_log)

# array with priors.
# for normal priors, par1 = mu, par2 = sigma.
# for inv-gamma priors, par1 = q, par2 = r

# Brute estimates ---------------------------------------------------------

# using samples from the first-step abc, get a point estimate of fixef
# effects. Those are:

# intercepts, mean and sd
# alpha, beta and sd for log_steps ~ log_area, with
#   predicted steps for non-spread fires (to use as starting values)
# alpha, beta and sd for (slope, wind, steps) ~ fwi

par_start <- vector("list", n_coef + 2)
names(par_start) <- c(par_names[-n_coef],
                      "steps_area", "steps_fwi", "area_steps")
nsim <- 5000

# make matrices for intercepts
for(v in par_names[1:5]) {
  mm <- matrix(NA, nsim, 2)
  colnames(mm) <- c("alpha", "sigma")
  par_start[[v]] <- mm
}

# make matrices for slope, wind and steps
for(v in c("slope", "wind", "steps_fwi", "area_steps")) {
  mm <- matrix(NA, nsim, 3)
  colnames(mm) <- c("alpha", "beta", "sigma")
  par_start[[v]] <- mm
}

# make matrices for steps ~ area
for(v in "steps_area") {
  mm1 <- matrix(NA, nsim, 3)
  colnames(mm1) <- c("alpha", "beta", "sigma")

  mm2 <- matrix(NA, nsim, nfires_nonspread)
  colnames(mm2) <- fire_ids_nonspread

  par_start[[v]] <- list(fixef = mm1, ranef = mm2)
}

# Turn samples into array
samples_temp <- lapply(samples_list_unc, function(x) x[1:nsim, ])
samples_arr <- abind::abind(samples_temp, along = 3)

# # Compute estimates
# for(i in 1:nsim) {
#   # i = 1
#   print(i)
#
#   # intercepts
#   for(v in par_names[1:5]) {
#     par_start[[v]][i, "alpha"] <- mean(samples_arr[i, v, ])
#     par_start[[v]][i, "sigma"] <- sd(samples_arr[i, v, ])
#   }
#
#   # slope, wind
#   for(v in c("slope", "wind")) {
#     parvals <- samples_arr[i, v, ]
#     mod <- lm(parvals ~ fwi, data = fires_data_spread)
#     par_start[[v]][i, ] <- c(coef(mod), sigma(mod))
#   }
#
#   # steps ~ area model
#   # used to provide initial values for the non-spread steps parameters
#   ss <- samples_arr[i, "steps", ]
#   mod <- lm(ss ~ area_ha_log, data = fires_data_spread)
#   par_start[["steps_area"]][["fixef"]][i, ] <- c(coef(mod), sigma(mod)) # will not be used
#   # simulate steps in non-spread fires (not predict!)
#   # before I was predicting them, not simulating
#   #
#   mu <- predict(mod, fires_data_nonspread)
#   sigma <- sigma(mod)
#   steps_sim <- rnorm(J2, mu, sigma)
#   par_start[["steps_area"]][["ranef"]][i, ] <- steps_sim
#
#   # estimate steps ~ fwi and area ~ steps using all fires
#   ss_full <- c(ss, steps_sim)
#   # steps ~ fwi
#   mod <- lm(ss_full ~ fwi_all)
#   par_start[["steps_fwi"]][i, ] <- c(coef(mod), sigma(mod))
#   # area ~ steps
#   mod <- lm(ss_full ~ area_all)
#   par_start[["area_steps"]][i, ] <- c(coef(mod), sigma(mod))
# }
# saveRDS(par_start, file.path("files", "hierarchical_model", "par_start.rds")) # predict steps
# par_start <- readRDS(file.path("files", "hierarchical_model", "par_start.rds"))
# saveRDS(par_start, file.path("files", "hierarchical_model", "par_start2.rds")) # simulate steps
par_start <- readRDS(file.path("files", "hierarchical_model", "par_start2.rds"))


# Priors for hyperparameters ----------------------------------------------

# fprior for fixed-effects prior
fprior <- array(NA, dim = c(n_coef + 1, 3, 2),
                dimnames = list(
                  par_names = c(par_names, "area"),
                  par_class = c("a", "b", "s"),
                  prior_par = c("mu_q", "s_r")
                ))

# a and b for logit-bounded parameters:
#   normal(0, 10)
fprior[, c("a", "b"), "mu_q"] <- 0
fprior[, c("a", "b"), "s_r"] <- 10

# but for steps and area (log link), place the mean of the intercept
# close to the MLE.
fprior[8:9, "a", "mu_q"] <- c(mean(par_start$steps_fwi[, "alpha"]),
                              mean(par_start$area_steps[, "alpha"]))

# s for all parameters has inv-gamma prior, for conjugacy.
#   inv-gamma(x, q = 1, r = 1000)
fprior[, "s", "mu_q"] <- 1
fprior[, "s", "s_r"] <- 1000

# Compute kdes for random effects (spread parameters) ----------------------

# ranef_kde <- lapply(1:nfires_spread, function(i) {
#   print(i)
#   sup_steps_log <- attr(samples_list[[i]], "support")[, "steps"] |> log
#   xx <- samples_temp[[i]] # unconstrained
#   kde <- kdevine(
#     xx,
#     xmin = c(rep(-Inf, n_coef-1), sup_steps_log[1]),
#     xmac = c(rep(Inf, n_coef-1), sup_steps_log[2]),
#     copula.type = "kde",
#     cores = 6
#   )
#   return(kde)
# })
# saveRDS(ranef_kde, file.path("files", "hierarchical_model", "ranef_kde.rds"))
# ranef_kde <- readRDS(file.path("files", "hierarchical_model", "ranef_kde.rds"))

# Expand random effects samples -------------------------------------------

# The ABC-sampling took around ~ 5000 samples from each fire. To expand those
# samples up to a larger number, we will use the fitted KDEs. But as simulation
# with kdevine is expensive, we will simulate just once, and then, resample from
# the larger number of samples, as in the classic Lunn method.

# registerDoMC(16)
#
# # 20 batches of 5000 samples.
# B <- 20
# bsize <- 5000
# slist <- vector("list", B)
#
# for(b in 1:B) {
#
#   ids_list <- sapply(1:J1, function(x) x)
#   sims <- foreach(ii = ids_list) %dopar% {
#     rkde(ii, bsize)
#   }
#
#   sims <- abind::abind(sims, along = 3)
#   fname <- paste("samples_batch_", b, ".rds", sep = "")
#   ff <- file.path("files", "hierarchical_model", fname)
#   saveRDS(sims, ff)
#
#   slist[[b]] <- sims
# }
#
# outout <- abind::abind(slist, along = 1)
# saveRDS(outout, file.path("files", "hierarchical_model",
#                           "samples_stage1_expanded.rds"))

a1 <- readRDS(file.path("files", "hierarchical_model",
                        "samples_stage1_expanded.rds"))
dimnames(a1) <- list(
  iter = 1:nrow(a1),
  param = par_names,
  fire = fire_ids
)

l1 <- lapply(1:J1, function(i) {
  # i = 1
  mm <- a1[, , i]
  mm_prior <- mm
  for(c in 1:7) {
    mm_prior[, c] <- dlogis(mm[, c], log = T)
  }
  mm_prior[, 8] <- expunif_lpdf(mm[, 8],
                                l = steps_bounds[i, 1],
                                u = steps_bounds[i, 2])
  lprior <- rowSums(mm_prior)
  keep <- is.finite(lprior) & !(is.na(lprior))
  return(mm[keep, ])
})

availables <- lapply(l1, function(x) nrow(x)) |> unlist()
ns1 <- 90000 # min(availables) # 96385

l2 <- lapply(l1, function(x) x[1:ns1, ])
dim(l2[[1]])

# [parname, fire_id,           iter]
ranef_proposal <- abind::abind(l2, along = 3)
ranef_proposal <- aperm(ranef_proposal, c(2, 3, 1))

# Try a few runs ----------------------------------------------------------

nsim <- 1000
run1 <- mcmc(nsim = nsim)

rri <- range(as.vector(run1$fixef[1:5, c(1, 3), ]))
par(mfrow = c(3, 2))
for(v in 1:5) {
  matplot(1:nsim, y = t(run1$fixef[v, c(1, 3), ]), type = "l",
          ylab = NA, xlab = NA,
          main = par_names[v], lty = 1, col = c(1, 3),
          ylim = rri)

}
par(mfrow = c(1, 1))

rri <- range(as.vector(run1$fixef[6:9, , ]))
par(mfrow = c(2, 2))
for(v in 6:9) {
  matplot(1:nsim, y = t(run1$fixef[v, , ]), type = "l",
          ylab = NA, xlab = NA, ylim = rri,
          main = par_names_all[v], lty = 1, col = c(1, 2, 3))
}
par(mfrow = c(1, 1))

# acceptance rate (note how gibbs makes beautiful rates)
apply(run1$fixef, 1:2, function(x) sum(abs(diff(x)) > 0)) / (nsim-1)
unname(apply(run1$steps_extra, 1, function(x) sum(abs(diff(x)) > 0)) / (nsim-1))


acf(run1$fixef["wet", "a", ], lag.max = 200)
acf(run1$fixef["wet", "s", ], lag.max = 200)

acf(run1$fixef["steps", "a", ], lag.max = 200)
acf(run1$fixef["steps", "s", ], lag.max = 200)
acf(run1$fixef["steps", "s", ], lag.max = 200)


# MCMC adaptation ---------------------------------------------------------

# # run a long chain to get good starting points
# nsim <- 1e6
# run0 <- mcmc(nsim = nsim, sd_jump_out = T)
# a0 <- acceptance(run0)
# a0$fixef
#
# run0_thin <- thin(run0, 100)
# saveRDS(run0_thin, file.path("files", "hierarchical_model", "run0_thin.rds"))
run0_thin <- readRDS(file.path("files", "hierarchical_model", "run0_thin.rds"))
nsim <- dim(run0_thin$fixef)[3]
# # Visualize... steady state?
# rri <- range(as.vector(run0_thin$fixef[6:9, , ]))
# par(mfrow = c(2, 2))
# for(v in 6:9) {
#   matplot(1:nsim, y = t(run0_thin$fixef[v, , ]), type = "l",
#           ylab = NA, xlab = NA, ylim = rri,
#           main = par_names_all[v], lty = 1, col = c(1, 2, 3))
# }
# par(mfrow = c(1, 1))
# acf(run0_thin$fixef["slope", "a", ], lag.max = 1000)
## OK

# Initial proposal sigma:
sd_jump1 <- list(
  fixef = sqrt(apply(run0_thin$fixef[6:9, , ], 1:2, sd) ^ 2 * 2),
  steps_extra = sqrt(apply(run0_thin$steps_extra, 1, sd) ^ 2 * 2)
)
# posterior sd after a long run, just for curiosity
sd_post <- list(
  fixef = apply(run0_thin$fixef[6:9, , ], 1:2, sd),
  steps_extra = apply(run0_thin$steps_extra, 1, sd)
)

sns <- 5000 # iter by step

# first run, to get acceptance
run1 <- mcmc(nsim = ns, sd_jump1, samples = run0_thin, sd_jump_out = T)
a1 <- acceptance(run1)

# acceptance trial steps
K <- 50

# make space for acceptance and sigma
accept_track <- list(
  fixef = run1$fixef[6:9, , 1:K],
  steps_extra = run1$steps_extra[, 1:K]
)
sigma_track <- list(
  fixef = run1$fixef[6:9, , 1:K],
  steps_extra = run1$steps_extra[, 1:K]
)
accept_track$fixef[, , ] <- NA
accept_track$steps_extra[, ] <- NA
sigma_track$fixef[, , ] <- NA
sigma_track$steps_extra[, ] <- NA
# just make space

# fill sigma and accept in the first run
accept_track$fixef[, , 1] <- a1$fixef
accept_track$steps_extra[, 1] <- a1$steps_extra
sigma_track$fixef[, , 1] <- sd_jump1$fixef
sigma_track$steps_extra[, 1] <- sd_jump1$steps_extra

# Fill 4 sigma values below and above the first, to fit the first regression
# using 5 data points
factors <- seq(0.1, 2, by = 0.2)
lf <- length(factors)
for(k in 2:(lf+1)) {
  print(k)
  sigma_track$fixef[, , k] <- sd_jump1$fixef * factors[k-1]
  sigma_track$steps_extra[, k] <- sd_jump1$steps_extra * factors[k-1]

  # run MCMC
  sss <- list(fixef = sigma_track$fixef[, , k],
              steps_extra = sigma_track$steps_extra[, k])
  run_k <- mcmc(ns, sd_jump = sss, samples = run0_thin)

  # compute and store acceptance
  a_k <- acceptance(run_k)
  accept_track$fixef[, , k] <- a_k$fixef
  accept_track$steps_extra[, k] <- a_k$steps_extra
}

# iterative fitting
# for(k in (lf+2):K) {
for(k in 20:K) {
  # k = 12
  print(k)
  # use previous runs to fit a regression of sigma ~ accept, and choose the
  # predicted sigma for accept = 0.44

  # fixef
  for(j in 1:3) {   # par type
    for(i in 1:4) { # spread parameter

      # i = 4; j = 1 # TEST

      if(!(j == 3 & i < 4)) {
        # get previous data
        dd <- data.frame(
          aa = accept_track$fixef[i, j, 1:(k-1)],
          ss = sigma_track$fixef[i, j, 1:(k-1)]
        )

        # remove sigma too close to zero
        dd <- dd[dd$ss >= 1e-4, ]
        dd$lss = log(dd$ss)

        # fit regression
        if(k < 20) {
          mm <- lm(lss ~ aa + I(aa ^ 2), data = dd)
        } else {
          mm <- gam(lss ~ s(aa, k = 6), data = dd, method = "REML")
        }

        ss_pred <- predict(mm, newdata = data.frame(aa = 0.44),
                           se.fit = F) |> exp()
        sigma_track$fixef[i, j, k] <- ifelse(ss_pred < 1e-4, 1e-4, ss_pred)
      }
    }
  }

  # steps_extra
  for(j in 1:J2) {
    aa <- accept_track$steps_extra[j, 1:(k-1)]
    lss <- log(sigma_track$steps_extra[j, 1:(k-1)])

    # get previous data
    dd <- data.frame(
      aa = accept_track$steps_extra[j, 1:(k-1)],
      ss = sigma_track$steps_extra[j, 1:(k-1)]
    )

    # remove sigma too close to zero
    dd <- dd[dd$ss >= 1e-4, ]
    dd$lss = log(dd$ss)

    # fit regression
    if(k < 20) {
      mm <- lm(lss ~ aa + I(aa ^ 2), data = dd)
    } else {
      mm <- gam(lss ~ s(aa, k = 6), data = dd, method = "REML")
    }

    ss_pred <- predict(mm, newdata = data.frame(aa = 0.44),
                       se.fit = F) |> exp()
    sigma_track$steps_extra[j, k] <- ifelse(ss_pred < 1e-4, 1e-4, ss_pred)
  }

  # run MCMC
  sss <- list(fixef = sigma_track$fixef[, , k],
              steps_extra = sigma_track$steps_extra[, k])
  run_k <- mcmc(ns, sd_jump = sss, samples = run0_thin)

  # compute and store acceptance
  a_k <- acceptance(run_k)
  accept_track$fixef[, , k] <- a_k$fixef
  accept_track$steps_extra[, k] <- a_k$steps_extra
}

# Fit reverse model (acceptance ~ sigma) and choose that value.
sd_jump_tune <- list(
  fixef = sigma_track$fixef[, , K],
  steps_extra = sigma_track$steps_extra[, K]
)

# fixef
for(j in 1:3) {   # par type
  for(i in 1:4) { # spread parameter

    # i = 4; j = 1 # TEST

    if(!(j == 3 & i < 4)) {
      # get data
      dd <- data.frame(
        aa = accept_track$fixef[i, j, 1:K],
        ss = sigma_track$fixef[i, j, 1:K]
      )

      # fit gam
      mm <- gam(aa ~ s(ss, k = 6), family = betar(),
                data = dd, method = "REML")

      fn <- function(ss) {
        apred <- predict(mm, newdata = data.frame(ss = ss), type = "response")
        return((apred - 0.44) ^ 2)
      }

      opt <- optim(sigma_track$fixef[i, j, K], fn, method = "Brent",
                   lower = min(dd$ss), upper = max(dd$ss))

      sd_jump_tune$fixef[i, j] <- opt$par

      # Visualize
      tit <- paste(rownames(sd_jump1$fixef)[i],
                   colnames(sd_jump1$fixef)[j], sep = "; ")
      ppp <- data.frame(ss = seq(min(dd$ss), max(dd$ss), length.out = 100))
      ppp$y <- predict(mm, ppp, type = "response")
      plot(aa ~ ss, data = dd, main = tit, xlab = "Sigma", ylab = "Acceptance")
      lines(y ~ ss, data = ppp)
      abline(v = opt$par, col = 2, lty = 2)
      abline(h = 0.44, col = 4, lty = 2)
    }
  }
}

# steps_extra
for(j in 1:J2) {
  dd <- data.frame(
    aa = accept_track$steps_extra[j, 1:K],
    ss = sigma_track$steps_extra[j, 1:K]
  )

  # fit gam
  mm <- gam(aa ~ s(ss, k = 6), family = betar(),
            data = dd, method = "REML")

  fn <- function(ss) {
    apred <- predict(mm, newdata = data.frame(ss = ss), type = "response")
    return((apred - 0.44) ^ 2)
  }

  opt <- optim(sigma_track$steps_extra[j, K], fn, method = "Brent",
               lower = min(dd$ss), upper = max(dd$ss))

  sd_jump_tune$steps_extra[j] <- opt$par
}

# compare posterior sd with tune sd
plot(sd_jump_tune$steps_extra ~ sd_post$steps_extra)
mm <- lm(sd_jump_tune$steps_extra ~ sd_post$steps_extra - 1) # b = 2.25
abline(c(0, coef(mm))) # coef(mm) = 2.3

fftune <- as.vector(sd_jump_tune$fixef)
ffpost <- as.vector(sd_post$fixef)
plot(fftune ~ ffpost)
mm <- lm(fftune ~ ffpost - 1)
abline(c(0, coef(mm))) # coef(mm) = 1.2

# saveRDS(sd_jump_tune, file.path("files", "hierarchical_model", "sd_jump_tune"))
sd_jump_tune <- readRDS(file.path("files", "hierarchical_model", "sd_jump_tune"))



# MCMC timing -------------------------------------------------------------

library(microbenchmark)

m1 <- microbenchmark(
  mcmc(nsim = 100, sd_jump = sd_jump_tune, samples = run0_thin),
  times = 30
)
