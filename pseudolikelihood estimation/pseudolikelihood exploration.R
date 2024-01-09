# This code inherits from _compact_support.R, but here I clean it and reduce the
# code to test only the exploration algorithm.

# Tasks:
# * Test how a global search approach deals with a bimodal surface.

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)

library(Rcpp)
library(terra)
library(abind)         # manage arrays

library(randtoolbox)   # sobol sequences
library(GauPro)        # fit gaussian process
library(mgcv)          # fit spline to choose non-bounded particles
library(MGMM)          # fit MVN distribution
library(Matrix)        # nearPD, searches the closest PD matrix.
library(stutils)       # tryNULL, for search of nearPD
library(logitnorm)

library(sn)
library(trialr) # rlkjcorr
library(truncnorm)

library(foreach)       # parallelization
library(doMC)          # doRNG had problems with cpp functions

library(FireSpread)    # spread and similarity functions

library(microbenchmark)

source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for rast_from_mat and a few constants

# source("estimation_functions.R") # prior_dist and other stuff

# Multicore settings -----------------------------------------------------

n_cores <- 15
registerDoMC(n_cores)

# Constants --------------------------------------------------------------

data_dir <- file.path("data", "focal fires data", "landscapes_ig-known")
filenames <- list.files(data_dir)

# dir to save output
target_dir <- file.path("files", "pseudolikelihood_estimation")

# constants for fire spread simulation
upper_limit <- 1
n_coef <- 6
par_names <- c("intercept", "vfi", "tfi", "slope", "wind", "steps")

# number of particles to choose by wave
n_pw <- 1000

# number of fires to simulate by particle
n_sim <- 20

# Functions ---------------------------------------------------------------

normalize <- function(x) x / sum(x)

# prior distribution to simulate parameters or to compute them from a sobol
# sequence (type = "quantile", which computes the icdf.)
prior_dist <- function(mu_int = 0, sd_int = 20,
                       r_slope = 0.04,
                       r_fi = 0.15,
                       r_wind = 0.3,
                       steps_lower = 10,
                       steps_upper = 1000,
                       type = "rng", # or "quantile" or "density" or "prob"
                       n = 1,
                       p = NULL,   # probabilities to compute quantiles
                       x = NULL) { # samples on which density will be evaluated

  if(type == "rng") {
    b <- cbind(
      "intercept" = rnorm(n, mu_int, sd_int),
      "vfi" = rexp(n, r_fi),
      "tfi" = rexp(n, r_fi),
      "slope" = rexp(n, r_slope),
      "wind" = rexp(n, r_wind),
      "steps" = round(runif(n, steps_lower, steps_upper))
    )

    if(n == 1) {
      parnames <- colnames(b)
      b <- as.numeric(b)
      names(b) <- parnames
    }
    return(b)
  }

  if(type == "quantile") {

    if(is.null(p)) stop("Provide probability to compute quantiles.")
    if(length(p) %% n_coef != 0) stop("p must be a multiple of the number of parameters.")

    if(length(dim(p)) < 2) p <- matrix(p, ncol = n_coef)
    if(length(dim(p)) == 2) {
      if(ncol(p) != n_coef) p <- matrix(as.numeric(p), ncol = n_coef)
    }

    q <- cbind(
      "intercept" = qnorm(p[, 1], mu_int, sd_int),
      "vfi" = qexp(p[, 2], r_fi),
      "tfi" = qexp(p[, 3], r_fi),
      "slope" = qexp(p[, 4], r_slope),
      "wind" = qexp(p[, 5], r_wind),
      "steps" = round(qunif(p[, 6], steps_lower, steps_upper))
    )

    return(q)
  }

  if(type == "density") {

    if(is.null(x)) stop("Provide x to compute density.")
    if(is.null(dim(x))) x <- matrix(as.numeric(x), nrow = 1)

    if(ncol(x) %% n_coef != 0) stop("Columns in x must be a multiple of the number of parameters.")

    dmarginals <- cbind(
      "intercept" = dnorm(x[, 1], mu_int, sd_int, log = TRUE),
      "vfi" = dexp(x[, 2], r_fi, log = TRUE),
      "tfi" = dexp(x[, 3], r_fi, log = TRUE),
      "slope" = dexp(x[, 4], r_slope, log = TRUE),
      "wind" = dexp(x[, 5], r_wind, log = TRUE),
      "steps" = dunif(x[, 6], steps_lower, steps_upper, log = TRUE)
    )

    djoint <- apply(dmarginals, 1, sum) %>% as.numeric() %>% exp()
    return(djoint)
  }

  if(type == "prob") {

    if(is.null(x)) stop("Provide x to compute probability")
    if(is.null(dim(x))) x <- matrix(as.numeric(x), nrow = 1)

    if(ncol(x) %% n_coef != 0) stop("Columns in x must be a multiple of the number of parameters.")

    probs <- cbind(
      "intercept" = pnorm(x[, 1], mu_int, sd_int),
      "vfi" = pexp(x[, 2], r_fi),
      "tfi" = pexp(x[, 3], r_fi),
      "slope" = pexp(x[, 4], r_slope),
      "wind" = pexp(x[, 5], r_wind),
      "steps" = punif(x[, 6], steps_lower, steps_upper)
    )

    return(probs)
  }
}

# function to assess if simulated fires reached the landscape bounds.
# returns the number of cells burned at the edge.
edge_count <- function(burned_layer) {
  r <- nrow(burned_layer); c <- ncol(burned_layer)
  rows_edge <- sum(burned_layer[c(1, r), ])
  cols_edge <- sum(burned_layer[-c(1, r), c(1, c)])
  return(rows_edge + cols_edge)
}

# Simulate fires and compare then with the observed one using the overlap_spatial
# function. The landscape argument includes all data to simulate fire, and also
# burned and burned_ids layers.
# Returns a matrix with the mean and variance across replicates of:
#   overlap_spatial,
#   size_diff: size differences, as simulated - observed,
#   edge: number of pixels burned at the edge of the landscape.
# The last two are used to rule out bounded fires.
similarity_simulate_particle <- function(particle, n_sim = 20,
                                         fire_data = NULL) {

  ## testo
  # particle <- particles_sim(N = 1)
  ## end testo
  metrics <- matrix(NA, n_sim, 4)
  colnames(metrics) <- c("overlap", "size_diff", "edge", "steps_used")

  # simulate and compute metrics
  for(i in 1:n_sim) { # simulated fires by particle
    fire_sim <- simulate_fire_compare(
      landscape = fire_data$landscape,
      burnable = fire_data$burnable,
      ignition_cells = fire_data$ig_rowcol,
      coef = particle[1:(n_coef-1)],
      upper_limit = upper_limit,
      steps = particle[n_coef]
    )

    metrics[i, "overlap"] <- overlap_spatial(
      fire_sim, fire_data[c("burned_layer", "burned_ids")]
    )

    metrics[i, "size_diff"] <- ncol(fire_sim$burned_ids) -
      ncol(fire_data$burned_ids)

    metrics[i, "edge"] <- edge_count(fire_sim$burned_layer)

    metrics[i, "steps_used"] <- fire_sim$steps_used
  }

  mo <- mean(metrics[, "overlap"])
  # extract values
  ll_summ <- c(
    overlap = mo,
    overlap_log = log(mo),
    overlap_logit = qlogis(mo),
    overlap_var = var(metrics[, "overlap"]),       # variance in original scale
    ll = qlogis(mo),
    var = var(metrics[, "overlap"]), # same as overlap_var, but kept for compatibility
    size_diff = mean(metrics[, "size_diff"]),
    edge =  mean(metrics[, "edge"]),
    edge_bin_prop = mean(metrics[, "edge"] > 0),
                    # proportion of fires that reached the edge.
    steps_used = mean(metrics[, "steps_used"])
  )

  return(ll_summ)
}

# Emulate the loglik over a list of particles in parallel. Returns the metrics
# in rows.
# It requires the matrix of particles to simulate (coefficient values).
similarity_simulate_parallel <- function(particles_mat = NULL,
                                         fire_data = NULL) {

  # turn particle matrix into list for parallel evaluation
  particles_list <- lapply(1:nrow(particles_mat),
                           function(x) particles_mat[x, ])

  # simulate in parallel
  result <- foreach(pp = particles_list) %dopar% {
    similarity_simulate_particle(pp, fire_data = fire_data)
  }

  # inherit previous dimnames
  names_single <- names(result[[1]])

  # rbind list result
  res <- do.call("rbind", result) %>% as.data.frame()

  names(res) <- names_single

  # add useful columns
  res$wave <- NA
  res$par_values <- particles_mat
  colnames(res$par_values) <- par_names

  return(res)
}

# the same function, but computing only the overlap. x is the parameter vector.
similarity_simulate_simpler <- function(x, n_sim = 20,
                                        fire_data = NULL) {

  ## testo
  # particle <- particles_sim(N = 1)
  ## end testo
  metrics <- numeric(n_sim)

  # simulate and compute metrics
  for(i in 1:n_sim) { # simulated fires by particle
    fire_sim <- simulate_fire_compare(
      landscape = fire_data$landscape,
      burnable = fire_data$burnable,
      ignition_cells = fire_data$ig_rowcol,
      coef = particle[1:(n_coef-1)],
      upper_limit = upper_limit,
      steps = particle[n_coef]
    )

    metrics[i] <- overlap_spatial(fire_sim,
                                  fire_data[c("burned_layer", "burned_ids")])
  }

  return(mean(metrics))
}


# get_bounds: computes the metrics for the largest and smallest fires possible
get_bounds <- function(fire_data, n_sim = 20) {

  coef_burn_all <- c(1e6, rep(0, n_coef - 1))
  coef_burn_none <- c(-1e6, rep(0, n_coef - 1))

  small_fire <- similarity_simulate_particle(coef_burn_none, n_sim = n_sim,
                                             fire_data = fire_data)

  large_fire <- similarity_simulate_particle(coef_burn_all, n_sim = n_sim,
                                             fire_data = fire_data)

  sim_bounds <- rbind(small_fire, large_fire)
  rownames(sim_bounds) <- c("smallest", "largest")

  return(sim_bounds)
}

# function to make new data varying only one predictor.
# the mle, if provided, must be named.
make_newdata <- function(varying = "intercept", data, mle = NULL) {

  values_list_mean <- lapply(par_names, function(v) {
    if(v != varying) {
      res <- mean(data$par_values[, v])
    } else {
      res <- seq(min(data$par_values[, v]),
                 max(data$par_values[, v]),
                 length.out = 150)
    }
    return(res)
  })

  values_list_mle <- lapply(par_names, function(v) {
    if(v != varying) {
      res <- mle[v]
    } else {
      res <- seq(min(data$par_values[, v]),
                 max(data$par_values[, v]),
                 length.out = 150)
    }
    return(res)
  })

  names(values_list_mean) <- names(values_list_mle) <- par_names

  g_mle <- expand.grid(values_list_mle)
  g_mean <- expand.grid(values_list_mean)
  g_mle$pred_type <- "mle"
  g_mean$pred_type <- "mean"

  new_data <- rbind(g_mle, g_mean)

  # add columns indicating which is the varying predictor and which are its
  # values, useful to plot later
  new_data$varying_var <- varying
  new_data$varying_val <- new_data[, varying]

  return(new_data)
}

# partial predictor function
partial_predictions <- function(varying = "intercept", data, loglik_model, mle) {

  ## TEST
  # varying = "all"
  # data = data_pred
  # loglik_model = loglik_model; mle = mle
  ###

  if(varying != "all") {
    new_data <- make_newdata(varying = varying, data = data, mle = mle)

    if(any(class(loglik_model) == "GauPro")) {
      pred <- loglik_model$pred(new_data[, par_names], se.fit = TRUE)

      new_data$mle <- pred$mean
      new_data$upper <- pred$mean + qnorm(0.975) * pred$se
      new_data$lower <- pred$mean - qnorm(0.975) * pred$se
    }

    if(any(class(loglik_model) == "gam")) {
      pred <- predict(loglik_model, newdata = new_data[, par_names],
                      se.fit = TRUE)

      new_data$mle <- pred$fit
      new_data$upper <- pred$fit + qnorm(0.975) * pred$se.fit
      new_data$lower <- pred$fit - qnorm(0.975) * pred$se.fit
    }


  } else {
    new_data <- do.call("rbind", lapply(par_names, function(x) {
      make_newdata(varying = x, data = data, mle = mle)
    }))

    if(any(class(loglik_model) == "GauPro")) {
      pred <- loglik_model$pred(new_data[, par_names], se.fit = TRUE)

      new_data$mle <- pred$mean
      new_data$upper <- pred$mean + qnorm(0.975) * pred$se
      new_data$lower <- pred$mean - qnorm(0.975) * pred$se
    }

    if(any(class(loglik_model) == "gam")) {
      pred <- predict(loglik_model, newdata = new_data[, par_names],
                      se.fit = TRUE)

      new_data$mle <- pred$fit
      new_data$upper <- pred$fit + qnorm(0.975) * pred$se.fit
      new_data$lower <- pred$fit - qnorm(0.975) * pred$se.fit
    }
  }

  return(new_data)
}

# loglik_plot: function to explore the advance of the gp.
loglik_plot <- function(fitting_wave,
                        varying = "all", # parameter to vary in the 1d plot
                        color_point = "in_bounds",
                        response = "overlap",
                        latest_wave = FALSE) {

  ### test
  # fitting_wave = w1
  # varying = "all"
  # color_point = "in_bounds"
  # response = "overlap"
  ###

  loglik_model <- fitting_wave$gps[[length(fitting_wave$gps)]]
  data <- fitting_wave$particles_data
  data$in_bounds <- factor(as.character(data$in_bounds), levels = c("0", "1"))

  # get current wave and add to data
  current <- which(data$wave == max(data$wave))
  old <- which(data$wave < max(data$wave))

  data$wave_plot <- "previous"
  data$wave_plot[current] <- "current"

  if(latest_wave) {
    data_pred <- data[current, ]
  } else {
    data_pred <- data
  }

  # get mle for partial predictions
  mle <- fitting_wave$loglik_optim$par

  # compute partial predictions
  preds <- partial_predictions(varying = varying,
                               data = data_pred,
                               loglik_model = loglik_model, mle = mle)

  if(varying != "all") {
    # bring out of the matrix the par_values for plotting
    data_plot <- do.call(data.frame, data)
    cols_change <- grep("par_values", names(data_plot))
    names(data_plot)[cols_change] <- par_names

    # plot
    p <-
      ggplot() +

      # data
      geom_point(data = data_plot,
                 mapping = aes_string(x = varying, y = response,
                                      color = color_point),
                 size = 2, alpha = 0.5) +
      scale_color_viridis(end = 0.7, discrete = TRUE) +

      # loglik_model predictions
      ggnewscale::new_scale_color() +
      scale_color_viridis(end = 0.7, discrete = TRUE, option = "A") +
      ggnewscale::new_scale_fill() +
      scale_fill_viridis(end = 0.7, discrete = TRUE, option = "A") +

      geom_ribbon(data = preds, mapping = aes_string(
        x = varying, y = "mle", ymin = "lower", ymax = "upper",
        fill = "pred_type"
      ), alpha = 0.3, color = NA) +

      geom_line(data = preds, mapping = aes_string(x = varying, y = "mle",
                                                   color = "pred_type")) +

      theme(panel.grid.minor = element_blank(),
            legend.position = "right") +#c(0.85, 0.15)) +
      ylab(response) +
      xlab("parameter value")

    print(p)
    return(list(plot = p, data_points = data_plot, data_preds = preds))

  } else {

    # when all parameters are varied, data must be replicated to be plotted
    # against every parameter.

    data_expand <- do.call("rbind", lapply(par_names, function(n) {
      data$varying_val <- data$par_values[, n]
      data$varying_var <- n
      return(data)
    }))

    data_expand$varying_var <- factor(data_expand$varying_var,
                                      levels = par_names)

    preds$varying_var <- factor(preds$varying_var,
                                levels = par_names)

    p <-
      ggplot() +

      # data
      geom_point(data = data_expand,
                 mapping = aes_string(x = "varying_val", y = response,
                                      color = color_point),
                 size = 2, alpha = 0.5) +
      scale_color_viridis(end = 0.7, discrete = TRUE) +

      # loglik_model predictions
      ggnewscale::new_scale_color() +
      scale_color_viridis(end = 0.7, discrete = TRUE, option = "A") +
      ggnewscale::new_scale_fill() +
      scale_fill_viridis(end = 0.7, discrete = TRUE, option = "A") +

      geom_ribbon(data = preds, mapping = aes(
        x = varying_val, y = mle, ymin = lower, ymax = upper, fill = pred_type
      ), alpha = 0.3, color = NA) +

      geom_line(data = preds, mapping = aes(x = varying_val, y = mle,
                                            color = pred_type)) +

      facet_wrap(vars(varying_var), scales = "free_x") +
      theme(panel.grid.minor = element_blank(),
            legend.position = "right") +#c(0.85, 0.15)) +
      ylab(response) +
      xlab("parameter value")

    print(p)
    return(list(plot = p, data_points = data_expand, data_preds = preds))
  }
}

wave_plot <- function(data, response = "overlap", alpha = 0.3, best = NULL) {

  if(!is.null(best)) {
    data <- data[order(data[, response], decreasing = TRUE), ]
  }

  yy <- range(data[, response])
  yy[2] <- yy[2] * 1.05

  # title for first plot
  tit <- deparse(substitute(data))

  par(mfrow = c(3, 2))
  for(i in 1:n_coef) {
    mm <- ifelse(i == 1, tit, NA)
    plot(data[, response] ~ data$par_values[, i], ylab = response,
         xlab = par_names[i], ylim = yy, main = mm,
         pch = 19, col = rgb(0, 0, 0, alpha))

    if(!is.null(best)) {
      points(data[1:best, response] ~ data$par_values[1:best, i],
           pch = 19, col = rgb(0, 0, 1, alpha))
    }
  }
  par(mfrow = c(1, 1))
}


# https://rpubs.com/binhho660/922614
#' rmvn by hand, Covariance parameterization, or maybe from a sobol sequence
#' @param n number of samples to draw
#' @param mu p-vector mean or p-columns matrix of means.
#' @param Sigma covariance matrix (p x p)
#' @param sobol logical, indicate whether to use or not a sobol sequence.
#' @return matrix of dimension p x n of samples
rmvn_sobol <- function(n, mu, Sigma, sobol = F){
  if(is.null(dim(mu))) {
    mu <- matrix(mu, nrow = 1)
  }

  p <- ncol(mu)

  if(sobol) {
    if(nrow(mu) > 1) stop("Sobol sequence not recommended with more than one center.")
    P <- sobol(n, p, init = F)
    Z <- apply(P, 2, qnorm) %>% t()
  } else {
    Z <- matrix(rnorm(p*n), p, n)
  }

  L <- t(chol(Sigma)) # By default R's chol function returns upper cholesky factor
  X <- L %*% Z

  if(nrow(mu) == 1) { # global centre
    mu <- as.numeric(mu)
    X <- sweep(X, 1, mu, FUN = `+`)
  } else {            # focal centre
    X <- X + t(mu)
  }
  return(X)
}


# positive-definiteness check made by makeSECdistr {sn}. It differs from
# is.positive.definite(), from {matrixcalc}
is_positive_definite <- function(m) {
  min(eigen(m, symmetric = TRUE, only.values = TRUE)$values) > 0
}

# make_positive_definite makes a vcov positive definite. It's supposed that if
# the correlation matrix is pd, then the vcov contructed with any set of
# scales would be. But it isn't so in practice.
# https://stats.stackexchange.com/questions/23317/create-positive-definite-3x3-covariance-matrix-given-specified-correlation-value

# If the positive-definiteness fails after max_iter trials, the result is NA.
# This is useful for optim(), because otherwise the iteration number for finding
# a positive-definite matrix took too long. It's safer to just ignore problematic
# parameter vector
make_positive_definite <- function(m, corr = FALSE, max_iter = 1e4) {
  if(is_positive_definite(m)) {
    return(m)
  } else {
    temp <- tryNULL(nearPD(m, corr = corr, keepDiag = TRUE, maxit = max_iter))
    if(is.null(temp)) return(NA) else {
      mat <- as.matrix(temp$mat)
      if(!is_positive_definite(mat)) {
        return(NA)
      } else return(mat)
    }
  }
}

    # mvn_fitted <- gam(
    #   list(intercept ~ 1, vfi ~ 1, tfi ~ 1, slope ~ 1, wind ~ 1, steps ~ 1),
    #   family = mvn(n_coef), data = as.data.frame(xun_seeds)
    # )
    # mus <- coef(mvn_fitted)[1:n_coef]
    # V <- solve(crossprod(mvn_fitted$family$data$R))

    # # Fit MVN by hand
    # V <- cov(xun_seeds)
    # mus <- apply(xun_seeds, 2, mean)
    #
    # # Ensure positive definiteness of enlarged fitted covariance matrix
    # Vlarge <- V * var_factor
    # Vlarge <- make_positive_definite(as.matrix(forceSymmetric(Vlarge)))
    #
    # # if we did not find a positive definite matrix, use a diagonal vcov.
    # if(is.null(dim(Vlarge))) { # returns "NA" when could not find it
    #   Vlarge <- diag(apply(xun_seeds, 2, var) * var_factor)
    #   # again force numerical PDness
    #   diag(Vlarge) <- diag(Vlarge) + 1e-12
    # }
    #
    # # sample new particles using enlarged covariance
    # xun_new <- apply(xun_seeds, 1, FUN = function(x) {
    #   mgcv::rmvn(n = 1, mu = mus, V = Vlarge)
    # }) %>% t




# invlogit_scaled: translates unconstrained variables (logit) to a flat compatc
#   support, by applying invlogit and then scaling.
# the support matrix has dimensions [c(lower, upper), var_names]
invlogit_scaled <- function(xun, support) {
  x <- xun
  for(j in 1:ncol(xun)) {
    x[, j] <- plogis(xun[, j]) * (support[2, j] - support[1, j]) + support[1, j]
  }
  return(x)
}

# logit_scaled: translates variables with compact support to the logit space,
#   by previously scaling them between 0 and 1.
# the support matrix has dimensions [c(lower, upper), var_names]
logit_scaled <- function(x, support) {
  xun <- x
  for(j in 1:n_coef) {
    xun[, j] <- qlogis((x[, j] - support[1, j]) / (support[2, j] - support[1, j]))
  }
  return(xun)
}

# function to scale [0, 1] particles to a given support
scale_params <- function(x, support) {
  xs <- x
  for(i in 1:n_coef) {
    xs[, i] <- x[, i] * (support[2, i] - support[1, i]) + support[1, i]
  }
  return(xs)
}


# explore_likelihood searches particles in a space with compact support.
# The algorithm starts with a dataset of simulated likelihood over a sobol
# sequence of particles, and then it resamples the best ones according to the
# likelihood. Good particles are reproduced with a MVN kernel.
# data: dataset with previous simulated particles.
# n: number of new particles to simulate.
# var_factor: multiplier of the covariance matrix computed from the best
#   particle. Ideally, between 1 and 2.
# n_best: number of best particles to reproduce. if "all", all data available
#   are resampled.
# p_best: the proportion of best particles to resample. Overrides n_best.
# weight_power: use likelihood ^ weight_power to resample good particles.
# support: matrix with lower limits (first row) and upper limits (second row)
#   of all variables.
# spread_data: data to simulate fires.
# centre: "local" or "global". New particles are simulated from a MVN centred
#   at each resampled particle (focal), or at the global mean (global). If global,
#   sobol may be set to TRUE, if focal, sobol is set to FALSE.
explore_likelihood <- function(data, n = 600, var_factor = 2, n_best = "all",
                               p_best = 0.1,
                               weight_power = 1, support, spread_data,
                               centre = "local", sobol = TRUE) {

  # which particles will be resampled?
  if(!is.null(p_best)) {
    n_best <- round(nrow(data) * p_best)
    if(n_best < 500) {
      n_best <- 500
    }
  }
  if(n_best == "all") {
    n_best <- nrow(data)
  } else {
    data <- data[order(data$overlap, decreasing = TRUE), ]
  }

  # resample the best particles
  ids_rep <- sample(1:n_best, size = n, replace = T,
                    prob = data$overlap[1:n_best] ^ weight_power)
  dexp <- data[ids_rep, ]

  # transform to unconstrained scale
  dexp$par_values_raw <- logit_scaled(dexp$par_values, support)

  # MVN kernel
  # first, remove infinite values and redefine n
  finite <- apply(dexp$par_values_raw, 1, function(x) all(is.finite(x))) %>% t
  dexp <- dexp[finite, ]
  n <- nrow(dexp)

  V <- make_positive_definite(cov(dexp$par_values_raw))
  if(anyNA(V)) V <- diag(apply(dexp$par_values_raw, 2, var))
  Vlarge <- V * var_factor

  if(centre == "local") {
    mus <- dexp$par_values_raw # centred at good particles
    sobol <- FALSE
  }
  if(centre == "global") {
    mus <- apply(dexp$par_values_raw, 2, mean) # centred at global mean
  }

  candidates_raw <- rmvn_sobol(n, mus, Vlarge, sobol) %>% t
  # transform to original scale
  colnames(candidates_raw) <- par_names
  candidates <- invlogit_scaled(candidates_raw, support)

  # message("Simulating fires")
  sim_result <- similarity_simulate_parallel(particles_mat = candidates,
                                             fire_data = spread_data)

  # define wave
  sim_result$wave <- max(data$wave) + 1
  colnames(sim_result$par_values) <- par_names

  # merge old and new datasets
  res <- rbind(
    data[, colnames(sim_result)],
    sim_result
  )

  return(res)
}

# iterate waves of explore_likelihood()
explore_likelihood_iterate <- function(
    n_waves = 10,
    data, n = 500, var_factor = 1.5, n_best = "all", p_best = NULL,
    weight_power = 1, support, spread_data,
    centre = "local", sobol = TRUE
) {

  last_wave <- max(data$wave)
  m <- paste("Search starts at max overlap = ",
             round(max(data$overlap), 4), sep = "")
  message(m)

  # initialize res as previous data
  res <- data

  # create var_factor sequence for the variable case
  var_factors <- rep(var_factor, ceiling(n_waves / length(var_factor)))

  # # create centre sequence for the variable case
  # each <- ceiling(length(var_factor) / 2) # var_factors are for each centre
  # centres_dup <- rep(centre, each = each)
  # centres <- rep(centres_dup, ceiling(n_waves / length(var_factor)))

  for(i in 1:n_waves) {
    res <- explore_likelihood(data = res, n = n, var_factor = var_factors[i],
                              n_best = n_best, p_best = p_best,
                              weight_power = weight_power,
                              support = support, spread_data = spread_data,
                              centre = centre, sobol = sobol)

    m <- paste("Search wave ", last_wave + i, "; max overlap = ",
               round(max(res$overlap), 4), sep = "")
    message(m)
  }

  return(res)
}


# the same functions but made for testing a known function

explore_test <- function(data, n = 600, var_factor = 2, n_best = "all",
                         p_best = 0.1, fn = NULL,
                         support, centre = "local",
                         sobol = TRUE,
                         ...) {

  # data = sim1
  # n = 600
  # var_factor = 1
  # n_best = "all"
  # p_best = NULL
  # fn = dsn_uni
  # support = support
  # centre = "local"
  # sobol = TRUE
  # dp = dp

  # which particles will be resampled?
  if(!is.null(p_best)) {
    n_best <- round(nrow(data) * p_best)
    if(n_best < 500) {
      n_best <- 500
    }
  }
  if(n_best == "all") {
    n_best <- nrow(data)
  } else {
    data <- data[order(data$y, decreasing = TRUE), ]
  }

  # resample the best particles
  ids_rep <- sample(1:n_best, size = n, replace = T,
                    prob = data$y[1:n_best])
  dexp <- data[ids_rep, ]

  # transform to unconstrained scale
  dexp$par_values_raw <- logit_scaled(dexp$par_values, support)

  # MVN kernel
  # first, remove infinite values and redefine n
  finite <- apply(dexp$par_values_raw, 1, function(x) all(is.finite(x))) %>% t
  dexp <- dexp[finite, ]
  n <- nrow(dexp)

  V <- make_positive_definite(cov(dexp$par_values_raw))
  if(anyNA(V)) V <- diag(apply(dexp$par_values_raw, 2, var))
  Vlarge <- V * var_factor

  # check we can find the cholesky factor
  tt <- tryNULL(chol.default(Vlarge))
  if(is.null(tt)) {
    sss <- apply(dexp$par_values_raw, 2, var)
    if(any(sss < 0.1)) {sss <- sss + 0.5}
    Vlarge <- diag(sss * var_factor)
  }


  if(centre == "local") {
    mus <- dexp$par_values_raw # centred at good particles
    sobol <- FALSE
  }
  if(centre == "global") {
    mus <- apply(dexp$par_values_raw, 2, mean) # centred at global mean
  }

  candidates_raw <- rmvn_sobol(n, mus, Vlarge, sobol) %>% t
  # transform to original scale
  colnames(candidates_raw) <- par_names
  candidates <- invlogit_scaled(candidates_raw, support)

  # evaluate target function
  sim_result <- fn(dexp$par_values, ...)

  # define wave
  sim_result$wave <- max(data$wave) + 1

  # merge old and new datasets
  res <- rbind(
    data[, colnames(sim_result)],
    sim_result
  )

  return(res)
}

# iterate waves of explore_likelihood()
explore_test_iterate <- function(
    n_waves = 10,
    data, n = 500, var_factor = 1.5, n_best = "all", p_best = NULL,
    weight_power = 1, support,
    centre = "global", sobol = TRUE,
    fn = NULL, ...
) {

  last_wave <- max(data$wave)
  m <- paste("Search starts at max y = ",
             round(max(data$y), 4), sep = "")
  message(m)

  # initialize res as previous data
  res <- data

  # create var_factor sequence for the variable case
  var_factors <- rep(var_factor, ceiling(n_waves / length(var_factor)))

  # # create centre sequence for the variable case
  # each <- ceiling(length(var_factor) / 2) # var_factors are for each centre
  # centres_dup <- rep(centre, each = each)
  # centres <- rep(centres_dup, ceiling(n_waves / length(var_factor)))

  for(i in 1:n_waves) {
    res <- explore_test(data = res, n = n, var_factor = var_factors[i],
                        n_best = n_best, p_best = p_best,
                        support = support,
                        centre = centre, sobol = sobol,
                        fn = fn, ...)

    m <- paste("Search wave ", last_wave + i, "; max y = ",
               max(res$y), sep = "")
    message(m)
  }

  return(res)
}





# Data for spread ---------------------------------------------------------

fire_file <- "2008_3.rds" #"2012_53.rds"#2008_3.rds" # "1999_25j.rds"
full_data <- readRDS(file.path(data_dir, fire_file))
# subset data needed for spread (to be cloned across workers)
spread_data <- full_data[c("landscape", "burnable", "ig_rowcol",
                           "burned_layer", "burned_ids")]

bb <- get_bounds(fire_data = spread_data)
ext_alpha <- 50
ext_beta <- 50
params_lower <- c("intercept" = -ext_alpha,
                  "vfi" = 0,
                  "tfi" = 0,
                  "slope" = 0,
                  "wind" = 0,
                  "steps" = 10)
params_upper <- c("intercept" = ext_alpha,
                  "vfi" = ext_beta,
                  "tfi" = ext_beta,
                  "slope" = ext_beta * 2,
                  "wind" = ext_beta,
                  "steps" = bb["largest", "steps_used"])
sup <- rbind(params_lower, params_upper)


# Exploring likelihood for one fire, varying method ---------------------

# first wave
ss <- sobol(n = 600, dim = n_coef)
particles <- scale_params(ss, sup)
sim1 <- similarity_simulate_parallel(particles, spread_data)
sim1$wave <- 0
colnames(sim1$par_values) <- par_names
# wave_plot(sim1)
# range(sim1$overlap)
# pairs(sim1$par_values[sim1$overlap > 0.3, ])


# Cosas que hay que probar:
# * definir n_best. Resample all data points?
# * var_factor?
# * weight_power?
# * centre focal or centre global?
# * use sobol for global centre?


# con centre global #

# n_best es nrow(sim1), focal global, using sobol, var_factor = 1.5, weight_po = 1
w1 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 1.5, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w1, alpha = 0.01)
wave_plot(w1[w1$wave <= 10, ], alpha = 0.08)

# n_best es nrow(sim1), focal global, using sobol, var_factor = 1.5, weight_po = 2
w2 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 1.5, weight_power = 2, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w2, alpha = 0.01)
wave_plot(w2[w2$wave <= 10, ], alpha = 0.08)
# muy similar con pow = 1 o 2.


# variamos var_factor.

# n_best es nrow(sim1), focal global, using sobol, var_factor = 1, weight_po = 1
w3 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 1, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w3, alpha = 0.01)
wave_plot(w3[w3$wave <= 10, ], alpha = 0.08)
# con var_factor = 1 parece acotar un poco el espacio, le come la parte de
# intercepts bajos

# n_best es nrow(sim1), focal global, using sobol, var_factor = 2, weight_po = 1
w4 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 2, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w4, alpha = 0.01)
wave_plot(w4[w4$wave <= 10, ], alpha = 0.08)
# con var_factor = 2, las propuestas se van al borde del soporte y se acumulan
# ahí. En este caso dificulta caracterizar el máximo.
# Por otro lado, acá veo que hay valores altos para intercepts mínimos, y antes
# no estaba llegando a esa parte. Entonces quizás sí convenga en algún momento
# aplicar var_factor = 2. O una mezcla?
# Y el var_factor = 1 se concentra mucho en el máximo.

# conclu con centro global:
# weigth_pow = 1 va bien, y quizás se podría combinar una búsqueda variando el
# var_factor en {1, 1.5, 2}.



# con centre focal #

# variar varianza

# n_best es nrow(sim1), using sobol, var_factor = 1.5, weight_po = 1
w5 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 1.5, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "local", sobol = TRUE
)
wave_plot(w5, alpha = 0.01)
# se queda en los extremos mal, habría que bajarle la varianza.
# claro, focal se optimiza con varianzas menores que el global.

# n_best es nrow(sim1), using sobol, var_factor = 1, weight_po = 1
w6 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 1, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "local", sobol = TRUE
)
wave_plot(w6, alpha = 0.01)
# muy mal, achicar más.

# n_best es nrow(sim1), var_factor = 0.7, weight_po = 1
w7 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 0.25, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "local", sobol = TRUE
)
wave_plot(w7, alpha = 0.05)
# con var_factor = 0.25 se queda en los máximos. Probemos 0.5

# n_best es nrow(sim1), var_factor = 0.7, weight_po = 1
w77 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 0.35, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "local", sobol = TRUE
)
wave_plot(w77, alpha = 0.05)
# con var_factor = 0.35 se queda en los máximos. Probemos 0.5

# n_best es nrow(sim1), using sobol, var_factor = 0.5, weight_po = 1
w8 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 0.45, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "local", sobol = TRUE
)
wave_plot(w8, alpha = 0.05)

w88 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 0.5, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "local", sobol = TRUE
)
wave_plot(w88, alpha = 0.05)

w888 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 0.75, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "local", sobol = TRUE
)
wave_plot(w888, alpha = 0.05)

colnames(w8$par_values) <- par_names
wave_plot(w8[w8$par_values[, "wind"]  < 30, ], alpha = 0.5)
# se acumula bastante en los extremos, pero anda muy bien.
# Quizás es combinable con 0.25

# n_best es nrow(sim1), using sobol, var_factor = 0.5, weight_po = 1
w9 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 0.35, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "local", sobol = TRUE
)
wave_plot(w9, alpha = 0.01)
# se compaorta un poco mejor que 0.25, pero aún se acumulan datos en los
# extremos.

# focal y global tienen un desempeño similar. En ambos hay que tunear la
# var_factor para que no se acumule ni en el centro ni en los extremos.

# Quizás ahora haya que probar qué pasa si se remuestrean todos los datos,
# no solo los mejores. Quizás eso genere más datos con overlap intermedio.
# Otra opción es remuestrear siempre el 10 % mejor, o algún otro porcentaje.


## OJO, acá siempre estuve usando como n_best = nrow(sim1),
## en vez de usar nrow(data). Quizás por eso se acumulan en las puntas.

# variamos n_best, usando la w3 como ref, que tiene propuesta global con
# var_factor = 1


# n_best es nrow(data), focal global, using sobol, var_factor = 1, weight_po = 1
w10 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = "all",
  var_factor = 1, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w10, alpha = 0.01)
# no llega al máximo como antes, se queda re abajo. Probar usando percentiles.

# n_best = 0.25 * nrow(data), focal global, using sobol, var_factor = 1, weight_po = 1
w11 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, p_best = 0.25,
  var_factor = 1, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w11, alpha = 0.01)
wave_plot(w3, alpha = 0.01)
# con p_best = 0.25 cumula muchos más datos arriba, quizás de una forma un poco
# más asimétrica que con n_best = 600.
# Si lo que nos importa más es el pico, quizás está bien usar p_best = alto.


# n_best = 0.5 * nrow(data), focal global, using sobol, var_factor = 1, weight_po = 1
w12 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, p_best = 0.5,
  var_factor = 1, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w12, alpha = 0.01)
# no, con p_best = 0.5 es malísimo.

# n_best = 0.1 * nrow(data), focal global, using sobol, var_factor = 1, weight_po = 1
w13 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, p_best = 0.1,
  var_factor = 1, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w13, alpha = 0.01)
# demasiados puntos en el tope


# combining var_factors, using almost w11 as reference (p_best = 0.20)

w14 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, p_best = 0.1, # pucha, quería 0.2
  var_factor = c(1, 1, 2),
  weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w14, alpha = 0.01)
# llega al max re rápido


w15 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, p_best = 0.25,
  var_factor = c(1, 2),
  weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w15, alpha = 0.05)
# llega al max re rápido

# con tres valores de factor
w16 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, p_best = 0.25,
  var_factor = c(1.5, 2, 1),
  weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w16, alpha = 0.01)

# con p 0.15
w17 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, p_best = 0.15,
  var_factor = c(1.5, 2, 1),
  weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w17, alpha = 0.01)


# Bien, esta, la 17, parece bastante buena. Y qué pasa si simulamos 10000 más?
w17b <- explore_likelihood_iterate(
  n_waves = 20, data = w17, n = 500, p_best = 0.15,
  var_factor = c(1.5, 2, 1),
  weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w17b, alpha = 0.05)
# con 10600 anduvo bien. Al agregar esas 10000 no mejoró casi nada.
wave_plot(w17[1:6600, ], alpha = 0.2)
# pero con menos parece no completar todo el espacio




## local, combining var_factors
wl1 <-  explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = c(0.4, 0.3, 0.7), weight_power = 1, support = sup,
  spread_data = spread_data, centre = "local", sobol = F
)
wave_plot(wl1, alpha = 0.05)


wl2 <-  explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = c(0.25, 0.35, 0.7), weight_power = 1, support = sup,
  spread_data = spread_data, centre = "local", sobol = F
)
wave_plot(wl2, alpha = 0.05)

wl3 <-  explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, p_best = 0.15,
  var_factor = c(0.25, 0.35, 0.7), weight_power = 1, support = sup,
  spread_data = spread_data, centre = "local", sobol = F
)
wave_plot(wl3, alpha = 0.05)


# la global combinada (w17) parece ser la mejor para explorar el espacio.
# Queda probar con otras funciones conocidas cómo se desempeñan el global
# y el local.
# quizás si hay una región con 3 modos no es tan buena. O sí?


# Probamos lo mismo con otros incendios ---------------------------------

fire_file <- "2012_58.rds"#"1999_28.rds"#"1999_25j.rds"#"2008_3.rds"
full_data <- readRDS(file.path(data_dir, fire_file))
# subset data needed for spread (to be cloned across workers)
spread_data <- full_data[c("landscape", "burnable", "ig_rowcol",
                           "burned_layer", "burned_ids")]

bb <- get_bounds(fire_data = spread_data)
ext_alpha <- 20
ext_beta <- 15
params_lower <- c("intercept" = -ext_alpha,
                  "vfi" = 0,
                  "tfi" = 0,
                  "slope" = 0,
                  "wind" = 0,
                  "steps" = 10)
params_upper <- c("intercept" = ext_alpha,
                  "vfi" = ext_beta,
                  "tfi" = ext_beta,
                  "slope" = ext_beta * 2,
                  "wind" = ext_beta,
                  "steps" = bb["largest", "steps_used"])

sup <- rbind(params_lower, params_upper)

# first wave
ss <- sobol(n = 600, dim = n_coef)
particles <- scale_params(ss, sup)
sim1 <- similarity_simulate_parallel(particles, spread_data)
sim1$wave <- 0

# Balcon del gutierrez
w17 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, p_best = 0.15,
  var_factor = c(1.5, 2, 1),
  weight_power = 1, support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(w17, alpha = 0.05)

# con balcón del gutierrez (1999_25j) son re planas, y anda bien con muchos steps.
# Debe ser porque la zona quemable que tiene es muy limitada.

# 1999_28 también tiene likes planas. Sólo importa una combn de viento e intercepts.
# Claramente steps se lleva todo.

# Empiezan con overlaps altísimos.
# 2012_58 muestra eff de la elevation.
# Y con steps altos, qué vars parecen importar?
wave_plot(w17[w17$par_values[, 6] > 250, ], alpha = 0.3)
# ahí se ve bien claro el efecto del intercept

saveRDS(w17, "files/pseudolikelihood_estimation/likelihood_points_2015_58.rds")

# Quizás para entender lo de drivers haga falta fijar el steps en el máximo
# posible, así le dejamos que frene por otros motivos.

# Pero digamos que el steps alto hace que tengas que mirar detalles más finos,
# porque lo grueso es determinado por el clima... ahí simplemente usamos un
# kernel de aceptación muy astringente y es como hacerle zoom a la puntita del
# overlap.
# Incluso así vamos a tener menor sensibilidad, pero no tantísima.



# Test exploration with known functions -----------------------------------
# unimodal skew-mvn -------------------------------------------------------

d <- 6
limit <- 10
support <- rbind(rep(-limit, d), rep(limit, d))

# simulate parameters to create a distribution
xis <- rnorm(d, sd = 3)
alphas <- rnorm(d, sd = 10)
sigmas <- rtruncnorm(d, a = 0, mean = 5, sd = 3)
rho <- rlkjcorr(1, d, eta = 0.2)
omega <- rho * (sigmas %*% t(sigmas))
dp <- list(xi = xis, Omega = omega, alpha = alphas)
smvn <- makeSECdistr(dp, family = "SN")
plot(smvn)

# density for unimodal skew-mvn
dsn_uni <- function(x, dp) { # x is a matrix, dp is the dp list
  dd <- data.frame(y = dmsn(x, dp = dp))
  dd$par_values <- x
  return(dd)
}

# sobol seq
ss <- sobol(n = 10000, dim = d)
particles <- scale_params(ss, support)
sim1 <- dsn_uni(particles, dp)
sim1$wave <- 0


t1 <- explore_test_iterate(n_waves = 10, data = sim1, n = 10000, var_factor = 1,
                           p_best = 0.15, fn = dsn_uni,
                           support = support, centre = "global",
                           sobol = TRUE, dp = dp)
nrow(t1)
t1sub <- t1[order(t1$y, decreasing = T)[1:10000], ]
wave_plot(t1sub, response = "y")
summary(t1$y)
# como la función es muy picuda, se concentra en un punto y no muestrea más que
# eso. habría que tunearle el kernel de reproducción, pero da fiaca.
# Vamos a comparar directamente el MCMC usando el surrogate y sin usarlo.

# TAREA -------------------------------------------------------------------

# Testear la exploración con funciones conocidas.
# Alternar búsqueda global con local, por si hay modos pequeños y finitos.
# Tunear la seq de varianzas locales.
# Ver cómo combinar local-global.
# Agrandar límites a [-50, 50].


# Comparing simulation-based y surrogate-based MCMC -----------------------











