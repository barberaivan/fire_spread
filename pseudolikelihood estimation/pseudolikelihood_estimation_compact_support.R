# This code inherits from _smc_02.R, and here I develop the search for a
# likelihood function with compact support.
# There remain unnecessary functions. The code should be clened.

# Ideas to test:

# * use MVN kernel to reproduce particles
# * tune sigma_factor of that kernel
# * try weighting the particles differently, for example:
#   prob = overlap ^ weight_power,
#   with weight_power ~ 2

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

  return(res)
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

  par(mfrow = c(3, 2))
  for(i in 1:n_coef) {
    plot(data[, response] ~ data$par_values[, i], ylab = response,
         xlab = par_names[i], ylim = yy,
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

  # create centre sequence for the variable case
  each <- length(var_factor) / 2 # var_factors are for each centre
  centres_dup <- rep(centre, each = each)
  centres <- rep(centres_dup, ceiling(n_waves / length(var_factor)))

  for(i in 1:n_waves) {
    res <- explore_likelihood(data = res, n = n, var_factor = var_factors[i],
                              n_best = n_best, p_best = p_best,
                              weight_power = weight_power,
                              support = support, spread_data = spread_data,
                              centre = centres[i], sobol = sobol)

    m <- paste("Search wave ", last_wave + i, "; max overlap = ",
               round(max(res$overlap), 4), sep = "")
    message(m)
  }

  return(res)
}


# now consider the kernel density to weight new particles
explore_likelihood_kernel <- function(data, n = 600, var_factor = 2,
                                      n_best = 500, p_best = NULL,
                                      support, spread_data) {

  # ## TEST INPUTS
  #
  # # first wave
  # ss <- sobol(n = 1000, dim = n_coef)
  # particles <- scale_params(ss, sup)
  # sim1 <- similarity_simulate_parallel(particles, spread_data)
  #
  # # prior at logit scale
  # particles_raw <- logit_scaled(particles, sup)
  # sim1$prior <- logistic_prior(particles_raw) # at logit scale!!
  # sim1$kernel <- 1
  # sim1$weight <- sim1$overlap * sim1$prior / sim1$kernel    ## prior * overlap / 1
  #
  # sim1$wave <- 0
  #
  # # fn arguments
  # data = sim1
  # n = 500
  # var_factor = 2
  # n_best = 500
  # p_best = NULL
  # weight_power = 1
  # support = sup
  # spread_data = spread_data
  # ## END TEST INPUTS

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
                    prob = data$weight[1:n_best])
  dexp <- data[ids_rep, ]
  # as seed particles are resampled, their weight is taken as 1 to compute
  # the kernew for their draws.

  # transform to unconstrained scale
  dexp$par_values_raw <- logit_scaled(dexp$par_values, support)

  # create new samples from seeds
  finite <- apply(dexp$par_values_raw, 1, function(x) all(is.finite(x))) %>% t
  dexp <- dexp[finite, ]
  n <- nrow(dexp)

  # Independent Normal Kernels
  vv <- apply(dexp$par_values_raw, 2, var) * var_factor
  candidates_raw <- sapply(1:n_coef, function(j) {
    rnorm(n) * vv[j] + dexp$par_values_raw[, j]
  })

  # transform to original scale
  candidates <- invlogit_scaled(candidates_raw, support)
  colnames(candidates) <- par_names

  # message("Simulating fires")
  sim_result <- similarity_simulate_parallel(particles_mat = candidates,
                                             fire_data = spread_data)

  # get prior, kernel and weight
  sim_result$prior <- logistic_prior(candidates_raw)
  sim_result$kernel <- normal_kernel(candidates_raw, dexp$par_values_raw,
                                     sqrt(vv))
  sim_result$weight <- sim_result$overlap * sim_result$prior / sim_result$kernel

  # define wave
  sim_result$wave <- max(data$wave) + 1

  # merge old and new datasets
  res <- rbind(
    data[, names(sim_result)],
    sim_result
  )

  return(res)
}


# iterate waves of explore_likelihood_kernel()
explore_likelihood_kernel_iterate <- function(
    n_waves = 10, data, n = 500, var_factor = 2, n_best = 500, p_best = NULL,
    support, spread_data
) {

  last_wave <- max(data$wave)
  m <- paste("Search starts at max overlap = ",
             round(max(data$overlap), 4), sep = "")
  message(m)

  # initialize res as previous data
  res <- data

  for(i in 1:n_waves) {
    res <- explore_likelihood_kernel(
      data = res, n = n, var_factor = var_factor, n_best = n_best,
      p_best = p_best, support = support, spread_data = spread_data
    )

    m <- paste("Search wave ", last_wave + i, "; max overlap = ",
               round(max(res$overlap), 4), sep = "")
    message(m)
  }

  return(res)
}

logistic_prior <- function(x) {
  margins <- apply(x, 2, dlogis, log = T)
  prob <- rowSums(margins) %>% exp
  return(prob)
}

normal_kernel <- function(draws, seeds, ss) {

  m <- matrix(NA, nrow(draws), nrow(seeds))
  for(i in 1:nrow(draws)) {
    for(j in 1:nrow(seeds)) {
      m[i, j] <- dnorm(draws[i, ], seeds[j, ], ss, log = TRUE) %>% sum %>% exp
    }
  }

  mm <- rowSums(m)
  # scale, so sum(p) == nrow(draws)
  mmm <- mm * length(mm) / sum(mm)

  return(mmm)
}


gam_formula <- formula(
  y ~

    # marginal effects
    s(intercept, bs = basis, m = m, k = k_side) +
    s(vfi, bs = basis, m = m, k = k_side) +
    s(tfi, bs = basis, m = m, k = k_side) +
    s(slope, bs = basis, m = m, k = k_side) +
    s(wind, bs = basis, m = m, k = k_side) +
    s(steps, bs = basis, m = m, k = k_side) +

    # interactions
    ti(intercept, vfi, k = k_int, bs = basis, m = m) +
    ti(intercept, tfi, k = k_int, bs = basis, m = m) +
    ti(intercept, slope, k = k_int, bs = basis, m = m) +
    ti(intercept, wind, k = k_int, bs = basis, m = m) +
    ti(intercept, steps, k = k_int, bs = basis, m = m) +

    ti(vfi, tfi, k = k_int, bs = basis, m = m) +
    ti(vfi, slope, k = k_int, bs = basis, m = m) +
    ti(vfi, wind, k = k_int, bs = basis, m = m) +
    ti(vfi, steps, k = k_int, bs = basis, m = m) +

    ti(tfi, slope, k = k_int, bs = basis, m = m) +
    ti(tfi, wind, k = k_int, bs = basis, m = m) +
    ti(tfi, steps, k = k_int, bs = basis, m = m) +

    ti(slope, wind, k = k_int, bs = basis, m = m) +
    ti(slope, steps, k = k_int, bs = basis, m = m) +

    ti(wind, steps, k = k_int, bs = basis, m = m)
)


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

# Búsqueda del máximo con previas planas pero cuadradas -------------------



# first wave
ss <- sobol(n = 600, dim = n_coef)
particles <- scale_params(ss, sup)
sim1 <- similarity_simulate_parallel(particles, spread_data)
sim1$wave <- 0
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
wave_plot(w7, alpha = 0.01)
# con var_factor = 0.25 se queda en los máximos. Probemos 0.5

# n_best es nrow(sim1), using sobol, var_factor = 0.5, weight_po = 1
w8 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, n_best = nrow(sim1),
  var_factor = 0.5, weight_power = 1, support = sup, spread_data = spread_data,
  centre = "local", sobol = TRUE
)
wave_plot(w8, alpha = 0.01)
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
wave_plot(w17, alpha = 0.1)


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



# Usando lógica de SMC para reproducir las partículas ---------------------

# El tema de que se acumulen en ciertas zonas tiene que ver mucho con el kernel
# que se usa. Pero algunas zonas tienen muchas partículas por culpa del kernel,
# y no por culpa de que ahí se acumularon más por ser zonas de alta likelihood.

# Entonces, para reproducir partículas, debería usarse el weight así:
# prior * overlap / kernel_prob,
# donde kernel_prob es la probabilidad de que la partícula haya salido
# del set de partículas anteriores, sumando las probs.

# first wave
ss <- sobol(n = 1000, dim = n_coef)
particles <- scale_params(ss, sup)
colnames(particles) <- par_names
sim1 <- similarity_simulate_parallel(particles, spread_data)
colnames(sim1$par_values)

# prior at logit scale
particles_raw <- logit_scaled(particles, sup)

sim1$prior <- logistic_prior(particles_raw) # at logit scale!!
sim1$kernel <- 1
sim1$weight <- sim1$overlap * sim1$prior    ## prior * overlap / 1
sim1$wave <- 0



# single_wave
w1 <- explore_likelihood_kernel(sim1, n = 500, p_best = 0.15, support = sup,
                                spread_data = spread_data)
wave_plot(w1)


w1 <- explore_likelihood_kernel_iterate(
  n_waves = 10, data = sim1, n = 500, n_best = 400, support = sup,
  spread_data = spread_data
)
wave_plot(w1, alpha = 0.05)


w2 <- explore_likelihood_kernel_iterate(
  n_waves = 10, data = w1, n = 500, n_best = 400, support = sup,
  spread_data = spread_data
)
wave_plot(w2, alpha = 0.01)

# no funciona, se acumulan en los extremos.

# pseudolikelihood model (GAM) ---------------------------------------------


w17 <- readRDS("files/pseudolikelihood_estimation/likelihood_points_2015_58.rds")
colnames(w17$par_values) <- par_names
# y ahora le ajustamos un GAM

# prepare data for gam
data_gam <- as.data.frame(cbind(w17$par_values,
                                y = w17$overlap)) # change logit or not here, name it "y"

# set gam parameters (to be evaluated by gam_formula)
k_side <- 30; k_int <- 5; m <- 2; basis <- "cr"

# # modelos logit-normal
# system.time(
#     mgam <- gam(gam_formula,
#               data = data_gam, method = "REML")
# )
# # usa re poca memoria, no sería un problema la RAM
# # tarda mucho. Me cansó.

# betar
mb1 <- microbenchmark(
  model_fit = {
    mgam_fast_beta <- bam(gam_formula, family = betar(),
                          data = data_gam, method = "fREML",
                          discrete = TRUE, nthreads = n_cores)
  }, times = 1
)

# logit scale
data_gam_logit <- as.data.frame(cbind(w17$par_values,
                                      y = w17$overlap_logit))
data_gam_log <- as.data.frame(cbind(w17$par_values,
                                    y = w17$overlap_log))
mb2 <- microbenchmark(
  model_logit = {
    mgam_fast_logit <- bam(gam_formula,
                           data = data_gam_logit, method = "fREML",
                           discrete = TRUE, nthreads = n_cores)
  },
  model_log = {
    mgam_fast_log <- bam(gam_formula,
                           data = data_gam_log, method = "fREML",
                           discrete = TRUE, nthreads = n_cores)
  },
  times = 1
)
mb2

# sd(residuals(mgam_fast_log, type = "response"))
sd_log <- sigma(mgam_fast_log)
sd_logit <- sigma(mgam_fast_logit)

ff <- fitted(mgam_fast_logit)
fitted_ov_invlogit_j <- numeric(nrow(data_gam_logit))
for(i in 1:length(fitted_ov_invlogit_j)) {
  fitted_ov_invlogit_j[i] <- momentsLogitnorm(ff[i], sd_logit)["mean"]
}

dcomp <- data.frame(
  overlap = w17$overlap,
  overlap_logit = w17$overlap_logit,
  overlap_log = w17$overlap_log,

  fitted_ov_invlogit = plogis(fitted(mgam_fast_logit)),
  fitted_ov_invlogit_j = fitted_ov_invlogit_j,
  fitted_ov_logit = fitted(mgam_fast_logit),

  fitted_ov_exp = exp(fitted(mgam_fast_log)),
  fitted_ov_exp_j = exp(fitted(mgam_fast_log) + 0.5 * sd_log ^ 2),
  fitted_ov_log = fitted(mgam_fast_log)
)

aa <- 0.01
par(mfrow = c(3, 2))
plot(fitted_ov_invlogit ~ overlap, dcomp, col = rgb(0, 0, 0, aa), pch = 19); abline(0, 1, lwd = 1.5, col = "red")
plot(fitted_ov_invlogit_j ~ overlap, dcomp, col = rgb(0, 0, 0, aa), pch = 19); abline(0, 1, lwd = 1.5, col = "red")
plot(fitted_ov_logit ~ overlap_logit, dcomp, col = rgb(0, 0, 0, aa), pch = 19); abline(0, 1, lwd = 1.5, col = "red")
plot(fitted_ov_exp ~ overlap, dcomp, col = rgb(0, 0, 0, aa), pch = 19); abline(0, 1, lwd = 1.5, col = "red")
plot(fitted_ov_log ~ overlap_log, dcomp, col = rgb(0, 0, 0, aa), pch = 19); abline(0, 1, lwd = 1.5, col = "red")
plot(fitted_ov_exp_j ~ overlap, dcomp, col = rgb(0, 0, 0, aa), pch = 19); abline(0, 1, lwd = 1.5, col = "red")
par(mfrow = c(1, 1))


mgam_fast_logit$sp


# add fitted values and residuals
w20$fitted <- plogis(fitted(mgam))
w20$fitted_logit <- fitted(mgam)
w20$res <- w20$fitted - w20$overlap
w20$res_logit <- w20$fitted_logit - w20$overlap_logit

# # for beta:
# w17$fitted <- fitted(mgam_fast)
# w17$res <- w17$fitted - w17$overlap

wave_plot(w17, "overlap")
wave_plot(w17, "fitted")
plot(fitted ~ overlap, w17, col = rgb(0, 0, 0, 0.2), pch = 19)
abline(0, 1, lwd = 1.5, col = "red")
# no anduvo usando te(), pero sí con te()

ccc <- cov2cor(vcov(mgam_fast))
range(ccc[lower.tri(ccc)]) # muchos par correlated

wave_plot(w17[w17$par_values[, 6] > 250, ], alpha = 0.3)
wave_plot(w17[w17$par_values[, 6] > 250, ], response = "fitted", alpha = 0.3)




# Pruebitas con gam -------------------------------------------------------

gam_formula_te <- formula(
  y ~ # te interactions
    te(intercept, vfi, k = k_int_te, bs = basis, m = m) +
    te(intercept, tfi, k = k_int_te, bs = basis, m = m) +
    te(intercept, slope, k = k_int_te, bs = basis, m = m) +
    te(intercept, wind, k = k_int_te, bs = basis, m = m) +
    te(intercept, steps, k = k_int_te, bs = basis, m = m) +

    te(vfi, tfi, k = k_int_te, bs = basis, m = m) +
    te(vfi, slope, k = k_int_te, bs = basis, m = m) +
    te(vfi, wind, k = k_int_te, bs = basis, m = m) +
    te(vfi, steps, k = k_int_te, bs = basis, m = m) +

    te(tfi, slope, k = k_int_te, bs = basis, m = m) +
    te(tfi, wind, k = k_int_te, bs = basis, m = m) +
    te(tfi, steps, k = k_int_te, bs = basis, m = m) +

    te(slope, wind, k = k_int_te, bs = basis, m = m) +
    te(slope, steps, k = k_int_te, bs = basis, m = m) +

    te(wind, steps, k = k_int_te, bs = basis, m = m)
)
# cambiar
mb3 <- microbenchmark(
  model_logit = {
    k_int_te <- 6
    mgam_fast_logit_te <- bam(gam_formula_te,
                           data = data_gam_logit, method = "fREML",
                           discrete = TRUE, nthreads = n_cores)
  },
  times = 1
)
mb3
summary(mgam_fast_logit_te)
plot(mgam_fast_logit_te)
plot(plogis(fitted(mgam_fast_logit_te)) ~ w17$overlap)
length(coef(mgam_fast_logit_te))
names(coef(mgam_fast_logit_te))

names(coef(mgam_fast_logit)[-1])

k_side <- 31
klist <- lapply(1:n_coef, function(i) {
  seq(sup[1, i], sup[2, i], length.out = k_side)
})
names(klist) <- par_names

k_int <- 11 # 5 * 5
# m_logit_knots <- bam(
#   gam_formula, data = data_gam_logit, method = "fREML", discrete = TRUE,
#   nthreads = n_cores, knots = klist
# )


m_logit_more <- bam(
  gam_formula, data = data_gam_logit, method = "fREML", discrete = TRUE,
  nthreads = n_cores, gamma = 0.01#, knots = klist
)
summary(m_logit_more)

plot(plogis(fitted(m_logit_more)) ~ w17$overlap,
     col = rgb(0, 0, 0, 0.01), pch = 19)

plot(plogis(fitted(m_logit_more)) - w17$overlap ~ w17$overlap,
     col = rgb(0, 0, 0, 0.05), pch = 19)

# Estudiar cómo proveer knots a ambas partes y luego
(10^2)*15


# con 10 bases por lado en las interacciones anda bien.
#
m_logit_more$smooth[[1]] %>% str # intercept
m_logit_more$smooth[[7]] %>% str # intercept,vfi
m_logit_more$smooth[[1]] %>% str # intercept

plot(m_logit_more$smooth[[1]]$xp)
plot(m_logit_more$smooth[[6]]$xp)

plot(m_logit_more$smooth[[11]]$xp)

# Igual parece que los knots están re bien puestos.


# hacer pred parciales, comparar posteriores, y bla.


# Partial predictions -----------------------------------------------------

k_side <- 30 +1
k_int <- 10 + 1 # 5 * 5
m = 1; basis = "cr"

models_list <- lapply(c(1e-12, 0.01, 0.5, 1), function(g) {
  m <- bam(
    gam_formula, data = data_gam_logit, method = "fREML", discrete = TRUE,
    nthreads = n_cores, gamma = g
  )

  plot(plogis(fitted(m)) ~ w17$overlap, col = rgb(0, 0, 0, 0.01), pch = 19,
       main = g)

  return(m)
})



# Uniform prior at logit scale? -------------------------------------------

rr <- runif(1e6)
qq <- qlogis(rr)
plot(density(qq))
curve(dlogis(x), add = T, col = "red")
curve(dnorm(x), add = T, col = "blue")
# flat prior in [0, 1] implies logistic distribution at the logit scale.

# this implies we can compute the prior at the flat prior transformed
# to the logit scale :)

qq <- rlogis(1e6)
rr <- plogis(qq)
plot(density(rr, from = 0, to = 1))
plot(quantile(rr, ppoints(100)) ~ ppoints(100)); abline(0, 1)



# TAREA -------------------------------------------------------------------

# Testear la exploración con funciones conocidas.
# Alternar búsqueda global con local, por si hay modos pequeños y finitos.
# Tunear la seq de varianzas locales.
# Ver cómo combinar local-global.
# Agrandar límites a [-50, 50].

