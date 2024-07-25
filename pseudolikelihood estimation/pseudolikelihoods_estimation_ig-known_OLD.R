# This code inherits from
# <single fire posterior - surrogate vs. simulator - extended.R>.

# There I tried many ways to fit a good surrogate model, and also tried to
# get a surrogate for a kernel-transformed similarity. But that did not work.

# The current solution is to define an overlap threshold and fit a binary GAM
# to predict where the overlap is above the threshold. To fit the mixed model,
# we would sample the probability, not simulating the Bernoulli event.

# The current problem is that the surrogate has high uncertainty in some regions
# of the parameter space, generating bias in the resulting posterior.
# Here I simulate more data and refit the model a few times to get reasonable
# uncertainty where it is too high (model-based design).

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(ggdensity) # geom_hdr

library(Rcpp)
library(terra)
library(abind)         # manage arrays

library(randtoolbox)   # sobol sequences
library(mgcv)          # fit GAM
library(MGMM)          # fit MVN distribution
library(Matrix)        # nearPD, searches the closest PD matrix.
library(stutils)       # tryNULL, for search of nearPD
library(logitnorm)

library(sn)
library(trialr) # rlkjcorr
library(truncnorm)

library(adaptMCMC)
library(posterior)     # manage posterior samples
library(tidybayes)     # not sure if this was used
library(bayesplot)     # visualize posteriors

library(foreach)       # parallelization
library(doMC)          # doRNG had problems with cpp functions

library(FireSpread)    # spread and similarity functions

library(microbenchmark)

source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for rast_from_mat and a few constants

# source("estimation_functions.R") # prior_dist and other stuff

# Constants --------------------------------------------------------------

data_dir <- file.path("data", "focal fires data", "landscapes_ig-known")
filenames <- list.files(data_dir)

# dir to save output
target_dir <- file.path("files", "pseudolikelihoods")

# constants for fire spread simulation
upper_limit <- 1
n_coef <- 6
par_names <- c("intercept", "vfi", "tfi", "slope", "wind", "steps")

# parameters bounds
ext_alpha <- 30
ext_beta <- 30

# number of fires to simulate by particle
n_sim <- 20

# gam to approximate likelihood

k_side <- 11
k_side_more <- 11
k_int <- 6; k_int_less <- 5
basis <- "cr"
fixed <- F

gam_formula <- formula(
  y ~

    # marginal effects
    s(intercept, bs = basis, k = k_side, fx = fixed) +
    s(vfi, bs = basis, k = k_side, fx = fixed) +
    s(tfi, bs = basis, k = k_side, fx = fixed) +
    s(slope, bs = basis, k = k_side, fx = fixed) +
    s(wind, bs = basis, k = k_side, fx = fixed) +
    s(steps, bs = basis, k = k_side_more, fx = fixed) +

    # interactions
    ti(intercept, vfi, k = k_int, bs = basis, fx = fixed) +
    ti(intercept, tfi, k = k_int, bs = basis, fx = fixed) +
    ti(intercept, slope, k = k_int, bs = basis, fx = fixed) +
    ti(intercept, wind, k = k_int, bs = basis, fx = fixed) +
    ti(intercept, steps, k = k_int, bs = basis, fx = fixed) +

    ti(vfi, tfi, k = k_int, bs = basis, fx = fixed) +
    ti(vfi, slope, k = k_int, bs = basis, fx = fixed) +
    ti(vfi, wind, k = k_int, bs = basis, fx = fixed) +
    ti(vfi, steps, k = k_int, bs = basis, fx = fixed) +

    ti(tfi, slope, k = k_int, bs = basis, fx = fixed) +
    ti(tfi, wind, k = k_int, bs = basis, fx = fixed) +
    ti(tfi, steps, k = k_int, bs = basis, fx = fixed) +

    ti(slope, wind, k = k_int, bs = basis, fx = fixed) +
    ti(slope, steps, k = k_int, bs = basis, fx = fixed) +

    ti(wind, steps, k = k_int, bs = basis, fx = fixed)
)

# smaller k for steps interactions
gam_formula_simpler <- formula(
  y ~

    # marginal effects
    s(intercept, bs = basis, k = k_side, fx = fixed) +
    s(vfi, bs = basis, k = k_side, fx = fixed) +
    s(tfi, bs = basis, k = k_side, fx = fixed) +
    s(slope, bs = basis, k = k_side, fx = fixed) +
    s(wind, bs = basis, k = k_side, fx = fixed) +
    s(steps, bs = basis, k = k_side_more, fx = fixed) +

    # interactions
    ti(intercept, vfi, k = k_int, bs = basis, fx = fixed) +
    ti(intercept, tfi, k = k_int, bs = basis, fx = fixed) +
    ti(intercept, slope, k = k_int, bs = basis, fx = fixed) +
    ti(intercept, wind, k = k_int, bs = basis, fx = fixed) +
    ti(intercept, steps, k = k_int_less, bs = basis, fx = fixed) +

    ti(vfi, tfi, k = k_int, bs = basis, fx = fixed) +
    ti(vfi, slope, k = k_int, bs = basis, fx = fixed) +
    ti(vfi, wind, k = k_int, bs = basis, fx = fixed) +
    ti(vfi, steps, k = k_int_less, bs = basis, fx = fixed) +

    ti(tfi, slope, k = k_int, bs = basis, fx = fixed) +
    ti(tfi, wind, k = k_int, bs = basis, fx = fixed) +
    ti(tfi, steps, k = k_int_less, bs = basis, fx = fixed) +

    ti(slope, wind, k = k_int, bs = basis, fx = fixed) +
    ti(slope, steps, k = k_int_less, bs = basis, fx = fixed) +

    ti(wind, steps, k = k_int_less, bs = basis, fx = fixed)
)


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
similarity_simulate_simpler <- function(particle, n_sim = 20,
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

wave_plot <- function(data, response = "overlap", alpha = 0.3, best = NULL,
                      x = "par_values", thres = NULL, bin = FALSE) {

  if(!is.null(best)) {
    data <- data[order(data[, response], decreasing = TRUE), ]
  }

  if(bin) {
    data$bin <- jitter(as.numeric(data$overlap >= thres), factor = 0.2)
    response <- "bin"
  }

  yy <- range(data[, response])
  yy[2] <- yy[2] * 1.05

  # title for first plot
  tit <- deparse(substitute(data))

  par(mfrow = c(3, 2))
  for(i in 1:n_coef) {
    mm <- ifelse(i == 1, tit, NA)
    plot(data[, response] ~ data[, x][, i], ylab = response,
         xlab = par_names[i], ylim = yy, main = mm,
         pch = 19, col = rgb(0, 0, 0, alpha))

    if(!is.null(best)) {
      points(data[1:best, response] ~ data[, x][1:best, i],
             pch = 19, col = rgb(0, 0, 1, alpha))
    }

    if(!is.null(thres) & response == "overlap") {
      abline(h = thres, col = "red")
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

invlogit_scaled_vec <- function(xun, support) {
  dmax <- support[2, ] - support[1, ]
  x <- plogis(xun) * dmax + support[1, ]
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

logit_scaled_vec <- function(x, support) {
  dobs <- x - support[1, ]
  dmax <- support[2, ] - support[1, ]
  x01 <- dobs / dmax
  return(qlogis(x01))
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
# accept_thres: overlap threshold below which particles will not be resampled.
#   overrides both n_best and p_best.
# prob_fn: function to transform overlap into probability. If not NULL,
#   the probability is used to resample particles.
# support: matrix with lower limits (first row) and upper limits (second row)
#   of all variables.
# spread_data: data to simulate fires.
# centre: "local" or "global". New particles are simulated from a MVN centred
#   at each resampled particle (focal), or at the global mean (global). If global,
#   sobol may be set to TRUE, if focal, sobol is set to FALSE.
explore_likelihood <- function(data, n = 600, var_factor = 1, n_best = "all",
                               p_best = 0.1, accept_thres = NULL,
                               prob_fn = NULL,
                               support, spread_data,
                               centre = "local", sobol = TRUE,
                               phase = NULL) {

  ### Testo
  # data = sim1; n = 500; p_best = 0.15; prob_fn = NULL;
  # support = sup; spread_data = spread_data;
  # centre = "global"; sobol = TRUE;
  # pow = 1; var_factor = 1
  ### endo testo

  # if threshold is used, ignore p_best
  if(!is.null(accept_thres)) {
    p_best <- NULL
    n_best <- sum(data$overlap >= accept_thres)
  }

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
  weight <- data$overlap

  # sometimes the best particles are too few, and the computation of the vcov
  # fails. In these cases we need to smooth a bit the weights, so more than 20
  # different particles are resampled.
  l <- 0
  k <- 0 # smoothing power (root)
  while(l < 30) {
    ids_rep <- sample(1:n_best, size = n, replace = T,
                      prob = weight[1:n_best] ^ (1 / (1 + k)))
    l <- length(unique(ids_rep))
    k <- k + 0.5
  }

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
  colnames(candidates) <- par_names

  # message("Simulating fires")
  sim_result <- similarity_simulate_parallel(particles_mat = candidates,
                                             fire_data = spread_data)

  # define wave
  sim_result$wave <- max(data$wave) + 1
  sim_result$phase <- phase
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
    data, n = 500, var_factor = 1.5,
    n_best = "all", p_best = NULL, accept_thres = NULL,
    support, spread_data,
    centre = "global", sobol = TRUE,
    phase = NULL,
    write = FALSE, fire_name = NULL
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
                              accept_thres = accept_thres,
                              support = support, spread_data = spread_data,
                              centre = centre, sobol = sobol,
                              phase = phase)

    m <- paste("Search wave ", last_wave + i, "; max overlap = ",
               round(max(res$overlap), 4), sep = "")

    if(write) {
      nn <- paste(fire_name, "-wave-", last_wave + i, ".rds", sep = "")
      saveRDS(res[res$wave == last_wave + i, ], file.path(target_dir, nn))
      # save only latest
    }

    message(m)
  }

  return(res)
}

# Function to perform adaptive MCMC in parallel by hand. It's just a wrapper
# for adaptMCMC::MCMC, because MCMC.parallel does not allow to fix the inits.
MCMC_parallel <- function(fun, n, adapt, scale, init_list, acc.rate = 0.234,
                          n_cores = 8, ...) {

  registerDoMC(n_cores)

  runs <- foreach(pp = init_list) %dopar% {
    MCMC(p = fun, n = n, adapt = adapt, scale = scale, init = pp,
         acc.rate = acc.rate, ...)
  }

  return(runs)
}


# make prob = 0 in particles with high uncertainty, measured as the width of the
# 95 % confidence interval (unc).
drop_uncertain <- function(unc, p, a = 73, b = 90) {
  return(p * plogis(a - b * unc))
}
curve(drop_uncertain(x, 1, 20, 90), ylim = c(0, 1))


# sample from the GAM taking independent draws (simulation)
# unc_thres [0, 1] determines the maximum uncertainty allowed in the model.
#   if higher, the particle is rejected.
rejection_sample <- function(iter, model, support) {

  ntry <- iter * 10
  got <- 0
  wave <- 0
  samples <- matrix(NA, ntry, n_coef)
  colnames(samples) <- par_names

  while(got < iter) {
    unif_draws <- matrix(runif(ntry * n_coef), ntry, n_coef)
    particles <- scale_params(unif_draws, support)
    colnames(particles) <- par_names
    part_df <- as.data.frame(particles)

    # compute probability
    prob <- predict(gam_bern, part_df, type = "response")

    # only evaluate uncertainty at particles with high probability.
    # (This step is time consuming)
    unc <- numeric(length(prob))
    high_prob_ids <- which(prob > 0.02)
    if(length(high_prob_ids) > 0) {
      pp <- predict(gam_bern, part_df[high_prob_ids, ], se.fit = T)
      lll <- plogis(pp$fit - qnorm(0.975) * pp$se.fit)
      uuu <- plogis(pp$fit + qnorm(0.975) * pp$se.fit)
      unc[high_prob_ids] <- uuu - lll
    }

    # drop probability if highly uncertain
    prob2 <- drop_uncertain(unc, prob)
    stay <- which(runif(ntry) <= prob2)

    l <- length(stay)
    if(l > 0) {
      samples[(got+1):(got+l), ] <- particles[stay, ]
    }

    wave <- wave + 1
    got <- got + l
    # print(paste("wave ", wave, ", got ", got, sep = ""))
  }

  samples <- samples[1:iter, ]
  return(samples)
}

# the same in parallel
rejection_sample_parallel <- function(iter, model, support, cores = 8) {
  registerDoMC(cores)
  ii <- as.list(1:cores)
  foreach(chain = ii) %dopar% {
    rejection_sample(iter, model, support)
  }
}


# extract samples from a MCMC.parallel run and return a draws_array
# (posterior class). (Also useful for MCMC_parallel result.)
tidy_samples <- function(samples, adapt, support) {
  nc <- length(samples)
  rr <- nrow(samples[[1]]$samples) - adapt
  cc <- ncol(samples[[1]]$samples)
  arr0 <- array(NA, dim = c(rr, cc, nc))
  for(c in 1:nc) {
    arr0[, , c] <- invlogit_scaled(samples[[c]]$samples[-(1:adapt), ], support)
  }
  arr <- aperm(arr0, c(1, 3, 2))
  draws_arr <- as_draws_array(arr)

  draws_arr2 <- rename_variables(draws_arr,
                   "intercept" = "...1",
                   "vfi" = "...2",
                   "tfi" = "...3",
                   "slope" = "...4",
                   "wind" = "...5",
                   "steps" = "...6")

  return(draws_arr2)
}

# the same as tidy samples to be used with the result from
# rejection_sample_parallel
tidy_samples_ind <- function(samples_list) {
  nc <- length(samples_list)
  rr <- nrow(samples_list[[1]])
  cc <- ncol(samples_list[[1]])
  arr0 <- array(NA, dim = c(rr, cc, nc))
  for(c in 1:nc) {
    arr0[, , c] <- samples_list[[c]]
  }
  arr <- aperm(arr0, c(1, 3, 2))
  draws_arr <- as_draws_array(arr)

  draws_arr2 <- rename_variables(draws_arr,
                                 "intercept" = "...1",
                                 "vfi" = "...2",
                                 "tfi" = "...3",
                                 "slope" = "...4",
                                 "wind" = "...5",
                                 "steps" = "...6")

  return(draws_arr2)
}

# Function to get a smaller subset of points to fit the GP. It returns a vector
# with the ids of the reduced points to use. Optionally, the best B points
# may be kept, performing the reduction only at lower values.
# The function performs a kmeans clustering and finds the point closest to the
# centroid in each cluster.
# The support is needed to transform the variables into the [0, 1] scale, so
# all variables contribute equally to the distance computation.
spread_points <- function(data, K = 2000, B = NULL, support, ...) {

  # make row_id to keep track of everything
  data$row_id <- 1:nrow(data)

  # keep the best B or not
  if(!is.null(B)) {
    data <- data[order(data$overlap, decreasing = T), ]
    best_ids <- data$row_id[1:B]
    data_k <- data[(B+1):nrow(data), ]
    K <- K - B
  } else {
    data_k <- data
    best_ids <- NULL # initialize at null
  }

  # data is scaled to [0, 1] so all distances weight the same
  data_k$par_values <- plogis(logit_scaled(data_k$par_values, support))
  kk <- kmeans(data_k$par_values, K, ...)

  # nearest points from centroids:
  closest_ids <- sapply(1:K, function(i) {
    rows_local <- data_k$row_id[kk$cluster == i]
    neighs <- data_k$par_values[data_k$row_id %in% rows_local, , drop = F]
    n <- nrow(neighs)
    dd <- dist(rbind(kk$centers[i, , drop = F], neighs)) %>% as.numeric
    id_near <- which.min(dd[1:n])
    return(rows_local[id_near])
  })

  res <- c(best_ids, closest_ids)
  return(res)
}

prob_fun <- function(x, sigma = 0.15, pow = 2, rescale = FALSE, ov_max = NULL) {
  if(!rescale) { # use a standard kernel
    p <- exp(-((1 - x) / sigma) ^ pow)
    return(p)
  } else { # use fire-wise kernel
    if(is.null(ov_max)) {
      ov_max <- max(x)
    }

    # recale overlap to [0, 1]
    x01 <- x / ov_max

    # turn into 1 the values larger than ov_max
    x01[x01 > 1] <- 1

    # turn into prob
    p <- exp(-((1 - x01) / sigma) ^ pow)
    return(p)
  }
}

# function to evaluate ABC-posterior, based on an overlap threshold
fn_like_sim_bin <- function(x, support = NULL, fire_data = NULL, thres = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  ov <- similarity_simulate_simpler(xx, n_sim = 20, fire_data = fire_data)
  if(ov >= thres) {
    return(sum(dlogis(x, log = TRUE)))
  } else {
    return(-Inf)
  }
}

# function to reduce support, so the rejection sampling is faster.
reduce_support <- function(x, support, prop = 0.75) {
  ranges <- apply(x, 2, range)
  widths <- apply(ranges, 2, diff)

  rnew <- ranges

  rnew[1, ] <- ranges[1, ] - widths * prop
  rnew[2, ] <- ranges[2, ] + widths * prop

  # constrain to the original support
  rnew[1, ] <- ifelse(rnew[1, ] < support[1, ], support[1, ], rnew[1, ])
  rnew[2, ] <- ifelse(rnew[2, ] > support[2, ], support[2, ], rnew[2, ])

  return(rnew)
}


# Loop over fires -------------------------------------------------------

# tests
# f <- 7# 20

for(f in 1:length(filenames)) {

  fire_file <- filenames[f]
  fire_name <- strsplit(fire_file, ".rds")[[1]]
  print(paste("Fire:", fire_name))
  full_data <- readRDS(file.path(data_dir, fire_file))
  # subset data needed for spread (to be cloned across workers)
  spread_data <- full_data[c("landscape", "burnable", "ig_rowcol",
                             "burned_layer", "burned_ids")]

  bb <- get_bounds(fire_data = spread_data)

  params_lower <- c("intercept" = -ext_alpha,
                    "vfi" = 0,
                    "tfi" = 0,
                    "slope" = 0,
                    "wind" = 0,
                    "steps" = 5)
  params_upper <- c("intercept" = ext_alpha,
                    "vfi" = ext_beta,
                    "tfi" = ext_beta,
                    "slope" = ext_beta * 2,
                    "wind" = ext_beta,
                    "steps" = bb["largest", "steps_used"])
  sup <- rbind(params_lower, params_upper)


  # Explore likelihood
  # Phase 1: sobol
  registerDoMC(15)

  ss <- sobol(n = 10000, dim = n_coef)
  particles <- scale_params(ss, sup)
  waves_1 <- rep(1:20, each = 500)
  # loop to save
  wave1 <- NULL
  print(paste("Phase: sobol. Fire: ", fire_name, sep = ""))
  for(w in 1:20) {
    # w = 1
    wlocal <- similarity_simulate_parallel(particles[waves_1 == w, ],
                                           spread_data)
    wlocal$wave <- w
    wlocal$phase <- "sobol"
    colnames(wlocal$par_values) <- par_names

    nn <- paste(fire_name, "-wave-", w, ".rds", sep = "")
    saveRDS(wlocal, file.path(target_dir, nn))

    wave1 <- rbind(wave1, wlocal)
    rr <- round(max(wave1$overlap), 4)
    print(paste("wave ", w, ", overlap max: ", rr, sep = ""))
  }

  # wave 2: reproduce 10 % best (500 * 10)
  print(paste("Phase: search higher. Fire: ", fire_name, sep = ""))
  wave2 <- explore_likelihood_iterate(
    n_waves = 10, data = wave1, n = 500, p_best = 0.10,
    var_factor = c(1.5, 2, 1),
    support = sup, spread_data = spread_data,
    centre = "global", sobol = TRUE, phase = "search_higher",
    write = TRUE, fire_name = fire_name
  )

  # define threshold
  n_accept <- 100
  pcum_thres <- 1 - n_accept / nrow(wave2)
  thres <- quantile(wave2$overlap, pcum_thres) %>% unname

  # wave 3: 5000 above threshold
  print(paste("Phase: above threshold. Fire: ", fire_name, sep = ""))
  wave3 <- explore_likelihood_iterate(
    n_waves = 10, data = wave2, n = 500, accept_thres = thres,
    var_factor = c(1.5, 2, 1),
    support = sup, spread_data = spread_data,
    centre = "global", sobol = TRUE, phase = "above_threshold",
    write = TRUE, fire_name = fire_name,
  )

  like_sim <- wave3
  attr(like_sim, "threshold") <- thres
  nn <- paste(fire_name, "simulations_data.rds", sep = "-")
  saveRDS(like_sim, file.path(target_dir, nn))

  # remove intermediate files
  all_files_saved <- list.files(target_dir)
  remove_this <- all_files_saved[grep(paste(fire_name, "-wave", sep = ""), all_files_saved)]
  invisible(lapply(remove_this, function(f) unlink(file.path(target_dir, f))))

  # Binary GAM
  data_gam_bern <- as.data.frame(cbind(like_sim$par_values,
                                       y = as.numeric(like_sim$overlap >= thres)))
  print(paste("Fitting GAM. Fire: ", fire_name, sep = ""))
  gam_bern <- bam(
    gam_formula, data = data_gam_bern, family = "binomial",
    method = "fREML", discrete = T, nthreads = 8
  )
  nn <- paste(fire_name, "likelihood_model_gam.rds", sep = "-")
  saveRDS(gam_bern, file.path(target_dir, nn))

  # Sample GAM-approximated posterior

  print(paste("Sampling posterior. Fire: ", fire_name, sep = ""))
  apply(like_sim$par_values[like_sim$overlap >= thres, ], 2, range)
  sup_reduced <- reduce_support(like_sim$par_values[like_sim$overlap >= thres, ],
                                sup, 20)
  r_gam <- rejection_sample_parallel(300, gam_bern, sup_reduced, cores = 15)

  draws <- do.call("rbind", r_gam) %>% as.data.frame

  nn <- paste(fire_name, "posterior_samples.rds", sep = "-")
  saveRDS(r_gam, file.path(target_dir, nn))

  # Plot posterior
  dflong <- pivot_longer(draws, all_of(1:n_coef), values_to = "par_value",
                         names_to = "par_names")
  dflong$par_names <- factor(dflong$par_names, levels = par_names)

  marginals_plot <-
    ggplot(dflong, aes(x = par_value)) +
      geom_density(alpha = 0, adjust = 1.5) +
      facet_wrap(vars(par_names), scales = "free", strip.position = "bottom") +
      theme(panel.grid.minor = element_blank(),
            strip.background = element_rect(color = "white", fill = "white"),
            strip.placement = "outside",
            axis.title.x = element_blank(),
            strip.text.x = element_text(vjust = 3)) +
      ggtitle(fire_name)
  ggsave(file.path(target_dir, paste(fire_name, "-plot_marginals.png", sep = "")),
         plot = marginals_plot, height = 12, width = 15, units = "cm")



  # posteriores de a pares

  combs <- as.data.frame(combn(par_names, 2) %>% t)
  names(combs) <- c("x", "y")
  ncomb <- nrow(combs)
  plist <- vector("list", ncomb)

  for(i in 1:ncomb) {
    # i = 1
    take <- combs[i, ] %>% as.character()
    dlocal <- draws[, take]

    p <-
    ggplot(dlocal, aes_string(x = take[1], y = take[2])) +
      geom_hdr(probs = seq(0.05, 0.95, by = 0.1),
               method = method_kde(adjust = 1.5)) +
      theme(panel.grid.minor = element_blank(),
            legend.position = "none")

    if(i == 9) {
      p <- p + theme(legend.position = "right")
    }

    if(i == 1) {
      p <- p + ggtitle(fire_name)
    }

    plist[[i]] <- p
  }

  pairs_plot <- egg::ggarrange(plots = plist, ncol = 3, draw = FALSE,
                               newpage = FALSE)
  ggsave(file.path(target_dir, paste(fire_name, "-plot_pairs.png", sep = "")),
         plot = pairs_plot, height = 18, width = 15, units = "cm")
  dev.off()

  # clean
  remove(spread_data, gam_bern, r_gam, like_sim, wave3)
  gc()
}

# Perfecto. Tardó 1.8 días aprox.
# Varias posteriores son bimodales; quizás convenga simular más algunos fuegos.
# O quizás ahora podría dedicarme a aceptar sólo las partículas mayores a
# max - 0.05



# Checking problematic posteriors -----------------------------------------

# Some posteriors are bimodal, which is suspectful. Simulate more points
# for those fires.

names_check <- c(
  "1999_28",
  "1999_2140469994_r",
  "2002_36",
  "2004_23",
  "2005_9",
  "2006_20E",
  "2009_3",
  "2009_28",
  "2009_2007583421",
  "2012_58",
  "2014_3",
  "2015_16", # no está tan mal
  "2015_42"
)

fnames <- list.files(target_dir)
dnames <- fnames[grep("simulations_data", fnames)] # data names
firenames_all <- sapply(dnames, function(n) {
  bb <- strsplit(n, split = "-simulations_data.rds")[[1]]
  return(bb)
})

names(dnames) <- unname(firenames_all)
dnames_check <- dnames[names_check]

dlist_check <- lapply(dnames_check, function(f) {
  readRDS(file.path(target_dir, f))
})


# Editting indivitual fires ------------------------------------------------

fire_name <- "2015_42" # vary this to check everyone manually
fire_file <- paste(fire_name, ".rds", sep = "")
full_data <- readRDS(file.path(data_dir, fire_file))

# subset data needed for spread (to be cloned across workers)
spread_data <- full_data[c("landscape", "burnable", "ig_rowcol",
                           "burned_layer", "burned_ids")]

bb <- get_bounds(fire_data = spread_data)

params_lower <- c("intercept" = -ext_alpha,
                  "vfi" = 0,
                  "tfi" = 0,
                  "slope" = 0,
                  "wind" = 0,
                  "steps" = 5)
params_upper <- c("intercept" = ext_alpha,
                  "vfi" = ext_beta,
                  "tfi" = ext_beta,
                  "slope" = ext_beta * 2,
                  "wind" = ext_beta,
                  "steps" = bb["largest", "steps_used"])
sup <- rbind(params_lower, params_upper)

like_sim <- dlist_check[[fire_name]]#readRDS(file.path(target_dir, "2002_36-simulations_data_extra.rds"))#
thres <- attr(like_sim, "threshold")

x11()
wave_plot(like_sim, thres = thres)
dev.off()

# # add more data (for some fires)
registerDoMC(15)
like_sim <- explore_likelihood_iterate(
  n_waves = 10, data = like_sim,
  n = 500, n_best = 250, #accept_thres = thres,
  var_factor = c(1.5, 2, 1),
  support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE, phase = "more_data",
  write = TRUE, fire_name = fire_name,
)
x11()
wave_plot(like_sim, thres = 0.55)
dev.off()
# more data with higher threshold
thres <- 0.55
like_sim <- explore_likelihood_iterate(
  n_waves = 15, data = like_sim,
  n = 500, accept_thres = thres,
  var_factor = c(1.5, 2, 1),
  support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE, phase = "above_threshold_2",
  write = TRUE, fire_name = fire_name,
)


# fit a simpler gam
data_gam_bern <- as.data.frame(cbind(like_sim$par_values,
                                     y = as.numeric(like_sim$overlap >= thres)))
gam_simpler <- bam(
  gam_formula_simpler, data = data_gam_bern, family = "binomial",
  method = "fREML", discrete = T, nthreads = 8
)

# # (just a normal gam for some fires, but keep the name)
gam_simpler <- bam(
  gam_formula, data = data_gam_bern, family = "binomial",
  method = "fREML", discrete = T, nthreads = 8
)

# compute densities in another way:
# huge sobol sequence, evaluate all particles, and use the weight in the density
# computation.

nlarge <- 2e6
unif_draws <- sobol(nlarge, n_coef, init = T)
particles <- scale_params(unif_draws, support = sup)
colnames(particles) <- par_names
part_df <- as.data.frame(particles)

# compute probability
prob <- predict(gam_simpler, part_df, type = "response")

# only evaluate uncertainty at particles with high probability.
# (This step is time consuming)
unc <- numeric(length(prob))
high_prob_ids <- which(prob > 0.02)
if(length(high_prob_ids) > 0) {
  pp <- predict(gam_simpler, part_df[high_prob_ids, ], se.fit = T)
  lll <- plogis(pp$fit - qnorm(0.975) * pp$se.fit)
  uuu <- plogis(pp$fit + qnorm(0.975) * pp$se.fit)
  unc[high_prob_ids] <- uuu - lll
}

# drop probability if highly uncertain
prob2 <- drop_uncertain(unc, prob)#, 10, 90)
prob2 <- normalize(prob2)

# compute densities in each margin
dmargs <- do.call("rbind", lapply(par_names, function(p) {
  d <- density(particles[, p], weights = prob2, n = 2 ^ 11,
               from = sup[1, p], to = sup[2, p], adjust = 3)
  dd <- data.frame(x = d$x, density = d$y, var = p)
  return(dd)
}))
dmargs$var <- factor(dmargs$var, levels = par_names)

# add points at density = 0
ddghost <- data.frame(y = rep(0, n_coef),
                      x = sup[1, ],
                      var = par_names)
ddghost$var <- factor(ddghost$var, levels = par_names)

marginals_plot <-
  ggplot(dmargs, aes(x = x, y = density)) +
  geom_line() +
  geom_point(aes(x, y), ddghost, alpha = 0) +
  facet_wrap(vars(var), scales = "free", strip.position = "bottom") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(color = "white", fill = "white"),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        strip.text.x = element_text(vjust = 3)) +
  ggtitle(fire_name)
marginals_plot

# replace image
ggsave(file.path(target_dir, paste(fire_name, "-plot_marginals.png", sep = "")),
       plot = marginals_plot, height = 12, width = 15, units = "cm")

# Es súper rápido estimar las marginales así.
# Para ver las bivariadas, podría remuestrear esos puntos en base a prob2.

boot_ids <- sample(1:nrow(part_df), size = 20000, prob = prob2, replace = T)
draws <- part_df[boot_ids, ]

# posteriores de a pares

combs <- as.data.frame(combn(par_names, 2) %>% t)
names(combs) <- c("x", "y")
ncomb <- nrow(combs)
plist <- vector("list", ncomb)

for(i in 1:ncomb) {
  # i = 1
  take <- combs[i, ] %>% as.character()
  dlocal <- draws[, take]

  p <-
    ggplot(dlocal, aes_string(x = take[1], y = take[2])) +
    geom_hdr(probs = seq(0.05, 0.95, by = 0.1),
             method = method_kde(adjust = 2)) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")

  if(i == 9) {
    p <- p + theme(legend.position = "right")
  }

  if(i == 1) {
    p <- p + ggtitle(fire_name)
  }

  plist[[i]] <- p
}

pairs_plot <- egg::ggarrange(plots = plist, ncol = 3, draw = FALSE,
                             newpage = FALSE)
# replace image
ggsave(file.path(target_dir, paste(fire_name, "-plot_pairs.png", sep = "")),
       plot = pairs_plot, height = 18, width = 15, units = "cm")
dev.off()

# replace gam
nn <- paste(fire_name, "likelihood_model_gam.rds", sep = "-")
saveRDS(gam_simpler, file.path(target_dir, nn))

# replace posterior samples
nn <- paste(fire_name, "posterior_samples.rds", sep = "-")
draws <- part_df[boot_ids, ]
saveRDS(draws, file.path(target_dir, nn))

# List of changes ---------------------------------------------------------

# 1999_28             It does not improve with more data, nor with a simpler GAM.
#                     With the chosen threshold it has a flat posterior, with
#                     high correlation between wind and slope.

# 1999_2140469994_r   Not so bad

# 2002_36							It does not improve with more data, nor with a simpler GAM.

# 2004_23             Not so bad, makes sense with data. It's a complicated
#                     posterior.

# 2005_9              Not so bad, makes sense with data. It's a complicated
#                     posterior.

# 2006_20E            Bad GAM approximation. An absurd bimodality in steps was
#                     resolved using a smaller k in the te() terms involving
#                     steps, k_int_smaller = 4.

# 2009_3							Same problem as above. It was a bad GAM. Resolved using
#                     k_int_smaller = 3.

# 2009_28             This one needed more data. I searched 10000 more particles
#                     above the threshold, then increased the threshold to 0.5
#                     and searched 10000 more above this new threshold.
#                     The same gam (not simpler) was fitted to the new dataset.

# 2009_2007583421     Not so bad.

# 2012_58             Bad GAM. Resolved using a simpler one, with k_int_less = 5.

# 2014_3              Needed more data. Search 5000 particles more with n_best = 250.
#                     Then, define new threshold at 0.31, and search 7500 particles
#                     above it. Fitted the same GAM as the first time (not simpler.)

# 2015_16             Not so bad. The effect from steps is not easily identified
#                     from this fire.

# 2015_42             This one needed more data. I searched 10000 more particles
#                     above the threshold, then increased the threshold to 0.5
#                     and searched 10000 more above this new threshold.
#                     The same gam (not simpler) was fitted to the new dataset.


# Done


