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

library(adaptMCMC)
library(posterior)  # manejar y resumir posteriores
library(tidybayes)  # no sé si lo uso
library(bayesplot)  # graficar posteriores

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

gam_formula <- formula(
  y ~

    # marginal effects
    s(intercept, bs = basis, k = k_side, fx = fixed) +
    s(vfi, bs = basis, k = k_side, fx = fixed) +
    s(tfi, bs = basis, k = k_side, fx = fixed) +
    s(slope, bs = basis, k = k_side, fx = fixed) +
    s(wind, bs = basis, k = k_side, fx = fixed) +
    s(steps, bs = basis, k = k_side, fx = fixed) +

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
                      x = "par_values") {

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
    plot(data[, response] ~ data[, x][, i], ylab = response,
         xlab = par_names[i], ylim = yy, main = mm,
         pch = 19, col = rgb(0, 0, 0, alpha))

    if(!is.null(best)) {
      points(data[1:best, response] ~ data[, x][1:best, i],
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
# prob_fn: function to transform overlap into probability. If not NULL,
#   the probability is used to resample particles.
# support: matrix with lower limits (first row) and upper limits (second row)
#   of all variables.
# spread_data: data to simulate fires.
# centre: "local" or "global". New particles are simulated from a MVN centred
#   at each resampled particle (focal), or at the global mean (global). If global,
#   sobol may be set to TRUE, if focal, sobol is set to FALSE.
explore_likelihood <- function(data, n = 600, var_factor = 1, n_best = "all",
                               p_best = 0.1, prob_fn = NULL,
                               support, spread_data,
                               centre = "local", sobol = TRUE,
                               ...) {

  ### Testo
  # data = sim1; n = 500; p_best = 0.15; prob_fn = prob_fun;
  # support = sup; spread_data = spread_data;
  # centre = "global"; sobol = TRUE;
  # sigma = 0.15; pow = 2
  ### endo testo

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
  if(is.null(prob_fn)) {
    weight <- data$overlap
  } else {
    weight <- data$prob
  }

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

  # Compute prob
  if(is.null(prob_fn)) {
    sim_result$prob <- sim_result$overlap
  } else {
    sim_result$prob <- prob_fn(sim_result$overlap, ...)
  }

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
    prob_fn = NULL,
    support, spread_data,
    centre = "global", sobol = TRUE,
    ...
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
                              prob_fn = prob_fn,
                              support = support, spread_data = spread_data,
                              centre = centre, sobol = sobol,
                              ...)

    m <- paste("Search wave ", last_wave + i, "; max overlap = ",
               round(max(res$overlap), 4), sep = "")
    message(m)
  }

  return(res)
}

# extract samples from a MCMC.parallel run and return a draws_array
# (posterior class)
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
  return(draws_arr)
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

# Data for spread ---------------------------------------------------------

fire_file <- "2008_3.rds" #"2012_53.rds"#2008_3.rds" # "1999_25j.rds"
full_data <- readRDS(file.path(data_dir, fire_file))
# subset data needed for spread (to be cloned across workers)
spread_data <- full_data[c("landscape", "burnable", "ig_rowcol",
                           "burned_layer", "burned_ids")]

bb <- get_bounds(fire_data = spread_data)
ext_alpha <- 30
ext_beta <- 30
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


# Exploring likelihood ---------------------------------------------------

# first wave
ss <- sobol(n = 1000, dim = n_coef)
particles <- scale_params(ss, sup)
sim1 <- similarity_simulate_parallel(particles, spread_data)
sim1$wave <- 0
colnames(sim1$par_values) <- par_names

# con p 0.15
like_sim <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, p_best = 0.15,
  var_factor = c(1.5, 2, 1),
  support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
wave_plot(like_sim, alpha = 0.01)


# Fitting surrogate (GAM) -------------------------------------------------

data_gam <- as.data.frame(cbind(like_sim$par_values,
                                y = like_sim$overlap_logit))
data_gam_beta <- as.data.frame(cbind(like_sim$par_values,
                                     y = like_sim$overlap))


k_side <- 31 # 30
k_int <- 11 # 10 * 10
gg <- 1
basis <- "cr"
fixed <- T

cl <- makeForkCluster(4)

# # normal model at logit scale
# bench1 <- microbenchmark(
#   "normal1" = like_model1 <- bam(
#     gam_formula, data = data_gam, method = "fREML", discrete = TRUE,
#     nthreads = parallel::detectCores() / 2, gamma = gg
#   ),
#   "normal2" = like_model2 <- bam(
#     gam_formula, data = data_gam, method = "fREML", discrete = FALSE,
#     cluster = cl, gamma = gg
#   ),
#   times = 1
# )
# bench1 # 29 vs 249 s
#
# par(mfrow = c(1, 2))
# plot(plogis(fitted(like_model1)) ~ like_sim$overlap, main = "normal1",
#      pch = 19, col = rgb(0, 0, 0, 0.05)); abline(0, 1, col = "red")
# plot(plogis(fitted(like_model2)) ~ like_sim$overlap, main = "normal2",
#      pch = 19, col = rgb(0, 0, 0, 0.05)); abline(0, 1, col = "red")
# par(mfrow = c(1, 1))
# # parece andar un poquito mejor sin discretization, pero igual es malo

# bam(y ~ s(steps, k = 11, bs = "cr", fx = TRUE), data = data_gam_beta, family = betar())


# beta model
bench2 <- microbenchmark(
  "beta1" = like_model_beta1 <- bam(
    gam_formula, data = data_gam_beta,#, family = betar()#,
    #method = "fREML"#, discrete = TRUE,
    nthreads = parallel::detectCores() / 2#, gamma = gg
  ),
  times = 1, control = list(warmup = 0)
)
bench2 # smaller model (k = 11 and 6, with no smoothing), 61 s
bench2 # smaller model (k = 31 and 11, with no smoothing),  s # TOO LONG, ABORTED
bench2 # normal model at overlap scale (k = 31 and 11, with no smoothing),  151 s
summary(like_model_beta1)


bench3 <- microbenchmark(
  "beta3" = {

    k_side <- 21 # 30
    k_int <- 7 # 10 * 10
    basis <- "cr"
    fixed <- T

    like_model_beta2 <- bam(
      gam_formula, data = data_gam_beta, family = betar(),
      method = "fREML", discrete = FALSE,
      cluster = cl
    )
  },
  times = 1, control = list(warmup = 0)
)
bench3 # 97 s


bench4 <- microbenchmark(
  "beta4" = {

    k_side <- 16 # 30
    k_int <- 6 # 10 * 10
    basis <- "cr"
    fixed <- T

    like_model_beta4 <- bam(
      gam_formula, data = data_gam_beta, family = betar(),
      method = "fREML", discrete = FALSE,
      cluster = cl
    )
  },
  times = 1, control = list(warmup = 0)
)
bench4 # 122
summary(like_model_beta3)

bench5 <- microbenchmark(
  "beta5" = {

    k_side <- 31 # 30
    k_int <- 11 # 10 * 10
    basis <- "cr"
    fixed <- T

    like_model_beta4 <- bam(
      gam_formula, data = data_gam_beta, family = betar(),
      method = "fREML", discrete = FALSE,
      cluster = cl
    )
  },
  times = 1, control = list(warmup = 0)
)
bench4 # 122
summary(like_model_beta3)

# parece que gam() tarda horrores y usa mucha RAM. Lo cancelé antes de que
# terminara.
# Será que tarda tanto porque le pedí gamma = 1e-10? Comparar luego.

# modelo en escala logit para las predictoras?
data_gam_beta_logit <- as.data.frame(cbind(like_sim$par_values_logit,
                                           y = like_sim$overlap))
bench6 <- microbenchmark(
  "beta6" = {

    k_side <- 21 # 30
    k_int <- 7 # 10 * 10
    basis <- "cr"
    fixed <- T

    like_model_beta6 <- bam(
      gam_formula, data = data_gam_beta_logit, family = betar(),
      method = "fREML", discrete = FALSE,
      cluster = cl
    )
  },
  times = 1, control = list(warmup = 0)
)
bench6 #  109 s


# check models

plot(fitted(like_model_beta1) ~ like_sim$overlap, main = "beta1", ylim = c(0, 0.65),
     pch = 19, col = rgb(0, 0, 0, 0.05)); abline(0, 1, col = "red")

plot(fitted(like_model_beta2) ~ like_sim$overlap, main = "beta2",
     pch = 19, col = rgb(0, 0, 0, 0.01)); abline(0, 1, col = "red")

plot(fitted(like_model_beta3) ~ like_sim$overlap, main = "beta3",
     pch = 19, col = rgb(0, 0, 0, 0.01)); abline(0, 1, col = "red")

plot(fitted(like_model_beta3) ~ fitted(like_model_beta2), main = "beta3",
     pch = 19, col = rgb(0, 0, 0, 0.01)); abline(0, 1, col = "red")

plot(fitted(like_model_beta4) ~ like_sim$overlap, main = "beta4",
     pch = 19, col = rgb(0, 0, 0, 0.01)); abline(0, 1, col = "red")

plot(fitted(like_model_beta6) ~ like_sim$overlap, main = "beta6",
     pch = 19, col = rgb(0, 0, 0, 0.01)); abline(0, 1, col = "red")

like_model <- like_model_beta3


# todos feos. ajustar sobre un kernel de overlap.



# Comparing simulation-based y surrogate-based MCMC -----------------------

# functions to evaluate log-likelihood, all made to sample at the unconstrained
# space

# simulate fire
fn_like_sim <- function(x, support = NULL, fire_data = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  ov <- similarity_simulate_simpler(xx, n_sim = 20, fire_data = fire_data)
  lp <- log(ov) + sum(dlogis(x, log = TRUE))
  return(lp)
}

# m1: fitted mean
fn_like_m1 <- function(x,  support = NULL, model = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  nd <- as.data.frame(matrix(xx, nrow = 1))
  colnames(nd) <- par_names
  mu <- predict(model, nd, type = "response")
  lp <- log(mu) + sum(dlogis(x, log = TRUE))
  return(lp)
}

# m2: fitted mean + noise
fn_like_m2 <- function(x,  support = NULL, model = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  nd <- as.data.frame(matrix(xx, nrow = 1))
  colnames(nd) <- par_names
  mu <- predict(model, nd, type = "response")
  phi <- exp(model$family$getTheta())
  a <- mu * phi
  b <- phi - a
  y <- rbeta(1, a, b)
  lp <- log(y) + sum(dlogis(x, log = TRUE))
  return(lp)
}

# m3: fitted mean + uncertainty over mean (no noise)
fn_like_m3 <- function(x,  support = NULL, model = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  nd <- as.data.frame(matrix(xx, nrow = 1))
  colnames(nd) <- par_names
  linpred <- predict(model, nd, se.fit = TRUE)
  mu <- plogis(rnorm(1, linpred$fit, linpred$se.fit))
  lp <- log(mu) + sum(dlogis(x, log = TRUE))
  return(lp)
}

# m4: fitted mean + uncertainty over mean + noise
fn_like_m4 <- function(x,  support = NULL, model = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  nd <- as.data.frame(matrix(xx, nrow = 1))
  colnames(nd) <- par_names
  linpred <- predict(model, nd, se.fit = TRUE)
  mu <- plogis(rnorm(1, linpred$fit, linpred$se.fit))
  phi <- exp(model$family$getTheta())
  a <- mu * phi
  b <- phi - a
  y <- rbeta(1, a, b)
  lp <- log(y) + sum(dlogis(x, log = TRUE))
  return(lp)
}


## MCMC

# initial values
best100 <- order(like_sim$overlap, decreasing = TRUE)[1:100]
make_init <- function() {
  logit_scaled_vec(like_sim$par_values[sample(best100, 1), ],
                   sup)
}
q70 <- quantile(like_sim$overlap, 0.7)
sigma_init <- cov(logit_scaled(like_sim$par_values[like_sim$overlap > q70, ],
                               support = sup))
sampling_iters <- 10000
adapt_iters <- 2000
nc <- 8

# sample 5 posteriors
r0 <- MCMC.parallel(fn_like_sim,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, fire_data = spread_data,
                    packages = "FireSpread")
saveRDS(r0, "files/pseudolikelihood_estimation/posterior-2008_3-simulation.rds")
gc()

r1 <- MCMC.parallel(fn_like_m1,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, model = like_model,
                    packages = "mgcv")
gc()

r2 <- MCMC.parallel(fn_like_m2,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, model = like_model,
                    packages = "mgcv")
gc()

r3 <- MCMC.parallel(fn_like_m3,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, model = like_model,
                    packages = "mgcv")
gc()

r4 <- MCMC.parallel(fn_like_m4,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, model = like_model,
                    packages = "mgcv")
gc()


# postprocessing
post_names <- c("simulation", "mu", "mu+noise", "mu+uncertainty", "mu+uncertainty+noise")
rall <- list(r0, r1, r2, r3, r4)
names(rall) <- post_names
arr_all <- lapply(rall, tidy_samples, adapt = adapt_iters, support = sup)
saveRDS(arr_all, "files/pseudolikelihood_estimation/posteriors-five-methods-2008_3.rds")

summs <- lapply(arr_all, function(a) {
  tt <- thin_draws(a, 10)
  return(summarise_draws(tt))
})
summs
# It seems like using just mu is the most efficient.

# Compare with densities
dflong <- do.call("rbind", lapply(arr_all, function(a) {
  tt <- thin_draws(a, 10)
  return(as_draws_df(tt, .nchains = 8))
}))
names(dflong) <- par_names
dflong$sampling <- factor(
  rep(post_names, each = nrow(dflong) / length(post_names)),
  levels = post_names
)

dflong <- dflong[, c(par_names, "sampling")]
dflonger <- pivot_longer(dflong, all_of(1:n_coef), values_to = "par_value",
                         names_to = "par_names")

ggplot(dflonger, aes(x = par_value, fill = sampling, color = sampling)) +
  geom_density(alpha = 0) +
  facet_wrap(vars(par_names), scales = "free") +
  theme(panel.grid.minor = element_blank()) +
  viridis::scale_color_viridis(discrete = TRUE)


# todas las model-based posteriors son parecidas.
# todas son bastante wiggly. Quizás sí convenga suavizar un poco.

# estimar de nuevo un gam con pocas bases y penalizado.
# Pero puede que la wiggliness sea por el MCMC en sí, que es ineficiente.
# Porque las densidades por simulación también son muy wiggly.

# Quizás el único problem exista con wind.


# Pero antes de volver a estimar otro gam sobre el overlap,
# probar muestrear con un overlap ponderado,
# ajustando un GAM a ese nuevo overlap.

like_sim$res <- residuals.gam(like_model, type = "response")
wave_plot(like_sim, "res")
# más resids donde hay más datos...

like_sim$par_values_logit <- logit_scaled(like_sim$par_values, sup)
apply(like_sim$par_values_logit, 2, range)

wave_plot(like_sim, x = "par_values_logit")



# Volviendo a los GPs -----------------------------------------------------

# Como los GAM andan medios flojos, vuelvo a los GP.
# Se pueden tomar puntos separados, que sean K puntos más cercanos a los
# K centroides de un KNN.

# Cuántos puntos parecen factibles?
library(GauPro)

# with GauPro
loglik_model <- gpkm(
  X = like_sim$par_values[1:2000, ],
  Z = like_sim$overlap[1:2000],
  parallel = F, useC = TRUE, nug.max = 2,
  kernel = "matern52"
)
# Estima 527 / 60 = 8.8 min


# # max fitted loglik value
# fitted_ll <- loglik_model$pred(particles_data_join$par_values[fit_this, ],
#                                se.fit = F)
# loglik_max_fitted <- max(fitted_ll)

# Get KNN

X <- matrix(rnorm(100), ncol = 2)

library(FNN)
# get nn
k1 <- get.knn(X, k = 10, algorithm = c("kd_tree"))

# compute centroids
k1 <- kmeans(X, 10, iter.max = 10, nstart = 1)

plot(X)
points(k1$centers, col = "red")


kk <- kmeans(like_sim$par_values, 2000)

plot(like_sim$par_values[, 1:2], col = rgb(0, 0, 0, 0.01), pch = 19)
plot(kk$centers[, 1:2], col = rgb(1, 0, 0, 0.1), pch = 19)
# cool

# closest points to centers:
kk$cluster
kk$betweenss
dist(kk$centers[1, , drop = F], like_sim$par_values[kk$cluster == 1, ])

a <- dist(rbind(kk$centers[1, , drop = F], like_sim$par_values[kk$cluster == 1, ]))
a[, 1]
c(a)[1:nrow()]
closest_ids <- sapply(1:nrow(kk$centers), function(i) {
  # print(i)
  ids_in <- which(kk$cluster == i)
  neighs <- like_sim$par_values[ids_in, , drop = F]
  n <- nrow(neighs)
  dd <- dist(rbind(kk$centers[i, , drop = F], neighs)) %>% as.numeric
  id <- which.min(dd[1:n])
  return(ids_in[id])
})

# check:
par(mfrow = c(1, 2))
plot(like_sim$par_values[, 1:2], col = rgb(0, 0, 0, 0.01), pch = 19,
     main = "centroids")
points(kk$centers[, 1:2], col = rgb(1, 0, 0, 0.05), pch = 19)

plot(like_sim$par_values[, 1:2], col = rgb(0, 0, 0, 0.01), pch = 19,
     main = "nearest to centroids")
points(like_sim$par_values[closest_ids, 1:2], col = rgb(0, 0, 1, 0.5), pch = 19)
par(mfrow = c(1, 1))

# with GauPro
gp1 <- gpkm(
  X = like_sim$par_values[closest_ids, ],
  Z = like_sim$overlap[closest_ids],
  parallel = F, useC = TRUE, nug.max = 2,
  kernel = "Gaussian"
)
# Estima 527 / 60 = 8.8 min


# quizás antes escalar los params a [0, 1], para que todas las dist pesen lo
# mismo.
# hacer k_means sucesivos en distintos rangos de overlap, asignando más
# datos a donde haya más overlap.
# Elegir los mejores 100 puntos? y luego, kmeans para abajo?
# Ajustar en escala logit, o en escala original?
# O en escala transformada?
# O logit de la transformada?


gp1fit <- gp1$pred(like_sim$par_values[closest_ids, ], se.fit = F)
plot(gp1fit ~ like_sim$overlap[closest_ids])
abline(0, 1, col = "red") # hermoso fit.

sss <- 1:nrow(like_sim)
test_ids <- sss[!(sss %in% closest_ids)]

gp1fit <- gp1$pred(like_sim$par_values[test_ids, ], se.fit = F)
plot(gp1fit ~ like_sim$overlap[test_ids], col = rgb(0, 0, 0, 0.5), pch = 19)
abline(0, 1, col = "red") # hermoso fit.

# bastante bien! mucho mejor que el gam aparently. Y hay room for improvement.

gp2 <- gpkm(
  X = like_sim$par_values[closest_ids, ],
  Z = like_sim$overlap_logit[closest_ids],
  parallel = F, useC = TRUE, nug.max = 2,
  kernel = "Gaussian"
)

gp2fit <- gp2$pred(like_sim$par_values[closest_ids, ], se.fit = F) %>% plogis
plot(gp2fit ~ like_sim$overlap[closest_ids], main = "train",
     col = rgb(0, 0, 0, 0.5), pch = 19)
abline(0, 1, col = "red") # hermoso fit.

gp2fit_test <- gp2$pred(like_sim$par_values[test_ids, ], se.fit = F) %>% plogis
plot(gp2fit_test ~ like_sim$overlap[test_ids], main = "test",
     col = rgb(0, 0, 0, 0.5), pch = 19)
abline(0, 1, col = "red") # hermoso fit.

# segundo test
filll <- sobol(n = 20000, dim = n_coef)
filll <- scale_params(filll, sup)

gp2fit_test_fill <- gp2$pred(filll, se.fit = F) %>% plogis
hist(gp2fit_test_fill)
abline(v = mean(like_sim$overlap[closest_ids]))

# si el GP tiene mucha incertidumbre, se va hacia la media. Y vemos que no
# se acumulan datos en la media.

sum(gp2fit_test_fill > 0.6) # 51 nomás, en 20000
51 / 20000 # reeee pocos.
hist(like_sim$overlap)

# mejorar el kmeans y selección de puntos, y probar distintas formas de muestrear.
# probar ajustar el GP en escala transformada.


# Overlap transform using kernel ------------------------------------------

# Rescale data to [0, 1] by
overlap_max <- max(like_sim$overlap)
like_sim$overlap_unit <- like_sim$overlap / (overlap_max + 1e-6)

# transform using gaussian kernel
scale <- 0.2
trans <- exp(-((1 - like_sim$overlap) / scale) ^ 2)
plot(trans ~ like_sim$overlap)
# plot(trans ~ like_sim$overlap_unit)
like_sim$prob <- trans
wave_plot(like_sim, "prob")

# quizás la exploración debería hacerse en la escala transformada?
# Sí, y para eso habría que definir un kernel general para todos los fuegos.
# Gaussiano me parece bien, pero la scale debería basarse en el valor
# de los overlaps residuales.


# Exploring transformed likelihood ----------------------------------------

# curve(prob_fun(x))
# curve(prob_fun(x, pow = 1.5), col = "red", add = T)
# for(p in 3:10) curve(prob_fun(x, pow = p), col = p, add = T)

# first wave
ss <- sobol(n = 1000, dim = n_coef)
particles <- scale_params(ss, sup)
sim1 <- similarity_simulate_parallel(particles, spread_data)
sim1$prob <- prob_fun(sim1$overlap)
sim1$wave <- 0
colnames(sim1$par_values) <- par_names

# con p 0.15
like_sim <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, p_best = 0.15, prob_fn = prob_fun,
  var_factor = c(1.5, 2, 1),
  support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE,
  sigma = 0.15, pow = 2
)
wave_plot(like_sim, "prob", alpha = 0.01)
wave_plot(like_sim, "overlap", alpha = 0.01)

like_sim2 <- explore_likelihood_iterate(
  n_waves = 20, data = sim1, n = 500, p_best = 0.15, prob_fn = NULL,
  var_factor = c(1.5, 2, 1),
  support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE
)
like_sim2$prob <- prob_fun(like_sim2$overlap)
wave_plot(like_sim2, "prob", alpha = 0.01)
wave_plot(like_sim2, "overlap", alpha = 0.01)

min(like_sim$prob)
min(like_sim2$prob)

min(like_sim$prob) %>% qlogis
min(like_sim2$prob)

curve(plogis(-8 + 30 * x), from = -5, to = 5)


# Fitting surrogate (GP) at prob scale ------------------------------------

# 2 GPs: one fitted at overlap scale, an other fitted at prob scale.
#   in the first, the search is done at overlap scale,
#   while in the second, the search is done at prob scale.

# If the posterior is defined by the prob, not overlap, the prob-gp posterior
# should be more similar to the simulated one.

points_ov <- spread_points(like_sim2, 2000, B = 100, support = sup)
points_prob <- spread_points(like_sim, 2000, B = NULL, support = sup,
                             iter.max = 1000, nstart = 30)
# OJO que a veces el coso no converge.

gp_ov <- gpkm(
  X = like_sim2$par_values[points_ov, ],
  Z = like_sim2$overlap_logit[points_ov],
  parallel = F, useC = TRUE, nug.max = 100,
  kernel = "Gaussian"
)

# prob model fitted at logit scale
like_sim$prob_logit <- qlogis(like_sim$prob)
gp_prob_logit <- gpkm(
  X = like_sim$par_values[points_prob, ],
  Z = like_sim$prob_logit[points_prob],
  parallel = F, useC = TRUE, nug.max = 100,
  kernel = "Gaussian"
)

# prob model fitted at log scale
like_sim$prob_log <- log(like_sim$prob)
gp_prob_log <- gpkm(
  X = like_sim$par_values[points_prob, ],
  Z = like_sim$prob_log[points_prob],
  parallel = F, useC = TRUE, nug.max = 100,
  kernel = "Gaussian"
)

# prob model fitted at prob scale
gp_prob <- gpkm(
  X = like_sim$par_values[points_prob, ],
  Z = like_sim$prob[points_prob],
  parallel = F, useC = TRUE, nug.max = 1,
  kernel = "Gaussian"
)

# check both GPs at test points
sss <- 1:nrow(like_sim)
points_ov_out <- sss[!(sss %in% points_ov)]
points_prob_out <- sss[!(sss %in% points_prob)]

# overlap model
gp_ov_test <- gp_ov$pred(like_sim2$par_values[points_ov_out, ], se.fit = F) %>% plogis
plot(gp_ov_test ~ like_sim2$overlap[points_ov_out], main = "overlap model",
     col = rgb(0, 0, 0, 0.5), pch = 19)
abline(0, 1, col = "red")

gp_ov_test_prob <- prob_fun(gp_ov_test)
plot(gp_ov_test_prob ~ like_sim2$prob[points_ov_out], main = "overlap model, prob",
     col = rgb(0, 0, 0, 0.5), pch = 19,
     ylim = c(0, 0.00025))
abline(0, 1, col = "red")
# ajusta muuuy mal en la escala prob. Claro, porque el overlap queda en un rango
# de la escala en donde pequeños errores tienen una consecuencia monstruosa.



# prob model, check prob scale but fitted at logit
gp_prob_logit_test <- gp_prob_logit$pred(like_sim$par_values[points_prob_out, ], se.fit = F) %>% plogis
plot(gp_prob_logit_test ~ like_sim$prob[points_prob_out], main = "prob model",
     col = rgb(0, 0, 0, 0.5), pch = 19)
abline(0, 1, col = "red")
# acá pasa lo mismo...

# prob model, logit scale
gp_prob_logit_test <- gp_prob_logit$pred(like_sim$par_values[points_prob_out, ], se.fit = F)
plot(gp_prob_logit_test ~ like_sim$prob_logit[points_prob_out], main = "prob model",
     col = rgb(0, 0, 0, 0.5), pch = 19)
abline(0, 1, col = "red")
# no está tan mal, pero claro, los errores se traducen muy fuertemente en
# la escala de prob.

# prob model, check prob scale but fitted at prob
gp_prob_test <- gp_prob$pred(like_sim$par_values[points_prob_out, ], se.fit = F)
plot(gp_prob_test ~ like_sim$prob[points_prob_out], main = "prob model",
     col = rgb(0, 0, 0, 0.5), pch = 19)
abline(0, 1, col = "red")

# and testing over training data? // perhaps its because spreading points
# leaves too few in important regions. Perhaps the sampling is not so proportional
# to probability.
gp_prob_test <- gp_prob$pred(like_sim$par_values[points_prob, ], se.fit = F)
plot(gp_prob_test ~ like_sim$prob[points_prob], main = "prob model",
     col = rgb(0, 0, 0, 0.5), pch = 19)
abline(0, 1, col = "red")



# prob model, check prob scale but fitted at log
gp_prob_test <- gp_prob_log$pred(like_sim$par_values[points_prob_out, ], se.fit = F) %>% exp
plot(gp_prob_test ~ like_sim$prob[points_prob_out], main = "prob model",
     col = rgb(0, 0, 0, 0.5), pch = 19)
abline(0, 1, col = "red")

wave_plot(like_sim, "prob")
wave_plot(like_sim[points_prob, ], "prob")



# Y un gam en escala de prob puntuda? -------------------------------------

# que es lo que dije que iba a hacer.
data_gam_beta <- as.data.frame(cbind(like_sim$par_values,
                                     y = like_sim$prob))

like_sim$prob_01 <- like_sim$prob / (max(like_sim$prob) + 1e-6)
data_gam <- as.data.frame(cbind(like_sim$par_values,
                                y = qlogis(like_sim$prob_01)))
hist(data_gam$y)
k_side <- 31 # 30
k_int <- 11 # 10 * 10
gg <- 1
basis <- "cr"
fixed <- F

# cl <- makeForkCluster(4)
gam_prob_01 <- bam(
  gam_formula, data = data_gam,
  method = "fREML", discrete = T, nthreads = 8
)
# beta es imposible, tarda demasiado.
# Será que Normal también tardaba demasiado?

plot(plogis(fitted(gam_prob_01)) ~ like_sim$prob_01,
     col = rgb(0, 0, 0, 0.5), pch = 19)
abline(0, 1, col = "red")

like_sim$prob_gam <- fitted(gam_prob)
wave_plot(like_sim, "prob")
wave_plot(like_sim, "prob_gam")


# Un KNN? -----------------------------------------------------------------

library(caret)

knn1 <- knnreg(like_sim$par_values, like_sim$prob_01, 5)
pred_y <- predict(knn1, data.frame(like_sim$par_values))
plot(pred_y ~ like_sim$prob_01); abline(0, 1, col = "red")


# fire-wise probability function ------------------------------------------

# the function should take 1 at max(overlap)_f, with a fixed scale around 0.25.
# Hence, the overlap function will not be so peaky.


like_sim2$prob_local <- prob_fun(like_sim2$overlap, 0.25, rescale = TRUE)
plot(like_sim2$prob_local ~ like_sim2$overlap)
wave_plot(like_sim2, "prob_local") # menos puntudo.

prob_local_fitted <- gp_ov$pred(like_sim2$par_values, se.fit = F) %>% plogis

plot(prob_local_fitted[points_ov_out] ~ like_sim2$prob_local[points_ov_out],
     main = "overlap model, prob local, test",
     col = rgb(0, 0, 0, 0.5), pch = 19)
abline(0, 1, col = "red")

plot(prob_local_fitted[points_ov] ~ like_sim2$prob_local[points_ov],
     main = "overlap model, prob local, train",
     col = rgb(0, 0, 0, 0.5), pch = 19)
abline(0, 1, col = "red")

# fit gp at the prob_local scale (plogis?)
wave_plot(like_sim2[points_ov, ], "prob_local") # menos puntudo.
gp_ov <- gpkm(
  X = like_sim2$par_values[points_ov, ],
  Z = like_sim2$overlap_logit[points_ov],
  parallel = F, useC = TRUE, nug.max = 100,
  kernel = "Gaussian"
)


# Comparing posteriors (solo en overlap scale) --------------------------

# functions to compute likelihood

# GP fitted at overlap scale
fn_like_gp_ov <- function(x, support = NULL, model = NULL) {
  xx <- matrix(invlogit_scaled_vec(x, support), nrow = 1)
  like <- plogis(model$pred(xx, se.fit = F))
  lp <- log(like) + sum(dlogis(x, log = TRUE))
  return(lp)
}


## MCMC

# initial values
best100 <- order(like_sim2$overlap, decreasing = TRUE)[1:100]
make_init <- function() {
  logit_scaled_vec(like_sim2$par_values[sample(best100, 1), ],
                   sup)
}
q70 <- quantile(like_sim2$prob, 0.7)
sigma_init <- cov(logit_scaled(like_sim2$par_values[like_sim2$prob > q70, ],
                               support = sup))
sampling_iters <- 10000
adapt_iters <- 2000
nc <- 8

# sample

# based on simulation was previously sampled
r0 <- readRDS("files/pseudolikelihood_estimation/posterior-2008_3-simulation.rds")

r1 <- MCMC.parallel(fn_like_gp_ov,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, model = gp_ov,
                    packages = "GauPro")
gc()

# fn_like_gp_ov(make_init(), support = sup, model = gp_ov)


# postprocessing
post_names <- c("simulation", "gp_ov")
rall <- list(r0, r1)
names(rall) <- post_names
arr_all <- lapply(rall, tidy_samples, adapt = adapt_iters, support = sup)
saveRDS(arr_all, "files/pseudolikelihood_estimation/posteriors-three-methods-2008_3-gp.rds")

summs <- lapply(arr_all, function(a) {
  tt <- thin_draws(a, 10)
  return(summarise_draws(tt))
})
summs
# It seems like using just mu is the most efficient.

# Compare with densities
dflong <- do.call("rbind", lapply(arr_all, function(a) {
  tt <- thin_draws(a, 10)
  return(as_draws_df(tt, .nchains = 8))
}))
names(dflong) <- par_names
dflong$sampling <- factor(
  rep(post_names, each = nrow(dflong) / length(post_names)),
  levels = post_names
)

dflong <- dflong[, c(par_names, "sampling")]
dflonger <- pivot_longer(dflong, all_of(1:n_coef), values_to = "par_value",
                         names_to = "par_names")

ggplot(dflonger, aes(x = par_value, fill = sampling, color = sampling)) +
  geom_density(alpha = 0) +
  facet_wrap(vars(par_names), scales = "free") +
  theme(panel.grid.minor = element_blank()) +
  viridis::scale_color_viridis(discrete = TRUE)


# 800 4 h
# 80000
# 400 h
# 400 / 24 = 17 días (si usamos este approach.)

# Comparing posteriors ----------------------------------------------------

# functions to compute likelihood

# simulate fire
fn_like_sim <- function(x, support = NULL, fire_data = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  ov <- similarity_simulate_simpler(xx, n_sim = 20, fire_data = fire_data)
  like <- prob_fun(ov)
  lp <- log(like) + sum(dlogis(x, log = TRUE))
  return(lp)
}

# GP fitted at overlap scale
fn_like_gp_ov <- function(x, support = NULL, model = NULL) {
  xx <- matrix(invlogit_scaled_vec(x, support), nrow = 1)
  ov <- plogis(model$pred(xx, se.fit = F))
  like <- prob_fun(ov)
  lp <- log(like) + sum(dlogis(x, log = TRUE))
  return(lp)
}

# GP fitted at prob scale
fn_like_gp_prob <- function(x, support = NULL, model = NULL) {
  xx <- matrix(invlogit_scaled_vec(x, support), nrow = 1)
  like <- plogis(model$pred(xx, se.fit = F))
  lp <- log(like) + sum(dlogis(x, log = TRUE))
  return(lp)
}

## MCMC

# initial values
best100 <- order(like_sim$overlap, decreasing = TRUE)[1:100]
make_init <- function() {
  logit_scaled_vec(like_sim$par_values[sample(best100, 1), ],
                   sup)
}
q70 <- quantile(like_sim$prob, 0.7)
sigma_init <- cov(logit_scaled(like_sim$par_values[like_sim$prob > q70, ],
                               support = sup))
sampling_iters <- 10000
adapt_iters <- 2000
nc <- 8

# sample 3 posteriors
r0 <- MCMC.parallel(fn_like_sim,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, fire_data = spread_data,
                    packages = "FireSpread")
saveRDS(r0, "files/pseudolikelihood_estimation/posterior-2008_3-simulation-prob.rds")
gc()

r1 <- MCMC.parallel(fn_like_gp_ov,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, model = gp_ov,
                    packages = "GauPro")
gc()

r2 <- MCMC.parallel(fn_like_gp_prob,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, model = gp_prob,
                    packages = "GauPro")
gc()

# postprocessing
post_names <- c("simulation", "gp_ov", "gp_prob")
rall <- list(r0, r1, r2)
names(rall) <- post_names
arr_all <- lapply(rall, tidy_samples, adapt = adapt_iters, support = sup)
saveRDS(arr_all, "files/pseudolikelihood_estimation/posteriors-three-methods-2008_3-gp-prob.rds")

summs <- lapply(arr_all, function(a) {
  tt <- thin_draws(a, 10)
  return(summarise_draws(tt))
})
summs
# It seems like using just mu is the most efficient.

# Compare with densities
dflong <- do.call("rbind", lapply(arr_all, function(a) {
  tt <- thin_draws(a, 10)
  return(as_draws_df(tt, .nchains = 8))
}))
names(dflong) <- par_names
dflong$sampling <- factor(
  rep(post_names, each = nrow(dflong) / length(post_names)),
  levels = post_names
)

dflong <- dflong[, c(par_names, "sampling")]
dflonger <- pivot_longer(dflong, all_of(1:n_coef), values_to = "par_value",
                         names_to = "par_names")

ggplot(dflonger, aes(x = par_value, fill = sampling, color = sampling)) +
  geom_density(alpha = 0) +
  facet_wrap(vars(par_names), scales = "free") +
  theme(panel.grid.minor = element_blank()) +
  viridis::scale_color_viridis(discrete = TRUE)


# Esto fue un gran fracaso!!! los GP en la escala prob funcionan muuuuuy mal.


# Aceptación por threshold (GAM bernoulli) ---------------------------------



# SEGUIR ACÁAAAAAA --------------------------------------------------------



# VOLVER A INTENTAR EL MODELO COMPLETAMENTE MIXTO.

# Si ajustamos un GP a los puntos de overlap más alto, quizás tenemos mejor
# precisión. Y cuando el GP tiene precisión baja, ahí rechazamos la partícula,
# asumiendo que está fuera de rango.

# Por otro lado, si un GP tiene media, por ej, 0.4, cuando deriva tiende
# a la media. Y si nuestro umbral es un valor mayor a la media, rara vez
# vamos a estar aceptando partículas fuera de rango.


hist(like_sim$overlap)
wave_plot(like_sim, alpha = .01)
sum(like_sim$overlap > 0.5)

# y un gam bernoulli??
data_gam_bern <- as.data.frame(cbind(like_sim$par_values,
                                     y = as.numeric(like_sim$overlap > 0.5)))
k_side <- 21 # 30
k_int <- 7 # 10 * 10
gg <- 1
basis <- "cr"
fixed <- F

# cl <- makeForkCluster(4)
gam_bern <- bam(
  gam_formula, data = data_gam_bern, family = "binomial",
  method = "fREML", discrete = T, nthreads = 8
)
summary(gam_bern)

boxplot(fitted(gam_bern) ~ data_gam_bern$y)
# Parece bastante bueno.
hist(fitted(gam_bern))
mean(data_gam_bern$y)
mean(fitted(gam_bern))

plot(fitted(gam_bern) ~ like_sim$overlap,
     col = rgb(0, 0, 0, 0.05), pch = 19)
abline(v = 0.5, col = "red")
# esto da buena idea del ajuste

like_sim$prob_gt05 <- fitted(gam_bern)
wave_plot(like_sim, alpha = 0.1)
wave_plot(like_sim, "prob_gt05", alpha = 0.1)

ppp <- predict(gam_bern, se.fit = T)
lll <- plogis(ppp$fit - qnorm(0.975) * ppp$se.fit)
uuu <- plogis(ppp$fit + qnorm(0.975) * ppp$se.fit)
like_sim$prob_gt05_ic_width <- uuu - lll

wave_plot(like_sim, "prob_gt05_ic_width", alpha = 0.1)

plot(prob_gt05_ic_width ~ overlap, data = like_sim, pch = 19,
     col = rgb(0, 0, 0, 0.05), ylab = "uncertainty")

plot(prob_gt05_ic_width ~ par_values[, "wind"], data = like_sim, pch = 19,
     col = rgb(0, 0, 0, 0.05), ylab = "uncertainty", xlab = "wind")

plot(prob_gt05_ic_width ~ prob_gt05, data = like_sim, pch = 19,
     col = rgb(0, 0, 0, 0.1), ylab = "uncertainty", xlab = "fitted")

# Juaaaaa claro, el mayor problema deben ser los overlap ~ 0 que tienen
# re pocos datos.

# mirar la distancia los k-vecinos más próximos de cada punto,
# y fijarse si eso explica un poco la uncertainty.

# reproducir partículas inciertas con un kernel pequeño



# En realidad, si tengo mucha incertidumbre pero la fitted es muy baja, no
# debería sobremuestrear una parte.
# El problema son las fitted altas con incertidumbre alta.
# O sea, juntar más partículas donde haya más incertidumbre pero también overlap
# más alto. Tendrá sentido??


# Comparando posteriores usando rechazo en 0.5 -------------------------------

# functions to compute likelihood

# simulate fire
fn_like_sim_bin <- function(x, support = NULL, fire_data = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  ov <- similarity_simulate_simpler(xx, n_sim = 20, fire_data = fire_data)
  if(ov > 0.5) {
    return(sum(dlogis(x, log = TRUE)))
  } else {
    return(-Inf)
  }
}

fn_like_gam_bin <- function(x,  support = NULL, model = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  nd <- as.data.frame(matrix(xx, nrow = 1))
  colnames(nd) <- par_names
  mu <- predict(model, nd, type = "response")
  y <- rbinom(1, 1, mu)
  if(y == 1) {
    return(sum(dlogis(x, log = TRUE)))
  } else {
    return(-Inf)
  }
}

fn_like_gam_bin_unc <- function(x,  support = NULL, model = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  nd <- as.data.frame(matrix(xx, nrow = 1))
  colnames(nd) <- par_names
  pred <- predict(model, nd, se.fit = T)
  mu <- plogis(rnorm(1, pred$fit, pred$se.fit))
  y <- rbinom(1, 1, mu)
  if(y == 1) {
    return(sum(dlogis(x, log = TRUE)))
  } else {
    return(-Inf)
  }
}

fn_like_gam_prob <- function(x,  support = NULL, model = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  nd <- as.data.frame(matrix(xx, nrow = 1))
  colnames(nd) <- par_names
  pred <- predict(model, nd, se.fit = T)
  mu <- predict(model, nd, type = "response")
  lp <- log(mu) + sum(dlogis(x, log = TRUE))
  return(lp)
}
## MCMC

# initial values
best100 <- order(like_sim2$overlap, decreasing = TRUE)[1:100]

make_init <- function() {
  logit_scaled_vec(like_sim2$par_values[sample(best100, 1), ],
                   sup)
}
sigma_init <- cov(logit_scaled(like_sim2$par_values[like_sim2$overlap > 0.5, ],
                               support = sup))
sampling_iters <- 10000
adapt_iters <- 2000
nc <- 8

# sample
r0 <- MCMC.parallel(fn_like_sim_bin,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, fire_data = spread_data,
                    packages = "FireSpread")
# saveRDS(r0, "files/pseudolikelihood_estimation/posterior-simulation-2008_3-binary.rds")

sampling_iters <- 20000
adapt_iters <- 2000
nc <- 15

r1 <- MCMC.parallel(fn_like_gam_bin,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, model = gam_bern,
                    packages = "mgcv")
gc()

r2 <- MCMC.parallel(fn_like_gam_bin_unc,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, model = gam_bern,
                    packages = "mgcv")
gc()

r3 <- MCMC.parallel(fn_like_gam_prob,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = nc, n.chain = nc,
                    scale = sigma_init, init = make_init(), acc.rate = 0.234,
                    support = sup, model = gam_bern,
                    packages = "mgcv")

# postprocessing
post_names <- c("simulation_bin", "gam_bin", "gam_bin_unc", "gam_bin_prob")
rall <- list(r0, r1, r2, r3)
names(rall) <- post_names
arr_all <- lapply(rall, tidy_samples, adapt = adapt_iters, support = sup)
# saveRDS(arr_all, "files/pseudolikelihood_estimation/posteriors-three-methods-2008_3-binary.rds")

arr_all <- readRDS("files/pseudolikelihood_estimation/posteriors-three-methods-2008_3-binary.rds")

summs <- lapply(arr_all, function(a) {
  tt <- thin_draws(a, 10)
  return(summarise_draws(tt))
})
summs
# It seems like using just mu is the most efficient.

# Compare with densities
dflong <- do.call("rbind", lapply(arr_all, function(a) {
  tt <- thin_draws(a, 10)
  return(as_draws_df(tt, .nchains = 8))
}))
names(dflong) <- par_names
dflong$sampling <- factor(
  rep(post_names, each = nrow(dflong) / length(post_names)),
  levels = post_names
)

dflong <- dflong[, c(par_names, "sampling")]
dflonger <- pivot_longer(dflong, all_of(1:n_coef), values_to = "par_value",
                         names_to = "par_names")
dflonger$par_names <- factor(dflonger$par_names, levels = par_names)

ggplot(dflonger, aes(x = par_value, fill = sampling, color = sampling)) +
  geom_density(alpha = 0) +
  facet_wrap(vars(par_names), scales = "free") +
  theme(panel.grid.minor = element_blank()) +
  viridis::scale_color_viridis(discrete = TRUE, end = 0.9, option = "A")

# Si esto anda bien, la "likelihood" será la probabilidad de superar un umbral.
# Pero mejor no considerar incertidumbre sobre la probabilidad estimada.
# Se ve que con _unc hay problemas porque la verdadera posterior no es
# multivariada normal. Y le cuesta mucho más muestrearlo.


# probar agrandando el modelo! más bases
# Aunque el modelo actual parece tener demasiado ruido.
# quizás con menos bases?


# NONO, ojo, parece que a muchos les cuesta converger, sobre todo a viento.
# hay que correr más cadenas y más iters con los gams. chequear eso antes de
# largarse a modificar el modelo.

# Chequeo posteriores problemáticas ---------------------------------------

library(bayesplot)

x <- rename_variables(arr_all$simulation_bin,
                      "intercept" = "...1",
                      "vfi" = "...2",
                      "tfi" = "...3",
                      "slope" = "...4",
                      "wind" = "...5",
                      "steps" = "...6")
variables(arr_all$simulation_bin)
variables(x)


arr_all <- lapply(arr_all, function(a) {
  rename_variables(a,
                   "intercept" = "...1",
                   "vfi" = "...2",
                   "tfi" = "...3",
                   "slope" = "...4",
                   "wind" = "...5",
                   "steps" = "...6")
})


for(i in 1:length(arr_all)) {
  p <-
    mcmc_dens_overlay(arr_all[[i]]) +
    scale_color_viridis(discrete = TRUE, option = "A", end = 0.8) +
    theme(panel.grid.minor = element_blank()) +
    ggtitle(names(arr_all)[i])
  print(p)
}


