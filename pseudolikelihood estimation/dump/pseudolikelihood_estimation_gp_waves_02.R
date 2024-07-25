# Here I run waves for 2008_3 at the overlap scale and correcting the bugs.
# Antes sd_int era 5, ahora lo puse en 20

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)

library(Rcpp)
library(terra)
library(abind)         # manage arrays

library(randtoolbox)   # sobol sequences
library(GauPro)        # fit gaussian process
library(mgcv)          # fit spline to choose non-bounded particles

library(foreach)       # parallelization
library(doMC)          # doRNG had problems with cpp functions

library(FireSpread)    # spread and similarity functions

source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for rast_from_mat and a few constants

# source("estimation_functions.R") # prior_dist and other stuff

# Multicore settings -----------------------------------------------------

n_cores <- 12
registerDoMC(n_cores)

# Data and constants -----------------------------------------------------

data_dir <- file.path("data", "focal fires data", "landscapes_ig-known")
filenames <- list.files(data_dir)

# dir to save output
target_dir <- file.path("files", "pseudolikelihood_estimation")

# constants for fire spread simulation
upper_limit <- 0.5
n_coef <- 5
par_names <- c("intercept", "vfi", "tfi", "slope", "wind")

# number of particles to choose by wave
n_pw <- 1000

# number of fires to simulate by particle
n_sim <- 10

# formula for in_bounds model (gam)
bounds_model_formula <- formula(
  in_bounds ~
    s(intercept, k = 8, bs = "cr") +
    s(vfi, k = 8, bs = "cr") +
    s(tfi, k = 8, bs = "cr") +
    s(slope, k = 8, bs = "cr") +
    s(wind, k = 8, bs = "cr")
)

################
# data for testing
i = 1
full_data <- readRDS(file.path(data_dir, "2008_3.rds"))
# full_data <- readRDS(file.path(data_dir, "1999_25j.rds"))
# subset data needed for spread (to be cloned across workers)
spread_data <- full_data[c("landscape", "burnable", "ig_rowcol",
                           "burned_layer", "burned_ids")]

# Functions ---------------------------------------------------------------

## add function to prepare data ------------------------------------------
# leave burned and burned_ids layers with the names needed by simulate_fire()


# prior distribution to simulate parameters or to compute them from a sobol
# sequence (type = "quantile", which computes the icdf.)
prior_dist <- function(mu_int = 0, sd_int = 20,
                       r_slope = 0.04,
                       r_fi = 0.15,
                       r_wind = 0.3,
                       type = "rng", # or "quantile"
                       n = 1,
                       p = NULL) {

  if(type == "rng") {
    b <- cbind(
      "intercept" = rnorm(n, mu_int, sd_int),
      "vfi" = rexp(n, r_fi),
      "tfi" = rexp(n, r_fi),
      "slope" = rexp(n, r_slope),
      "wind" = rexp(n, r_wind)
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
    if(length(p) %% 5 != 0) stop("p must be a multiple of the number of parameters (5).")

    if(length(dim(p)) < 2) p <- matrix(p, ncol = 5)
    if(length(dim(p)) == 2) {
      if(ncol(p) != 5) p <- matrix(as.numeric(p), ncol = 5)
    }

    q <- cbind(
      "intercept" = qnorm(p[, 1], mu_int, sd_int),
      "vfi" = qexp(p[, 2], r_fi),
      "tfi" = qexp(p[, 3], r_fi),
      "slope" = qexp(p[, 4], r_slope),
      "wind" = qexp(p[, 5], r_wind)
    )

    return(q)
  }
}

# function to make particles from a sobol sequence
particles_sim <- function(N = 100, d = n_coef) {
  prior_dist(type = "quantile",
             p = sobol(N, dim = d, seed = 123, init = TRUE))
}


# function to assess if simulated fires reached the landscape bounds
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
similarity_simulate_particle <- function(particle, n_sim = 10,
                                         fire_data = NULL) {

  ## testo
  # particle <- particles_sim(N = 1)
  ## end testo
  metrics <- matrix(NA, n_sim, 3)
  colnames(metrics) <- c("overlap", "size_diff", "edge")

  # simulate and compute metrics
  for(i in 1:n_sim) { # simulated fires by particle
    fire_sim <- simulate_fire_compare(
      landscape = fire_data$landscape,
      burnable = fire_data$burnable,
      ignition_cells = fire_data$ig_rowcol,
      coef = particle,
      upper_limit = upper_limit
    )

    metrics[i, "overlap"] <- overlap_spatial(
      fire_sim, fire_data[c("burned_layer", "burned_ids")]
    )

    metrics[i, "size_diff"] <- ncol(fire_sim$burned_ids) -
      ncol(fire_data$burned_ids)

    metrics[i, "edge"] <- edge_count(fire_sim$burned_layer)
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
    edge =  mean(metrics[, "edge"])
  )

  return(ll_summ)
}

# Emulate the loglik over a list of particles in parallel. Returns the metrics
# in rows.
similarity_simulate_parallel <- function(particle_ids = 1:100,
                                         fire_data = NULL) {

  # turn particle matrix into list for parallel evaluation
  particles_mat_local <- particles_all[particle_ids, , drop = F]
  particles_list <- lapply(1:nrow(particles_mat_local),
                           function(x) particles_mat_local[x, ])

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
  res$particle_id <- particle_ids
  res$wave <- NA
  res$par_values <- particles_mat_local

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


# in_bounds: evaluate whether the particles are in_bounds, returning a binary
# vector to fit in_bounds model.
# similarity_bounds argument as provided by the get_bounds function.
in_bounds <- function(data, similarity_bounds,
                      prop = 0.85, var_factor = 2,
                      edge_prop_upper = 0.05,
                      edge_abs_upper = NULL) {

  # # # test
  # data <- particles_data_new
  # similarity_bounds <- bb
  # var_factor = 10
  # prop = 0.85
  # var_factor = 100
  # edge_prop_upper = 0.05
  # # # end testo

  size_diff_max <- similarity_bounds["largest", "size_diff"]
  size_diff_min <- similarity_bounds["smallest", "size_diff"]

  size_diff_lower <- prop * size_diff_min
  size_diff_upper <- prop * size_diff_max

  # get variance in overlap
  low_var <- c(
    similarity_bounds["largest", "var"], # overlap variance, original scale
    similarity_bounds["smallest", "var"]
  ) %>% max

  lowest_var <- low_var * var_factor

  edge_upper <- ifelse(
    is.null(edge_abs_upper),
    edge_prop_upper * similarity_bounds["largest", "edge"],
    edge_abs_upper
  )

  too_large <- ((data$size_diff >= size_diff_upper) &
                  (data$var <= lowest_var)) |
    (data$edge > edge_upper)
  too_small <- (data$size_diff <= size_diff_lower) &
    (data$var <= lowest_var)

  keep <- as.numeric(!(too_large | too_small))

  return(keep)
}


# fit the in_bounds model, specifying the data and a formula
fit_in_bounds <- function(
    particles_data,
    form = bounds_model_formula
) {

  data <- cbind(in_bounds = particles_data$in_bounds,
                particles_data$par_values) %>% as.data.frame
  m <- gam(form, data = data, family = binomial(link = "logit"), method = "REML")

  return(m)
}

# function to simulate which particles are to be evaluated, based on the relative
# probability of being in bounds.
simulate_in_bounds <- function(in_model, particle_ids) {

  # get maximum fitted probability to relativize predictions.
  pmax <- predict(in_model, type = "response") %>% max

  # predict in_range_prob
  prob <- predict(in_model, as.data.frame(particles_all[particle_ids, ]),
                  type = "response") / pmax
  prob[prob > 1] <- 1

  use_bin <- rbinom(length(prob), size = 1, prob = prob)

  return(use_bin)
}


# function to fit a GP and optimize the surface.
gp_fit <- function(data, response = "overlap", gp_optimize = T) {

  message(paste("Fitting Gaussian Process, N =", nrow(data)))
  loglik_model <- gpkm(
    X = data$par_values,
    Z = data[, response],
    parallel = F, useC = TRUE, nug.max = 100,
    kernel = "matern52"
  )

  # define loglik_high associated with this new GP and the new loglik_threshold,
  # to be used in the next wave
  fitted_ll <- loglik_model$pred(data$par_values, se.fit = F)
  loglik_high <- max(fitted_ll)

  # find MLE in the fitted GP
  if(gp_optimize) {
    params_max_fitted <- data$par_values[which.max(fitted_ll), ]

    # bounds
    lowers <- apply(data$par_values, 2, min)
    uppers <- apply(data$par_values, 2, max)

    message("Optimizing pseudo-likelihood surface")
    op <- optim(params_max_fitted, gp_predict, loglik_model = loglik_model,
                se.fit = FALSE,
                control = list(fnscale = -1, maxit = 1e5),
                method = "L-BFGS-B", lower = lowers, upper = uppers)
  }

  # tidy result
  res <- list(
    loglik_model = loglik_model,
    response = response,
    op = if(gp_optimize) op else NULL,
    loglik_high = loglik_high,
    data = data # for plotting of independent fit
  )

  return(res)
}

# for GauPro gp
gp_predict <- function(x, loglik_model, se.fit = F) {
  nd <- matrix(x, nrow = 1)
  return(loglik_model$pred(nd, se.fit = se.fit))
}

# Function to filter particles based on the previously fitted gps
gp_filter <- function(gp_list = NULL, loglik_thresholds = NULL,
                      particle_ids = NULL) {

  # compute predictions for all particles at all GPs
  preds <- lapply(gp_list, function(x) {
    x$pred(particles_all[particle_ids, , drop = F],
           se.fit = T)
  })

  # extract mean and sd
  pred_means <- do.call("cbind", lapply(1:length(preds),
                                        function(x) preds[[x]]$mean))
  pred_sds <- do.call("cbind", lapply(1:length(preds),
                                      function(x) preds[[x]]$se))

  # evaluate plausibility in all previous GPs, using each ones' loglik_max
  keep_wide <- do.call("cbind", lapply(1:ncol(pred_means), function(x) {
    keep <- (pred_means[, x] + 3 * pred_sds[, x]) >= (loglik_thresholds[x])
    return(as.numeric(keep))
  }))

  keep_these <- which(rowSums(keep_wide) == length(preds))

  return(particle_ids[keep_these])
}

# This just searchs for more particles to simulate, returning a particle_id
# vector.
more_particles <- function(gp_list = NULL,
                           loglik_threshold = NULL,
                           in_model = NULL,
                           last_particle = NULL,
                           n_pw = 1000) {

  if(!is.null(gp_list)) {

    added <- 0
    new_ids <- NULL

    while(added < n_pw) {

      try_ids <- (last_particle + 1) :
        (last_particle + n_pw * (1 + length(gp_list)))

      # If you run out of particles, make more
      if(max(try_ids) > nrow(particles_all)) {
        particles_all <<- particles_sim(nrow(particles_all) + 1e5)
        # this is inefficient because recreates the sequence from scratch,
        # but it ensures the particles used are always the same even if
        # the computation is run over separated R sessions.
      }

      in_range_ids <- try_ids[simulate_in_bounds(in_model, try_ids) == 1]

      new_ids <- c(new_ids,
                   gp_filter(gp_list = gp_list,
                             loglik_threshold = loglik_threshold,
                             particle_ids = in_range_ids))

      added <- length(new_ids)
      last_particle <- try_ids[length(try_ids)]
    }

    new_particles <- new_ids[1:n_pw]
  } else {
    new_particles <- 1:n_pw
  }

  return(new_particles)
}


# Function to subset particles to fit a GP.
# It takes the best n_best plus n_others sampled at random.
# It returns the particle_ids to use for GP fitting
subset_particles <- function(data, n_fit = 800, prop_best = 0.3,
                             use_in_bounds = F) {

  if(use_in_bounds) data <- data[data$in_bounds == 1, ]

  if(nrow(data) < n_fit) {
    # warning("subset_particles: Small dataset, no need to subset.")
    return(data$particle_id)
  }

  n_best <- ceiling(n_fit * prop_best)
  n_others <- n_fit - n_best

  data <- data[order(data$ll, decreasing = TRUE), ]
  ids_best <- data$particle_id[1:n_best]
  ids_others <- sample(data$particle_id[(n_best+1) : nrow(data)],
                       size = n_others, replace = F)

  return(c(ids_best, ids_others))
}



# Function to get a loglik_threshold, which can be based on a maximum - tolerance
# or a low loglik value. If all is provided, the max (more astringent) is returned.
get_loglik_threshold <- function(loglik_high = NULL,
                                 loglik_tolerance = NULL,
                                 loglik_lower = NULL) {
  # define likelihood-based threshold for counting particles
  relative_lower <- if(is.null(loglik_tolerance) | is.null(loglik_high)) NULL else {
    loglik_high - loglik_tolerance
  }

  if(is.null(relative_lower) & is.null(loglik_lower)) {
    stop("No sufficient loglik information provided.")
  } else {
    loglik_threshold <- max(relative_lower, loglik_lower)
  }

  return(loglik_threshold)
}

# loglik_update_cumulative runs a wave of likelihood emulator fitting,
# as proposed by Wilkinson et al. 2014.

# returns
# gps: list of fitted Gaussian Processes, one by wave.
# in_model: gam model predicting the probability of particles laying in_bounds.
# loglik_optim: optim() output for the last GP (like a MLE).
# loglik_thresholds: the thresholds defined after every wave to choose good
#   particles.
# particles_data: data.frame with the following columns
#   ids: particle id,
#   wave: latest wave in which the particle was used to fit a GP,
#   ll: log of the mean of the simulated likelihoods (overlap) in each particle,
#   var: variance of the simulated likelihoods in each particle, in the original
#     scale (not log),
#   parameter_values: matrix with the parameter values.
#   in_bounds: integer [0, 1] indicating if the particles are in_bounds.

# arguments
# fire_data: data for simulate_fire_compare (list).
# similarity_bounds: data from get_bounds, to define bounded particles
# previous_wave: result from a previous wave (default is NULL, which
#   starts the estimation process),
# loglik_tolerance: new particles will be judged as good if they have high
#   probability of being > (loglik_high - loglik_tolerance).
#     See get_loglik_threshold.
# loglik_lower: alternatively, the particles are jugded as good if they are
#   > loglik_lower. If present, the loglik_threshold is the maximum between
#   (loglik_high - loglik_tolerance) and loglik_lower,
# n_pw: number of new good particles to include in the new wave,
# n_fit: maximum number of particles to use when fitting the GP. This is used
#   to avoid gps taking too long to fit when all simulated data is used. Instead,
#   the best particles are used.
# prop_best: Which proportion of n_fit should be the best particles? The remaining
#   are sampled at random.
# gp_optimize = should the fitted gp be optimized? Not recommended
#   if it's first fitted with many particles out of bounds. Useful for plotting
#   the MLE.

loglik_update_cumulative <- function(
    fire_data = NULL,
    similarity_bounds = NULL,
    previous_wave = NULL,
    loglik_tolerance = 3,
    loglik_lower = NULL,
    n_pw = 800,
    n_fit = 800,
    prop_best = 0.3,
    gp_optimize = TRUE,
    use_in_bounds = FALSE,
    n_sim = 20,
    response = "overlap"
) {

  #### testo
  # fire_data = spread_data
  # similarity_bounds = NULL
  # previous_wave = NULL
  # loglik_tolerance = 3
  # loglik_lower = NULL
  # n_pw = 800
  # n_fit = 800
  # prop_best = 0.3
  # gp_optimize = TRUE
  # use_in_bounds = FALSE
  # n_sim = 20
  # response = "overlap"
  ####

  # when a previous wave has been run
  if(!is.null(previous_wave)) {
    gp_list <- previous_wave$gps
    wave_last <- length(gp_list)    # number of previous waves
    gp_last <- gp_list[[wave_last]]
    message(paste("Wave", 1 + wave_last))

    ### agregar esto para que el threshold se calcule al comienzo
    # loglik_threshold_new <- get_loglik_threshold(loglik_high,
    #                                              loglik_tolerance,
    #                                              loglik_lower)
    ###



    loglik_thresholds <- previous_wave$loglik_thresholds

    in_model <- previous_wave$in_model
    particles_data <- previous_wave$particles_data

    # filter simulated particles according to last gp (this is the
    # post-simulation filter)
    message("Filtering old particles")
    ids_eval <- particles_data[particles_data$wave >= (wave_last - 1), "particle_id"]
    ids_keep <- gp_filter(list(gp_last), particle_ids = ids_eval,
                          loglik_thresholds = loglik_thresholds[wave_last])

    if(length(ids_keep) > 0) {
      ids_keep_index <- which(particles_data$particle_id %in% ids_keep)
      particles_data$wave[ids_keep_index] <- wave_last + 1
    }

    # get new particles meeting the same criterion
    message("Getting more particles")
    new_particles <- more_particles(gp_list = gp_list,
                                    loglik_threshold = loglik_thresholds,
                                    in_model = in_model,
                                    last_particle = max(particles_data$particle_id),
                                    n_pw = n_pw)

  } else {
    message("Wave 1")
    message("Getting more particles")
    new_particles <- more_particles(n_pw = n_pw)
  }

  # Define wave number
  this_wave <- ifelse(is.null(previous_wave),
                      1,
                      wave_last + 1)

  # Simulate loglik on new particles
  message("Simulating fires")

  particles_data_new <- similarity_simulate_parallel(
    particle_ids = new_particles,
    fire_data = fire_data
  )

  particles_data_new$wave <- this_wave

  # define in_bounds particles (binary)
  particles_data_new$in_bounds <- in_bounds(particles_data_new,
                                            similarity_bounds,
                                            prop = 0.85,
                                            var_factor = 2,
                                            edge_prop_upper = 0.05)

  # merge old and new particles data
  if (is.null(previous_wave)) {
    particles_data_join <- particles_data_new
  } else {
    particles_data_join <- rbind(particles_data, # the free object
                                 particles_data_new)
  }

  # fit or update in_bounds_model using all the data
  message("Fitting in_model")

  in_model <- fit_in_bounds(particles_data_join)

  # Define particles to fit the GP (too much takes too long)
  # If it's the first wave, use all particles to fit the loglik model.
  if(this_wave == 1) {
    rows_fit <- 1:nrow(particles_data_join)
  } else {
    particles_use <- subset_particles(
      particles_data_join[particles_data_join$wave == this_wave, ],
      n_fit, prop_best,
      use_in_bounds = use_in_bounds
    )

    # ensure you fit the GP with more than 100 points
    if(length(particles_use) < 100) {
      data2 <- particles_data_join[
        !(particles_data_join$particle_id %in% particles_use),
      ]

      data2 <- data2[order(data2$overlap_logit, decreasing = TRUE), ]
      particles_more <- data2$particle_id[1:200]
      particles_use <- c(particles_use, particles_more)
    }

    rows_fit <- which(particles_data_join$particle_id %in% particles_use)
  }

  # Fit GP
  gp_data <- particles_data_join[rows_fit, ]

  gp_result <- gp_fit(data = gp_data, response = response,
                      gp_optimize = gp_optimize)

  # get threshold
  loglik_threshold_new <- get_loglik_threshold(gp_result$loglik_high,
                                               loglik_tolerance,
                                               loglik_lower)
  # Tidy objects for return
  if(is.null(previous_wave)) {
    gp_list <- NULL
    loglik_thresholds <- NULL
  }

  result <- list(
    gps = c(gp_list, gp_result$loglik_model),
    in_model = in_model,
    loglik_optim = if(gp_optimize) gp_result$op else NULL,
    loglik_thresholds = c(loglik_thresholds, loglik_threshold_new),
    particles_data = particles_data_join,
    response = response
  )

  return(result)
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



# Function to run a few waves
waves_runner <- function(n_waves = 5,
                         loglik_tolerance = 3,
                         loglik_lower = NULL,
                         n_pw = 1000,
                         n_fit = 800,
                         prop_best = 0.3,
                         metric = NULL,
                         rep_id = NULL) {

  if(n_waves < 2) stop("Use this function to fit more than 1 wave.")

  waves_list <- vector("list", n_waves)

  # initialize first wave
  message(paste("Metric: ", metric, ", rep_id: ", rep_id, ", wave 1/",
                n_waves, sep = ""))

  waves_list[[1]] <- loglik_update_cumulative(
    loglik_tolerance = loglik_tolerance,
    loglik_lower = loglik_lower,
    n_pw = n_pw,
    n_fit = n_fit,
    prop_best = prop_best,
    metric = metric,
    rep_id = rep_id,
    gp_optimize = FALSE
  )

  # run the remaining
  for(w in 2:n_waves) {

    message(paste("Metric: ", metric, ", rep_id: ", rep_id, ", wave ", w,
                  "/", n_waves, sep = ""))

    do_optim <- ifelse(w == n_waves, TRUE, FALSE)

    waves_list[[w]] <- loglik_update_cumulative(
      previous_wave = waves_list[[w-1]],
      loglik_tolerance = loglik_tolerance,
      loglik_lower = loglik_lower,
      n_pw = n_pw,
      n_fit = n_fit,
      prop_best = prop_best,
      metric = metric,
      rep_id = rep_id,
      gp_optimize = do_optim
    )
  }

  last_wave <- waves_list[[n_waves]]
  rm(waves_list)
  gc()
  return(last_wave)
}

# Long sobol sequence -----------------------------------------------------

# for all estimations to run over the same particles (when possible), a very long
# sobol sequence will be produced, expecting not to need further "simulations".
# A larger sequence will be created from scratch if more particles are needed,
# to ensure the same particles are always used.
n_large <- 300000
particles_all <- particles_sim(n_large)

# Testing -----------------------------------------------------------------

bb <- get_bounds(fire_data = spread_data)

w1 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = NULL,
  loglik_tolerance = 0.03,
  loglik_lower = NULL,
  n_pw = 800,
  n_fit = 800,
  prop_best = 0.3,
  gp_optimize = TRUE,
  response = "overlap"
)
p1 <- loglik_plot(w1)

w2 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w1,
  n_pw = 800,
  n_fit = 800,
  loglik_tolerance = 0.03
)
p2 <- loglik_plot(w2)

w3 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w2,
  n_pw = 800,
  n_fit = 800,
  loglik_tolerance = 0.03
)
p3 <- loglik_plot(w3)


w3 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w2,
  n_pw = 800,
  n_fit = 800,
  loglik_tolerance = 0.03
)
p3 <- loglik_plot(w3)

w4 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w3,
  n_pw = 800,
  n_fit = 1200,
  loglik_tolerance = 0.03
)
p4 <- loglik_plot(w4) # hermoso
# mejora, y a medida que se pone más puntudo, tarda más en buscar nuevas
# partículas.

w5 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w4,
  n_pw = 800,
  n_fit = 1200,
  loglik_tolerance = 0.03
)
p5 <- loglik_plot(w5)


w6 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w5,
  n_pw = 800,
  n_fit = 800,
  loglik_tolerance = 0.03
)
p6 <- loglik_plot(w6)



# que tarde mucho en buscar nuevas partículas puede deberse a que
# la previa sea demasiado amplia, entonces intenta con un montón
# de valores que no tienen sentido, y solo algunos funcionan.
# Entonces, quizás no sea necesario poner una previa taaaan amplia.
# Al fin y al cabo, si hacemos rejection sampling, no sería necesario
# poder describir tan bien esa región en donde el overlap se hace casi plano
# y cercano a cero.

# bueno, hay que mejorar la busqueda de nuevas partículas porque tarda mucho más
# eso que simular fuegos.

# ¿Cómo se podría restringir?
# En base a las partículas que pasaron el "test", definir el rango de una
# nueva secuencia de sobol, pero que sea plana. Es decir, obtenemos los
# límites de los "datos" buenos. Esos pueden ser los que sobrevivieron a la última
# ola, o, si son muy pocos, los merjores 100. Obtenemos el rango como
# apply(par_values, 2, range)
# y lo ampliamos un poquito:
# lower <- range[1] - diff(range) * 0.25
# upper <- range[2] - diff(range) * 0.25
# (Para los params truncados, hacemos lower = ifelse(lower < 0, 0, lower)).
# y ahí creamos una nueva sequencia de sobol en la que buscamos nuevas partículas.

# Es impresionante lo que tarda esta garcha.

# Y no debería ser problema que se ponga muy puntudo,

# saveRDS(w8, file.path("files", "pseudolikelihood_estimation", "gp_waves-2008_3.rds"))


# Running waves -----------------------------------------------------------

# results_list <- vector("list", n_rep)
# names(results_list) <- 1:n_rep
# for(i in 1:n_rep) {
#   tmp <- vector("list", n_met)
#   names(tmp) <- metric_names
#   results_list[[i]] <- tmp
# }
# for(r in 1:n_rep) {
#   for(m in metric_names) {
#     results_list[[r]][[m]] <- waves_runner(
#       n_waves = 5,
#       loglik_tolerance = 3,
#       loglik_lower = NULL,
#       n_pw = 1000,
#       n_fit = 800,
#       prop_best = 0.3,
#       metric = m,
#       rep_id = r
#     )
#     filename <- paste("files/similarity_selection_result_rep-", r,
#                       "-", m, ".rds", sep = "")
#     saveRDS(results_list[[r]][[m]], filename)
#   }
# }

# f0 <- list.files(file.path("files"))
# f1 <- f0[grep("similarity_selection_result_", f0)]
# f2 <- strsplit(f1, split = "[-|.]")
# freps <- sapply(f2, function(i) i[2]) %>% as.numeric
# fmets <- sapply(f2, function(i) i[3])
# for(i in 1:length(freps)) {
#   r <- freps[i]
#   m <- fmets[i]
#   results_list[[r]][[m]] <- readRDS(file.path("files", f1[i]))
# }
# saveRDS(results_list, file.path("files", "similarity_selection_results.rds"))
results_list <- readRDS(file.path("files", "similarity_selection_results.rds"))
