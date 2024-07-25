# SMC notes ---------------------------------------------------------------

# acá heredo mucho código de los GP, porque aún no sabemos cuán útil será
# cambiar a SMC.
# Y en este script agrego el parámetro steps al modelo.

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
upper_limit <- 1
n_coef <- 6
par_names <- c("intercept", "vfi", "tfi", "slope", "wind", "steps")

# number of particles to choose by wave
n_pw <- 1000

# number of fires to simulate by particle
n_sim <- 20

# formula for in_bounds model (gam)
bounds_model_formula <- formula(
  in_bounds ~
    s(intercept, k = 8, bs = "cr") +
    s(vfi, k = 8, bs = "cr") +
    s(tfi, k = 8, bs = "cr") +
    s(slope, k = 8, bs = "cr") +
    s(wind, k = 8, bs = "cr") +
    s(steps, k = 8, bs = "cr")
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
                       steps_lower = 10,
                       steps_upper = 1000,
                       type = "rng", # or "quantile"
                       n = 1,
                       p = NULL) {

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
}

# function to make particles from a sobol sequence at the prior distribution.
particles_sim_prior <- function(N = 100, d = n_coef, sobol_init = FALSE,
                                mu_int = 0, sd_int = 20,
                                r_slope = 0.04,
                                r_fi = 0.15,
                                r_wind = 0.3,
                                steps_lower = 10,
                                steps_upper = 1000) {
  prior_dist(type = "quantile",
             p = sobol(N, dim = d, init = sobol_init),
             mu_int = mu_int, sd_int = sd_int,
             r_slope = r_slope,
             r_fi = r_fi,
             r_wind = r_wind,
             steps_lower = steps_lower,
             steps_upper = steps_upper)
}

# function to make a sobol sequence on a compact square distribution, defined
# by the bounds matrix. The matrix must have 2 rows (lower-upper) and d columns.
particles_sim_box <- function(N = 100, d = n_coef, sobol_init = FALSE,
                              bounds = NULL) {

  p <- sobol(N, dim = d, init = sobol_init)
  pbox <- p

  for(c in 1:d) {
    pbox[, c] <- qunif(p[, c], min = bounds["lower", c], max = bounds["upper", c])
  }

  colnames(pbox) <- par_names

  return(pbox)
}

# function factory: it creates a function that returns the sobol sequence on a
# box. This is to be passed to more_particles()
particles_sim_box_factory <- function(bounds = NULL) {

  sim_box <- function(N = 100, d = n_coef, sobol_init = FALSE) {
    particles_sim_box(N = N, d = d, sobol_init = sobol_init, bounds = bounds)
  }

  return(sim_box)
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
# Returns a binary vector indicating which rows of the particles_mat should be
# considered.
simulate_in_bounds <- function(in_model, particles_mat) {

  # get maximum fitted probability to relativize predictions.
  pmax <- predict(in_model, type = "response") %>% max

  # predict in_range_prob
  prob <- predict(in_model, as.data.frame(particles_mat),
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

# Function to filter particles based on the previously fitted gps.
# It returns a vector with the row index of the provided particles_mat.
gp_filter <- function(gp_list = NULL, loglik_thresholds = NULL,
                      particles_mat = NULL) {

  # compute predictions for all particles at all GPs
  preds <- lapply(gp_list, function(x) {
    x$pred(particles_mat, se.fit = T)
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

  return(keep_these)
}

# Search for more particles to simulate, returning a parameters matrix
# (particles_mat).
# target_dist_prng stands for "pseudo random number generator of the target
# distribution". It's a function that generates a sobol sequence over the
# desired parameter space. It could use the prior or a compact constrained
# space. It must have arguments N and init.
# sobol_init indicates whether the sobol sequence should be initialized.
# It should be set to TRUE when a never-used target_dist is to be used.
more_particles <- function(gp_list = NULL,
                           loglik_thresholds = NULL,
                           in_model = NULL,
                           n_pw = 1000,
                           target_dist_prng = NULL,
                           sobol_init = T
                           ) {

  if(!is.null(gp_list)) {

    added <- 0
    new_particles <- NULL
    cycle <- 1

    while(added < n_pw) {

      n_try <- n_pw * (1 + length(gp_list))

      sobol_init <- ifelse(cycle == 1, sobol_init, FALSE)
      particles_try <- target_dist_prng(N = n_try, sobol_init = sobol_init)

      # increase cycle step befor jumps
      cycle <- cycle + 1

      # evaluate in_bounds
      in_range_ids <- which(simulate_in_bounds(in_model, particles_try) == 1)

      if(length(in_range_ids) == 0) next

      # keep only in_range
      particles_try <- particles_try[in_range_ids, , drop = F]

      # evaluate with GPs. This returns a vector with row ids of passing particles.
      high_ll_ids <- gp_filter(gp_list = gp_list,
                               loglik_thresholds = loglik_thresholds,
                               particles_mat = particles_try)

      if(length(high_ll_ids) == 0) next

      # keep only good ones
      particles_keep <- particles_try[high_ll_ids, , drop = F]

      # merge with previous new particles
      new_particles <- rbind(new_particles, particles_keep)
      added <- nrow(new_particles)
    } # end while

  } else {
    # if it's the first wave
    new_particles <- particles_sim_prior(N = n_pw, sobol_init = T)
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


# box_bounds: function to get the parameter bounds over which to define a
# compact search space. This is used to constrain the search space with respect
# to the prior distribution. Hence, it takes the best particles to get the
# parameters ranges.

# n_best: how many of the best particles are to be used to define parameters bounds?
#   if NULL, all the particles from the latest wave are used (n_last), but if
#   n_last < n_min, the best n_min are used.
# expand_lower-upper: range proportion in which the observed bounds are expanded,
#   as +- diff(range(x)) * expand. bounded parameters are forced to have lower
#   limit at zero.
# in_bounds: consider only in_bounds particles?
box_bounds <- function(data, n_best = NULL, n_min = 50,
                       expand_lower = 0.1, expand_upper = 0.1,
                       in_bounds = F) {

  if(in_bounds) data <- data[data$in_bounds == 1, ]

  # subset data to compute range
  if(is.null(n_best)) {
    last_wave_ids <- which(data$wave == max(data$wave))
    n_last <- length(last_wave_ids)
    if(n_last >= n_min) {
      par_use <- data$par_values[last_wave_ids, ]
    } else {
      data_ord <- data[order(data$overlap, decreasing = T), ]
      par_use <- data_ord$par_values[1:n_min, ]
    }
  } else {
    n_best <- min(n_best, nrow(data))
    data_ord <- data[order(data$overlap, decreasing = T), ]
    par_use <- data_ord$par_values[1:n_best, ]
  }

  # compute ranges and expand them
  lowers_obs <- apply(par_use, 2, min)
  uppers_obs <- apply(par_use, 2, max)
  width <- uppers_obs - lowers_obs

  lowers_wide <- lowers_obs - width * expand_lower
  uppers_wide <- uppers_obs - width * expand_upper

  # contrain minimums to zero
  lowers_wide[-1] <- ifelse(lowers_wide[-1] < 0, 0, lowers_wide[-1])

  # matrix to return
  b <- rbind(lowers_wide, uppers_wide)
  colnames(b) <- colnames(data$par_values)
  rownames(b) <- c("lower", "upper")

  return(b)
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
# target_prior: should the prior be used to get new particles? If FALSE, uses a
#   search over a compact range of the parameters (the range of the best
#   particles found up to the current wave, see particles_sim_box).
# gp_optimize = should the fitted gp be optimized? Not recommended
#   if it's first fitted with many particles out of bounds. Useful for plotting
#   the MLE.

loglik_update_cumulative <- function(
    fire_data = NULL,
    similarity_bounds = NULL,
    previous_wave = NULL,
    loglik_tolerance = 0.03,
    loglik_lower = NULL,
    n_pw = 800,
    n_fit = 800,
    prop_best = 0.3,
    target_prior = TRUE,
    gp_optimize = TRUE,
    use_in_bounds = FALSE,
    n_sim = 20,
    response = "overlap",
    # arguments for box_bounds, needed when target_prior = F
    box_n_best = NULL, box_n_min = 50,
    box_expand_lower = 0.1, box_expand_upper = 0.1,
    box_in_bounds = F
) {

  #### testo

  # for wave 1

  # fire_data = spread_data
  # similarity_bounds = bb
  # previous_wave = w1
  # n_pw = 800
  # n_fit = 800
  # loglik_tolerance = 0.03
  # target_prior = F         # better search
  # box_n_min = 100

  # for wave 2

  # fire_data = spread_data
  # similarity_bounds = bb
  # previous_wave = w1
  # loglik_tolerance = 0.03
  # loglik_lower = NULL
  # n_pw = 800
  # n_fit = 800
  # prop_best = 0.3
  # target_prior = F
  # gp_optimize = TRUE
  # use_in_bounds = FALSE
  # n_sim = 20
  # response = "overlap"
  ## arguments for box_bounds needed when target_prior = F
  # box_n_best = NULL
  # box_n_min = 100
  # box_expand_lower = 0.1
  # box_expand_upper = 0.1
  # box_in_bounds = F

  ####

  # when a previous wave has been run
  if(!is.null(previous_wave)) {
    gp_list <- previous_wave$gps
    wave_last <- length(gp_list)    # number of previous waves
    gp_last <- gp_list[[wave_last]]
    message(paste("Wave", 1 + wave_last))

    # The previous wave brings K GPs, K loglik_highs, and K-1 loglik_thresholds.
    # This is aimed to allow the use of a loglik_tolerance before simulating
    # new fires. Hence, the new threshold for the last GP is computed now.

    # get loglik_highs to compute the latest threshold
    loglik_highs <- previous_wave$loglik_highs
    loglik_high_last <- loglik_highs[length(loglik_highs)]

    # get threshold
    loglik_threshold_new <- get_loglik_threshold(loglik_high_last,
                                                 loglik_tolerance,
                                                 loglik_lower)

    # merge with old thresholds
    loglik_thresholds <- c(previous_wave$loglik_thresholds,
                           loglik_threshold_new)

    in_model <- previous_wave$in_model
    particles_data <- previous_wave$particles_data
    # remove row_id
    part_col <- which(names(particles_data) == "particle_id")
    particles_data <- particles_data[, -part_col]

    # filter simulated particles according to last gp (this is the
    # post-simulation filter)
    message("Filtering old particles")
    rows_eval <- which(particles_data$wave >= (wave_last - 1))
    ids_keep <- gp_filter(list(gp_last),
                          loglik_thresholds = loglik_threshold_new,
                          particles_mat = particles_data$par_values[rows_eval, ])

    if(length(ids_keep) > 0) {
      rows_keep <- rows_eval[ids_keep]
      particles_data$wave[rows_keep] <- wave_last + 1
    }

    # get new particles meeting the same criterion
    message("Getting more particles")

    if(target_prior) {
      prng <- particles_sim_prior
      sobol_init <- F # because there is a previous wave
    } else {
      limits <- box_bounds(
        particles_data, n_best = box_n_best, n_min = box_n_min,
        expand_lower = box_expand_lower, expand_upper = box_expand_upper,
        in_bounds = box_in_bounds
      )
      prng <- particles_sim_box_factory(bounds = limits)
      sobol_init <- T
    }

    new_particles <- more_particles(gp_list = gp_list,
                                    loglik_threshold = loglik_thresholds,
                                    in_model = in_model,
                                    n_pw = n_pw,
                                    target_dist_prng = prng,
                                    sobol_init = sobol_init)

  } else {
    message("Wave 1")
    message("Getting more particles")
    new_particles <- more_particles(n_pw = n_pw,
                                    target_dist_prng = particles_prior_sim,
                                    sobol_init = T) # because it's the first wave!
  }

  # Define wave number
  this_wave <- ifelse(is.null(previous_wave),
                      1,
                      wave_last + 1)

  # Simulate loglik on new particles
  message("Simulating fires")

  particles_data_new <- similarity_simulate_parallel(
    particles_mat = new_particles,
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

    # create row id
    particles_data_join$particle_id <- 1:nrow(particles_data_join)

    # subset
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

      data2 <- data2[order(data2$overlap, decreasing = TRUE), ]
      particles_more <- data2$particle_id[1:200]
      particles_use <- c(particles_use, particles_more)
    }

    rows_fit <- which(particles_data_join$particle_id %in% particles_use)
  }

  # Fit GP
  gp_data <- particles_data_join[rows_fit, ]

  gp_result <- gp_fit(data = gp_data, response = response,
                      gp_optimize = gp_optimize)

  # Tidy objects for return
  if(is.null(previous_wave)) {
    gp_list <- NULL
    loglik_thresholds <- NULL
    loglik_highs <- NULL
  }

  result <- list(
    gps = c(gp_list, gp_result$loglik_model),
    in_model = in_model,
    loglik_optim = if(gp_optimize) gp_result$op else NULL,
    loglik_thresholds = loglik_thresholds,
    loglik_highs = c(loglik_highs, gp_result$loglik_high),
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

wave_plot <- function(data, alpha = 0.3) {
  par(mfrow = c(3, 2))
  for(i in 1:n_coef) {
    plot(data$overlap ~ data$par_values[, i], ylab = "overlap",
         xlab = par_names[i],
         pch = 19, col = rgb(0, 0, 0, alpha))
  }
  par(mfrow = c(1, 1))
}


# sequential monte carlo wave. It takes data from a previous set of simulated
# particles, reproduces them according to their overlap only and returns a new
# set of simulated particles.
# focal_sampling resamples the good particles. If FALSE, the new samples are
# generated from a global MVN distribution. In that case, it is convenient
# to set lowlik_lower relatively low.
smc <- function(data, n = 100, shrink_factor = 0.2,
                loglik_lower = 0.2,
                n_min = 50,
                n_best = NULL,
                rescale_prob = TRUE,
                focal_sampling = TRUE) {

  # filter the good data points
  dgood <- data[data$overlap >= loglik_lower, ]
  n_good <- nrow(dgood)
  if(n_good < n_min) {
    data2 <- data[order(data$overlap, decreasing = TRUE), ]
    dgood <- data2[1:n_min, ]
  }

  # resample the good particles (to increase the selectivity)
  if(rescale_prob) {
    p <- dgood$overlap - min(dgood$overlap)
  } else {
    p <- dgood$overlap
  }

  ids_rep <- sample(1:nrow(dgood), size = n, replace = T,
                    prob = p)

  # expanded dataset
  dexp <- dgood[ids_rep, ]

  # log the positives
  dexp$par_values_log <- dexp$par_values
  dexp$par_values_log[, 2:n_coef] <- log(dexp$par_values_log[, 2:n_coef])

  mvn_fitted <- FitGMM(dexp$par_values_log)

  # sample new particles using shrunk vcov
  V_shrunk <- mvn_fitted@Covariance * shrink_factor ^ 2
  # (shrink factor is at the sigma scale)
  candidates_log <- mgcv::rmvn(n = n,
                               mu = dexp$par_values_log, V = V_shrunk)

  colnames(candidates_log) <- par_names

  # unlog the coef
  candidates <- candidates_log
  candidates[, 2:n_coef] <- exp(candidates_log[, 2:n_coef])

  # message("Simulating fires")
  sim_result <- similarity_simulate_parallel(particles_mat = candidates,
                                             fire_data = spread_data)

  # define wave
  sim_result$wave <- max(dgood$wave) + 1

  # merge old and new datasets
  res <- rbind(
    data[, colnames(sim_result)],
    sim_result
  )

  return(res)
}


# now we have an iterator for the smc() function. It's just a loop.
smc_iterate <- function(n_waves = 10,
                        data, n = 100,
                        shrink_factor = 0.2,
                        loglik_lower = 0.2,
                        n_min = 50,
                        rescale_prob = TRUE) {

  res <- smc(data = data, n = n,
             shrink_factor = shrink_factor,
             loglik_lower = loglik_lower,
             n_min = n_min,
             rescale_prob = rescale_prob)
  m <- paste("SMC iteration 1; max overlap = ",
             round(max(data$overlap), 4), sep = "")
  message(m)

  for(i in 2:n_waves) {
    res <- smc(data = res, n = n,
               shrink_factor = shrink_factor,
               loglik_lower = loglik_lower,
               n_min = n_min,
               rescale_prob = rescale_prob)

    m <- paste("SMC iteration ", i, "; max overlap = ",
               round(max(res$overlap), 4), sep = "")
    message(m)
  }

  return(res)
}


# SMC from scratch --------------------------------------------------------

# initial settings

bb <- get_bounds(fire_data = spread_data)
p0 <- particles_sim_prior(N = 600, sobol_init = T,
                          sd_int = 5,
                          steps_lower = 10,
                          steps_upper = bb["largest", "steps_used"] * 2)


# w1
w1 <- similarity_simulate_parallel(p0, spread_data)
wave_plot(w1)

# w2
smc2 <- smc_iterate(n_waves = 6, data = w1)

# w3
smc3 <- smc_iterate(n_waves = 6, data = smc2, shrink_factor = 0.5)
smc3 <- smc_iterate(n_waves = 6, data = smc3, shrink_factor = 1)

wave_plot(smc3, alpha = 0.1)
pairs(smc3$par_values[smc3$overlap > 0.2, ],
      pch = 19, col = rgb(0, 0, 0, 0.05))
pairs(smc3$par_values,
      pch = 19, col = rgb(0, 0, 0, 0.05))

saveRDS(smc3, "files/pseudolikelihood_estimation/smc_waves_2008_3.rds")
saveRDS(bb, "files/pseudolikelihood_estimation/bounds_2008_3.rds")
# MVN from overlap (resampling) -------------------------------------------

dgood <- smc3[smc3$overlap >= 0.2, ]
prs <- dgood$overlap # not rescale
# resample the good particles
n_expand <- 20000
ids_rep <- sample(1:nrow(dgood), size = n_expand, replace = T, prob = prs)
dexp <- dgood[ids_rep, ]

# log the positives
dexp$par_values_log <- dexp$par_values
dexp$par_values_log[, 2:n_coef] <- log(dexp$par_values_log[, 2:n_coef])

# fit a skew-mvn to these fake data
f2 <- sn::selm.fit(
  x = matrix(rep(1, nrow(dexp)), ncol = 1),
  y = dexp$par_values_log,
  family = "SN"
)
f2$param
d2 <- sn::makeSECdistr(dp = f2$param$dp, family = "SN", compNames = par_names)
plot(d2)
r2 <- sn::rmsn(n = 3e4, dp = f2$param$dp)

par(mfrow = c(3, 2))
for(i in 1:n_coef) {
  plot(density(dexp$par_values_log[, i]), main = par_names[i])
  # lines(density(r1[, i]), col = 2)
  lines(density(r2[, i]), col = 3)
}
par(mfrow = c(1, 1))


# unlog
r2ori <- r2
r2ori[, 2:n_coef] <- exp(r2[, 2:n_coef])

par(mfrow = c(3, 2))
for(i in 1:n_coef) {
  plot(density(dexp$par_values[, i]), main = par_names[i])
  # lines(density(r1ori[, i]), col = 2)
  lines(density(r2ori[, i]), col = "red")
}
par(mfrow = c(1, 1))

summary(dgood$par_values)

# qqplot

r3 <- sn::rmsn(n = nrow(dexp), dp = f2$param$dp)
r3ori <- r3
r3ori[, 2:n_coef] <- exp(r3[, 2:n_coef])

par(mfrow = c(3, 2))
for(i in 1:n_coef) {
  # i = 1
  q_aprox <- quantile(r3ori[, i], probs = ppoints(n = 200))
  q_obs <- quantile(dexp$par_values[, i], probs = ppoints(n = 200))
  plot(q_obs ~ q_aprox, pch = 19, col = rgb(0, 0, 0, 0.5),
       main = par_names[i])
  abline(0, 1)
  abline(h = quantile(dexp$par_values[, i], probs = c(0.025, 0.975)),
         lty = 2, col = "red")
  abline(h = quantile(dexp$par_values[, i], probs = c(0.01, 0.99)),
         lty = 2, col = "blue")
}
par(mfrow = c(1, 1))


# Acá hemos encontrado un ajuste bimodal. Está complicado trabajar con eso.
# En parte es entendible que haya corr entre intercept y tfi, ya que el fuego
# es peque y eso, pero la verdad que no se explica por qué el overlap
# es tan alto en dos zonas de intercept.

# la bimodalidad la da el viento, con sus partículas > 6.
plot(dgood$edge_bin_prop ~ dgood$par_values[, "wind"])
plot(dgood$edge ~ dgood$par_values[, "wind"])
plot(dgood$par_values[, "steps"] ~ dgood$par_values[, "wind"])

# Seguir jugando con esto. Probar si "muestreando" la verdadera likelihood,
# sin sesgos, se llega a algo mejor, aunque no creo que cure esa bimodalidad.

# Luego chequear si muestreando una density con incertidumbre se llega a la
# misma density.



# kdevine para estimar likelihood desde muestras --------------------------------

library(kdevine)
# estimate density (use xmin to indicate positive support)
fit <- kdevine(dgood$par_values, xmin = c(-Inf, rep(0, n_coef - 1)))
# evaluate density estimate
fitted_density <- dkdevine(dgood$par_values, fit)
plot(fitted_density ~ dgood$overlap)
plot(fit)

draw <- rkdevine(3000, fit)

for(i in 1:n_coef) {
  plot(density(dexp$par_values[, i]), main = par_names[i])
  # lines(density(r1ori[, i]), col = 2)
  lines(density(r2ori[, i]), col = "red")
  lines(density(draw[, i]), col = "blue")
}
# bueno, parece que anda bien.
# Luego muestrear la likelihood como una verdadera distrib.
