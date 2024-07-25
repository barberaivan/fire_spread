# Here I run waves for 2008_3 at the overlap scale and correcting the bugs.
# Antes sd_int era 5, ahora lo puse en 20.


# Notes for efficient search of new particles -----------------------------

# Acá implemento una forma más eficiente de buscar nuevas partículas,
# en base a lo que vi en el archivo anterior (_02).

# Cuando se encuentra una pequeña zona con alto overlap, el GP se pone muy
# filoso, y la búsqueda de nuevas partículas se hace muy lenta. Esto se debe
# a que la previa es demasiado amplia para la escala del GP, entonces casi todas
# las partículas son rechazadas. Además, evaluar el GP no es de lo más rápido,
# menos si se estima con muchos datos.

# Para reducir este problema, a medida que el GP avanza, se restringe el espacio
# parámetros en donde se proponen nuevas partículas. En vez de ser una secuencia
# de sobol sobre la previa, se usa un espacio compacto rectangurar, es decir,
# una uniforme acotada. Los límites de este espacio restringido se definen
# algo así:

# Cuando comienza una ola, se evalúa qué partículas pasan a la siguiente ola.
# Si son 50 o más, se usan solo esas, y si son < 50, se eligen las mejores 50
# y listo.

# Obtenemos el rango como
# apply(par_values, 2, range)
# y lo ampliamos un poquito:
# lower <- range[1] - diff(range) * widen
# upper <- range[2] + diff(range) * widen,
# con widen ~ 0.25, a tunear.

# Para los params truncados, además hacemos
# lower <- ifelse(lower < 0, 0, lower)).

# Y definimos la nueva sequencia de sobol para uniformes en [lower, upper].

# Esto implica que ya no se use una larga secuencia de sobol incial de la cual
# vamos sacando partículas para simular. En cada búsqueda directamente generamos
# una secuencia y devolvemos esos parámetros para ser simulados.
# O sea que more_particles no debe devolver ids de partículas, solo una matriz
# de parámetros.


# Notes for better use of loglik_threshold --------------------------------

# En este script también modifico la forma en que se filtran las partículas
# para ver si sobreviven a la siguiente ola. La idea es que cuando uno corre
# una ola pueda definir el umbral que se usará en esa misma ola. Para eso, los
# pasos serían así:

# wave 0:
#   it does not require loglik_tolerance.
#   returns list with a GP with its respective loglik_high.
# wave 1:
#   may use a loglik_tolerance.
#   with that tolerance, it computes the loglik_threshold to filter new particles.
#   return a list with
#     2 GPs,
#     2 loglik_high,
#     1 loglik_threshold (associated to the first GP)

# Hence, the wave result will always have K GPs and K loglik_highs, but
# only K-1 loglik_thresholds, because the wave finishes by fitting the GP.
# The use of the thresholds is after that, in the following wave.

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

# function to make particles from a sobol sequence at the prior distribution.
particles_sim_prior <- function(N = 100, d = n_coef, sobol_init = FALSE) {
  prior_dist(type = "quantile",
             p = sobol(N, dim = d, init = sobol_init))
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

# Testing -----------------------------------------------------------------

bb <- get_bounds(fire_data = spread_data)

w1 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  n_pw = 800,
  n_fit = 800,
  prop_best = 0.3,
  target_prior = TRUE,
  gp_optimize = TRUE,
  # arguments for box_bounds, needed when target_prior = F
  box_n_best = NULL, box_n_min = 50,
  box_expand_lower = 0.1, box_expand_upper = 0.1,
  box_in_bounds = F
)
p1 <- loglik_plot(w1)

w2 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w1,
  n_pw = 800,
  n_fit = 800,
  loglik_tolerance = 0.05,
  target_prior = F,         # better search
  box_n_min = 100
)
p2 <- loglik_plot(w2)

w3 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w2,
  n_pw = 800,
  n_fit = 800,
  loglik_tolerance = 0.03,  # con 0.02 rechaza todo
  target_prior = F,
  box_n_min = 50
)
p3 <- loglik_plot(w3)

w4 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w3,
  n_pw = 800,
  n_fit = 800,
  loglik_tolerance = 0.03,  # con 0.02 rechaza todo
  target_prior = F,
  box_n_best = 50           # lowering this speedsup the search
)
p4 <- loglik_plot(w4)


w5 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w4,
  n_pw = 800,
  n_fit = 800,
  loglik_tolerance = 0.03,  # con 0.02 rechaza todo
  target_prior = F,
  box_n_best = 30           # lowering this speedsup the search
)
p5 <- loglik_plot(w5)

w6 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w5,
  n_pw = 800,
  n_fit = 800,
  loglik_tolerance = 0.03,  # con 0.02 rechaza todo
  target_prior = F,
  box_n_best = 20           # lowering this speedsup the search
)
p6 <- loglik_plot(w6)

# the search is much faster now, even with peaks in the GPs

w7 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w6,
  n_pw = 800,
  n_fit = 800,
  loglik_tolerance = 0.03,  # con 0.02 rechaza todo
  target_prior = F,
  box_n_best = 10           # lowering this speedsup the search
)
p7 <- loglik_plot(w7)

w8 <- loglik_update_cumulative(
  fire_data = spread_data,
  similarity_bounds = bb,
  previous_wave = w7,
  n_pw = 800,
  n_fit = 800,
  loglik_tolerance = 0.02,  # con 0.02 rechaza todo
  target_prior = F,
  box_n_best = 10           # lowering this speedsup the search
)
p8 <- loglik_plot(w8, color_point = "wave_plot")



# SMC ---------------------------------------------------------------------
# (o algo así)
# library(Rfast) # to fit a mvn dist
library(MGMM)

dgood <- w8$particles_data[w8$particles_data$overlap >= 0.2, ]
nrow(dgood) # 65, re pocos.

# resample the good particles
n_expand <- 500
ids_rep <- sample(1:nrow(dgood), size = n_expand, replace = T, prob = dgood$overlap)
dexp <- dgood[ids_rep, ]

# log the positives
dexp$par_values_log <- dexp$par_values
dexp$par_values_log[, 2:n_coef] <- log(dexp$par_values_log[, 2:n_coef])

pairs(dexp$par_values)
pairs(dexp$par_values_log)
# sí, hay correlation

f1 <- FitGMM(dexp$par_values_log)
print(f1)
r1 <- rGMM(n_expand, d = n_coef, k = 1, means = f1@Mean, cov = f1@Covariance)
# ajustar una mvn??
pairs(r1)
pairs(dexp$par_values_log, add = TRUE, col = 2)

for(i in 1:n_coef) {
  plot(density(dexp$par_values_log[, i]), main = par_names[i])
  lines(density(r1[, i]), col = 2)
}
# the normal approximation is not so bad.

# use the correlation matrix, but decrease sigma to 1 / 10 of the fitted values,
# for a more local approximation

cor1 <- cov2cor(f1@Covariance)
sigma <- sqrt(diag(f1@Covariance))
shrink_factor <- 0.25
sigma_shrunk <- sigma * shrink_factor
V_shrunk <- diag(sigma_shrunk) %*% cor1 %*% diag(sigma_shrunk)

# sample new particles using that vcov
candidates_log <- mgcv::rmvn(n = n_expand, mu = dexp$par_values_log, V = V_shrunk)

par(mfrow = c(3, 2))
for(i in 1:n_coef) {
  plot(density(dexp$par_values_log[, i]), main = par_names[i])
  lines(density(candidates_log[, i]), col = 2)
}
par(mfrow = c(1, 1))

colnames(candidates_log) <- par_names

# unlog the coef
candidates <- candidates_log
candidates[, 2:n_coef] <- exp(candidates_log[, 2:n_coef])

smc1 <- similarity_simulate_parallel(particles_mat = candidates,
                                     fire_data = spread_data)

summary(smc1$overlap)

summary(dgood$overlap)


# focal_sampling resamples the good particles. If false, the new samples are
# generated from a global MVN distribution. In that case, it is convenient
# to set lowlik_lower relatively low could be focal (reproduce particles)
# or global, when an MVN is fitted to all the good particles and samples
# from this global distribution are produced.
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

  ids_rep <- sample(1:nrow(dgood), size = n_expand, replace = T,
                    prob = p)

  # expanded dataset
  dexp <- dgood[ids_rep, ]

  # log the positives
  dexp$par_values_log <- dexp$par_values
  dexp$par_values_log[, 2:n_coef] <- log(dexp$par_values_log[, 2:n_coef])

  mvn_fitted <- FitGMM(dexp$par_values_log)

  # sample new particles using shrunk vcov
  V_shrunk <- f1@Covariance * shrink_factor ^ 2
  # (shrink factor is at the sigma scale)
  candidates_log <- mgcv::rmvn(n = n_expand,
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


smc1 <- smc(w8$particles_data); max(smc1$overlap)
smc2 <- smc(smc1); max(smc2$overlap)
smc3 <- smc(smc2); max(smc3$overlap)
smc4 <- smc(smc3); max(smc4$overlap)
smc5 <- smc(smc4); max(smc5$overlap)
smc6 <- smc(smc5); max(smc6$overlap)

smc7 <- smc_iterate(data = smc6)
nrow(smc7[smc7$overlap >= 0.23, ])
smc8 <- smc_iterate(data = smc7, shrink_factor = 0.3, loglik_lower = 0.23)
smc9 <- smc_iterate(data = smc8, shrink_factor = 0.45, loglik_lower = 0.23,
                    n_waves = 10)

smc10 <- smc_iterate(data = smc9, loglik_lower = 0.23)

plot(ecdf(smc9$overlap[smc9$overlap >= 0.2]))

par(mfrow = c(2, 3))
for(i in 1:n_coef) {
  plot(smc10$overlap ~ smc10$par_values[, i], pch = 19,
       col = rgb(0, 0, 0, 0.1), xlim = range(smc10$par_values[, i]),
       xlab = par_names[i], ylab = "overlap")
}
par(mfrow = c(1, 1))

nrow(smc10[smc10$overlap >= 0.25, ]) # 105
smc11 <- smc_iterate(data = smc10, loglik_lower = 0.25)
nrow(smc11[smc11$overlap >= 0.25, ]) # 135
smc12 <- smc_iterate(data = smc11, loglik_lower = 0.25, shrink_factor = 0.5)

nrow(smc12[smc12$overlap >= 0.20, ]) # 3256

par(mfrow = c(3, 2))
for(i in 1:n_coef) {

  if(i == 1) {
    rr <- c(-20, 30)
  } else {
    rr <- range(smc12$par_values[, i])
  }

  plot(smc12$overlap ~ smc12$par_values[, i], pch = 19,
       col = rgb(0, 0, 0, 0.03), xlim = rr,
       xlab = par_names[i], ylab = "overlap")
}
par(mfrow = c(1, 1))


# compare densities of truncated data with estimated MVN densities
dgood <- smc12[smc12$overlap >= 0.2, ]

prs <- dbest$overlap - 0.2 + 1e-4 # rescale overlap
prs <- dbest$overlap # not rescale
# resample the good particles
n_expand <- 20000
ids_rep <- sample(1:nrow(dgood), size = n_expand, replace = T, prob = prs)
dexp <- dgood[ids_rep, ]

# log the positives
dexp$par_values_log <- dexp$par_values
dexp$par_values_log[, 2:n_coef] <- log(dexp$par_values_log[, 2:n_coef])

f1 <- FitGMM(dexp$par_values_log)
r1 <- rGMM(n_expand, d = n_coef, k = 1, means = f1@Mean, cov = f1@Covariance)

# good. compare also with a skew-mvn
# library(sn)
f2 <- sn::selm.fit(
  x = matrix(rep(1, nrow(dexp)), ncol = 1),
  y = dexp$par_values_log,
  family = "SN"
)
f2$param
d2 <- sn::makeSECdistr(dp = f2$param$dp, family = "SN", compNames = par_names)
plot(d2)
r2 <- sn::rmsn(n = 3e4, dp = f2$param$dp)

par(mfrow = c(2, 3))
for(i in 1:n_coef) {
  plot(density(dexp$par_values_log[, i]), main = par_names[i])
  lines(density(r1[, i]), col = 2)
  lines(density(r2[, i]), col = 3)
}
par(mfrow = c(1, 1))


# unlog
r1ori <- r1
r1ori[, 2:n_coef] <- exp(r1[, 2:n_coef])

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
  q_aprox <- r3ori[order(r3ori[, i]), i]
  q_obs <- dexp$par_values[order(dexp$par_values[, i]), i]
  plot(q_obs ~ q_aprox, pch = 19, col = rgb(0, 0, 0, 0.03),
       main = par_names[i])
  abline(0, 1)
}
par(mfrow = c(1, 1))

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


# Esto es un muy buen ajuste para el 95 % de los datos!!
# Y si usamos este método, podemos muestrear la posterior del modelo mixto con
# Stan!!!