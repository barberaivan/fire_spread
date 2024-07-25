# Failed to estimate separate likelihoods for every fire. By the moment, try a 
# fixed effects model. 

# I dont know how to perform the simulations. Is there RAM enough to bear all fires
# in memory? Does it take too long to load and remove all fires in each wave?


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
library(microbenchmark)

source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for rast_from_mat and a few constants

# source("estimation_functions.R") # prior_dist and other stuff

# Multicore settings -----------------------------------------------------

n_cores <- 15
registerDoMC(n_cores)

# Data and constants -----------------------------------------------------

data_dir <- file.path("data", "focal fires data", "landscapes_ig-known")
filenames <- list.files(data_dir)

fire_ids <- sapply(filenames, function(x) {
  # x <- "fire_data_raw_CoVentana.tif"
  id <- strsplit(x, "[.]")[[1]][1]
  return(id)
}) %>% unname
n_fires <- length(fire_ids)

# dir to save output
target_dir <- file.path("files", "pseudolikelihood_estimation")

# constants for fire spread simulation
upper_limit <- 0.5
par_names <- c("intercept", "vfi", "tfi", "slope", "wind", "fwi")
n_coef <- length(par_names)

# number of particles to choose by wave
n_pw <- 1000 

# number of fires to simulate by particle
n_sim <- 20

# formula for in_bounds model (gam)
bounds_model_formula <- formula(
  cbind(in_count, out_count) ~ 
    s(intercept, k = 8, bs = "cr") + 
    s(vfi, k = 8, bs = "cr") + 
    s(tfi, k = 8, bs = "cr") +     
    s(slope, k = 8, bs = "cr") + 
    s(wind, k = 8, bs = "cr") + 
    s(fwi, k = 8, bs = "cr")
)

# Functions ---------------------------------------------------------------

# prior distribution to simulate parameters or to compute them from a sobol
# sequence (type = "quantile", which computes the icdf.)
prior_dist <- function(mu_int = 0, sd_int = 15,
                       r_slope = 0.04,
                       r_fi = 0.15,
                       r_wind = 0.3,
                       r_fwi = 0.15,
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
      "fwi" = rexp(n, r_fwi)
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
    
    if(length(dim(p)) < 2) p <- matrix(p, ncol = 6)
    if(length(dim(p)) == 2) {
      if(ncol(p) != 6) p <- matrix(as.numeric(p), ncol = 6)
    }
    
    q <- cbind(
      "intercept" = qnorm(p[, 1], mu_int, sd_int),
      "vfi" = qexp(p[, 2], r_fi),
      "tfi" = qexp(p[, 3], r_fi),
      "slope" = qexp(p[, 4], r_slope),
      "wind" = qexp(p[, 5], r_wind),
      "fwi" = qexp(p[, 6], r_fwi)
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
similarity_simulate_particle <- function(particle, n_sim = 20,
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
  
  # extract values
  ll_summ <- c(
    overlap = mean(metrics[, "overlap"]),
    overlap_var = var(metrics[, "overlap"]),       
    size_diff = mean(metrics[, "size_diff"]),
    edge =  mean(metrics[, "edge"])
  )
  
  return(ll_summ)
}

# Emulate the loglik over a list of particles in parallel. Returns the metrics
# in rows.
similarity_simulate_parallel <- function(particle_ids = 1:100,
                                         fire_data = NULL,
                                         wave = NULL,
                                         n_sim = 20) {
  
  # turn particle matrix into list for parallel evaluation
  particles_mat_local <- particles_all[particle_ids, , drop = F]
  
  # define the fire's intercept as intercept + b_fwi * ((fwi[f] - mean) / sd)
  fwi_z <- (fire_data[["fwi"]]$fwi_expquad_day - fwi_mean) / fwi_sd
  
  new_intercept <- particles_mat_local[, "intercept"] + 
                   particles_mat_local[, "fwi"] * fwi_z
  
  particles_mat2 <- particles_mat_local[, -6]
  particles_mat2[, "intercept"] <- new_intercept
  
  particles_list <- lapply(1:nrow(particles_mat2), 
                           function(x) particles_mat2[x, ])
  
  # simulate in parallel
  result <- foreach(pp = particles_list) %dopar% {
    similarity_simulate_particle(pp, fire_data = fire_data, n_sim = n_sim)
  }

  # inherit previous dimnames
  names_single <- names(result[[1]])
  
  # rbind list result
  res <- do.call("rbind", result) %>% as.data.frame()
  colnames(res) <- names_single
  rownames(res) <- particle_ids
  # write intermediate output to disk
  saveRDS(res, file.path("files", "pseudolikelihood_estimation",
                         paste("wave_", wave, "_", fire_data$fire_id, ".rds",
                               sep = "")))
  
  # tidy this later  
  # # add useful columns
  # res$particle_id <- particle_ids
  # res$wave <- NA
  # res$par_values <- particles_mat_local
  
  return(res)
}

# function to repeat the similarity simulate parallel, returning an array
# (3D: particles, metrics, fires) and writing the result to disk.
# this loads and removes every fire, and prints "simulating tal fire"
simulation_wave <- function(particle_ids, wave, n_sim = 20) {
  
  sims_list <- lapply(1:n_fires, function(f) {
    message(paste("simulating fire", fire_ids[f]))
    full_data <- readRDS(file.path(data_dir, filenames[f]))
    # subset data needed for spread (to be cloned across workers)
    spread_data <- full_data[c("landscape", "burnable", "ig_rowcol",
                               "burned_layer", "burned_ids",
                               "fwi", "fire_id")]
    
    res <- similarity_simulate_parallel(particle_ids, spread_data, wave,
                                        n_sim = n_sim)
    rm(full_data, spread_data)
    return(res)
  })
  
  # inherit previous dimnames
  dn_single <- dimnames(sims_list[[1]])
  
  # turn into array
  res <- abind(sims_list, along = length(dn_single) + 1)
  
  # set dimnames
  dimnames(res) = c(
    dn_single,
    list(fire_id = fire_ids)
  )
  
  # write raw output
  saveRDS(res, file.path("files", "pseudolikelihood_estimation",
                         paste("wave_", wave, "_all_fires.rds", sep = "")))
  
  # remove temporary results (by fire) from disk
  ff_all <- list.files(file.path("files", "pseudolikelihood_estimation"))
  ff_delete <- ff_all[grep("all_fires", ff_all, invert = TRUE)]
  lapply(ff_delete, function(f) {
    unlink(file.path("files", "pseudolikelihood_estimation", f))
  })
  
  gc()
  
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
    similarity_bounds["largest", "overlap_var"], # overlap variance, original scale
    similarity_bounds["smallest", "overlap_var"]
  ) %>% max
    
  lowest_var <- low_var * var_factor 
  
  edge_upper <- ifelse(
    is.null(edge_abs_upper),
    edge_prop_upper * similarity_bounds["largest", "edge"],
    edge_abs_upper
  )
  
  too_large <- ((data$size_diff >= size_diff_upper) & 
                (data$overlap_var <= lowest_var)) | 
               (data$edge > edge_upper)
  too_small <- (data$size_diff <= size_diff_lower) & 
               (data$overlap_var <= lowest_var)
  
  keep <- as.numeric(!(too_large | too_small))
  
  return(keep)  
} 

# function to summarize overlap and evaluate in_bounds in all fires and particles.
# It takes as input the simulations array (3D: particles, metrics, fires).
# it should return the mean overlap across fires and the number of fires
# in_bounds, and the total (for the in_model.)
# here it should give the "overlap_logit" and "par_values" column.

# It needs a list with the similarity bounds for every fire, which will be 
# constant
summarize_simulations <- function(sim_array,
                                  prop = 0.85, var_factor = 2,
                                  edge_prop_upper = 0.05,
                                  edge_abs_upper = NULL) {
  
  # the 3rd dimension is fire_id
  
  # compute matrix with in_bounds vector for every fire.
  in_mat <- do.call("cbind", lapply(1:n_fires, function(f) {
    inn <- in_bounds(sim_array[, , f] %>% as.data.frame,
                     similarity_bounds = similarity_bounds[[f]],
                     prop = prop, var_factor = var_factor,
                     edge_prop_upper = edge_prop_upper,
                     edge_abs_upper = edge_abs_upper)
  }))
  
  res1 <- data.frame(particle_id = as.numeric(dimnames(sim_array)[[1]]),
                     in_count = rowSums(in_mat),
                     out_count = n_fires - rowSums(in_mat),
                     in_prop = rowSums(in_mat) / n_fires)
  
  mmm <- apply(sim_array, 1:2, mean) %>% as.data.frame()
  
  res2 <- cbind(res1, mmm)
  res2$overlap_logit <- qlogis(mmm$overlap)
  res2$overlap_log <- log(mmm$overlap)
  res2$par_values <- particles_all[res1$particle_id, ]
  
  return(res2)
}


# fit the in_bounds model, specifying the data and a formula
fit_in_bounds <- function(
    particles_data, 
    form = bounds_model_formula
) {
  
  data <- cbind(in_count = particles_data$in_count, 
                out_count = particles_data$out_count,
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
  
  if(is.numeric(n_fit)) {
      if(nrow(data) < n_fit) {
      # warning("subset_particles: Small dataset, no need to subset.")
      return(data$particle_id)
    }
    
    n_best <- ceiling(n_fit * prop_best)
    n_others <- n_fit - n_best
    
    data <- data[order(data$overlap_logit, decreasing = TRUE), ]
    ids_best <- data$particle_id[1:n_best]
    ids_others <- sample(data$particle_id[(n_best+1) : nrow(data)], 
                         size = n_others, replace = F)
    
    return(c(ids_best, ids_others))
  }
  
  if(n_fit == "all") return(data$particle_id)
}

# regular spacing after the best ones
subset_particles2 <- function(data, n_fit = 800, prop_best = 0.3,
                              use_in_bounds = F) {
  
  if(use_in_bounds) data <- data[data$in_bounds == 1, ]
  
  if(is.numeric(n_fit)) {
    if(nrow(data) < n_fit) {
      # warning("subset_particles: Small dataset, no need to subset.")
      return(data$particle_id)
    }
    
    n_best <- ceiling(n_fit * prop_best)
    n_others <- n_fit - n_best
    
    data <- data[order(data$overlap_logit, decreasing = TRUE), ]
    ids_best <- data$particle_id[1:n_best]
    # ids_others <- sample(data$particle_id[(n_best+1) : nrow(data)], 
    #                      size = n_others, replace = F)
    ids_others <- floor(seq(n_best + 1, nrow(data), length.out = n_others))                      
    
    return(c(ids_best, ids_others))
  }
  
  if(n_fit == "all") return(data$particle_id)
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
#   overlap_logit: log of the mean of the simulated likelihoods (overlap) in each particle,
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
  
  # # ## testo
  # previous_wave = w3
  # loglik_tolerance = 1
  # loglik_lower = -2
  # n_pw = 800
  # n_fit = "all"
  # prop_best = 0.3
  # gp_optimize = TRUE
  # use_in_bounds = FALSE
  # n_sim = 20
  # # ## 
  
  # when a previous wave has been run
  if(!is.null(previous_wave)) {
    gp_list <- previous_wave$gps
    wave_last <- length(gp_list)    # number of previous waves
    gp_last <- gp_list[[wave_last]]
    message(paste("Wave", 1 + wave_last))
    
    ### agregar esto para que el threshold se calcule al comienzo
    loglik_threshold_new <- get_loglik_threshold(loglik_high,
                                                 loglik_tolerance,
                                                 loglik_lower)
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
  
  sims <- simulation_wave(particle_ids = new_particles,
                          wave = this_wave,
                          n_sim = n_sim)
  
  particles_data_new <- summarize_simulations(
    sim_array = sims
  ) 
  
  particles_data_new$wave <- this_wave

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
  
  # write result to disk
  saveRDS(result, file.path("files", "pseudolikelihood_estimation",
                            paste("wave_", this_wave, "_all_fires_gp-fitted.rds",
                                  sep = "")))
  
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
    
    pred <- loglik_model$pred(new_data[, par_names], se.fit = TRUE)
    
    new_data$mle <- pred$mean
    new_data$upper <- pred$mean + qnorm(0.975) * pred$se
    new_data$lower <- pred$mean - qnorm(0.975) * pred$se
    
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
                        varying = "intercept", # parameter to vary in the 1d plot 
                        color_point = c("wave_plot", "in_prop"),
                        latest_wave = FALSE,
                        response = "overlap_logit") { 
  
  # ## test
  # fitting_wave <- w2
  # varying = "all"
  # color_point = "in_prop"
  # latest_wave = F
  
  discrete_color <- ifelse(color_point == "in_prop", F, T)
  
  loglik_model <- fitting_wave$gps[[length(fitting_wave$gps)]]
  data <- fitting_wave$particles_data

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
      scale_color_viridis(end = 0.7, discrete = discrete_color) +      
      
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
      scale_color_viridis(end = 0.7, discrete = discrete_color) +      
      
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

# function to convert the result from gp_fit to a loglik_update result 
# to use loglik_plot
gp2plot <- function(gp_result) {
  w <- list(
    gps = list(gp_result$loglik_model),
    particles_data = gp_result$data,
    loglik_optim = gp_result$op
  )
  return(w)
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

# Define more constants based on functions (result saved) -----------------

# fires_list <- lapply(filenames, function(fname) {
#   fd <- readRDS(file.path(data_dir, fname))
#   return(fd)
# })
# names(fires_list) <- fire_ids
# # get fwi mean and sd to standardize
# fwi_values <- sapply(fires_list, function(x) {
#   x$fwi$fwi_expquad_day
# })
# names(fwi_values) <- fire_ids
# saveRDS(fwi_values, file.path("data", "focal fires data", "fwi_values_ig-known.rds"))
# # compute similarity bounds
# similarity_bounds <- lapply(fires_list, function(x) {
#   get_bounds(x)
# })
# names(similarity_bounds) <- fire_ids
# saveRDS(similarity_bounds, file.path("data", "focal fires data", "similarity_bounds_ig-known.rds"))
# # free memory
# rm(fires_list); gc()

similarity_bounds <- readRDS(file.path("data", "focal fires data", "similarity_bounds_ig-known.rds"))
fwi_values <- readRDS(file.path("data", "focal fires data", "fwi_values_ig-known.rds"))
fwi_mean <- mean(fwi_values)
fwi_sd <- sd(fwi_values)


# n_sim definition --------------------------------------------------------

# # cholila
# full_data <- readRDS(file.path(data_dir, "2015_50.rds"))
# # subset data needed for spread (to be cloned across workers)
# spread_data <- full_data[c("landscape", "burnable", "ig_rowcol",
#                            "burned_layer", "burned_ids",
#                            "fwi", "fire_id")]
# 
# n_cores <- 16 # 12, 14?
# registerDoMC(n_cores)
# 
# particles_try <- 1:50
# mbm <- microbenchmark(
#   n10 = {
#     r10 <- similarity_simulate_parallel(particles_try, spread_data, wave = 1,
#                                         n_sim = 10)
#   },
#   n20 = {
#     r20 <- similarity_simulate_parallel(particles_try, spread_data, wave = 1,
#                                         n_sim = 20)
#   },
#   times = 1,
#   unit = "seconds"
# )
# gc()
# mbm
# # n10 101.7542 using 12 workers, ram reached 75 %
# # n20 216.6474
# # 100 / 500 # 0.2 s by fire
# # 216 / (50 * 20) # 0.21 s / fire, using 
# 
# # n10 101.5733 using 14 workers, ram reached ~70 %
# # n20 197.8763
# 
# # n10 104 using 16 workers, ram reached ~80 %
# # n20 203
# # 100 / 500 # 0.2 s by fire
# # 216 / (50 * 20) # 0.21 s / fire, using 
# 
# # con n_sim = 20
# 
# # cholila tarda 
# # 20 * 1000 * 0.2 = 4000 s / 3600 = 1.11 h
# nsim <- 20; npart <- 1000; stime = 0.2 * 1.5 / 60
# (mins <- 20 * 800 * stime)
# (hs <- mins / 60)
# 
# 
# # varianza?
# plot(r10[, "overlap_var"] ~ r20[, "overlap_var"]); abline(0, 1)
# plot(r10[, "overlap"] ~ r20[, "overlap"]); abline(0, 1)
# 
# # usando entre 12 y 16 cores los tiempos casi no cambiaron. usaremos
# # 15 cores
# # cholila debería tardar aprox 1:20 h en simular 800 partículas con n_sim = 20
# 
# 
# # vamos a probar correr 50 partículas pero con todos los incendios, a ver
# # si hay problemas.
# time50 <- microbenchmark(
#   sim = {s1 <- simulation_wave(1:50, 1, n_sim = 20)},
#   times = 1
# )
# # 984 s para simular 50 partículas, 20 reps
# # O sea, 984/50 = 19.68 s por partícula
# # si tiro 1000 serán 19680 s = 5.45 h
# # quizás mejor 800
# # 800 * 19.68 / 3600 # 4.37 h

# correré de a 800

# Running waves ---------------------------------------------------------

# shrink the loglik_tolerance while advancing waves.
# if the fitted GP is too wiggly and the tolerance is small, 
# it will be hard to get new particles, because the 
# plausible region will be small. If that happens and there is still
# need to improve the similarity, perhaps fitting a GP with all the data and 
# filtering only using that large, last gp could help.

w1 <- loglik_update_cumulative(
    previous_wave = NULL,
    loglik_tolerance = 3,
    loglik_lower = NULL,
    n_pw = 800,
    n_fit = "all",
    n_sim = n_sim
)
p1 <- loglik_plot(w1, varying = "all", color_point = "in_prop")

w2 <- loglik_update_cumulative(
    previous_wave = w1,
    loglik_tolerance = 3,
    loglik_lower = NULL,
    n_pw = 800,
    n_fit = "all",
    prop_best = 0.3,
    gp_optimize = TRUE,
    use_in_bounds = FALSE,
    n_sim = n_sim
)
p2 <- loglik_plot(w2, varying = "all", color_point = "in_prop")

w3 <- loglik_update_cumulative(
  previous_wave = w2,
  loglik_tolerance = 2,
  loglik_lower = NULL,
  n_pw = 800,
  n_fit = "all",
  prop_best = 0.3,
  gp_optimize = TRUE,
  use_in_bounds = FALSE,
  n_sim = n_sim
)
p3 <- loglik_plot(w3, varying = "all", color_point = "in_prop")

w4 <- loglik_update_cumulative(
  previous_wave = w3,
  loglik_tolerance = 1,
  loglik_lower = -2.5,
  n_pw = 800,
  n_fit = "all",
  prop_best = 0.3,
  gp_optimize = TRUE,
  use_in_bounds = FALSE,
  n_sim = n_sim
)
p4 <- loglik_plot(w4, varying = "all", color_point = "in_prop",
                  response = "overlap")

# w4 ya corrió, en escala logit.
# Ahora lo cargo, reajusto el GP en escala overlap, guardo los cambios en la 
# w4 (incluye cambiar el threshold) y luego de eso va la ola 5. 
# Para eso antes testear usando 10 partículas y n_sim = 2

# haciendo esto el gp se ajusto plano, no se por que.
# ajustarlo usando todos los datos a ver si mejora.
w4 <- readRDS("/home/ivan/Insync/Fire spread modelling/fire_spread/files/pseudolikelihood_estimation/wave_4_all_fires_gp-fitted.rds")

# rows_fit <- w4$particles_data$wave == max(w4$particles_data$wave)

gp4 <- gpkm(
  X = w4$particles_data$par_values, 
  Z = w4$particles_data$overlap, 
  parallel = F, useC = TRUE, nug.max = 1,
  kernel = "matern52"
)
# 900 s dice que tarda el culiau, con 2400 datos, 15 min

# define loglik_high associated with this new GP and the new loglik_threshold,
# to be used in the next wave
loglik_model <- gp4
fitted_ll <- loglik_model$pred(w4$particles_data$par_values, #[rows_fit, ]
                               se.fit = F)
loglik_high <- max(fitted_ll)
loglik_threshold_new <- get_loglik_threshold(loglik_high,
                                             loglik_tolerance = 0.05,
                                             loglik_lower = 0.12)

# find MLE in the fitted GP
# max fitted parameters value
id_fitted <- w4$particles_data$particle_id#[rows_fit]
id_max <- id_fitted[which.max(fitted_ll)]
id_filter <- which(w4$particles_data$particle_id == id_max)
params_max_fitted <- w4$particles_data$par_values[id_filter, ]

# bounds
lowers <- apply(w4$particles_data$par_values, 2, min)
uppers <- apply(w4$particles_data$par_values, 2, max) 

message("Optimizing pseudo-likelihood surface")
op <- optim(params_max_fitted, gp_predict, loglik_model = loglik_model,
            se.fit = FALSE,
            control = list(fnscale = -1, maxit = 1e5),
            method = "L-BFGS-B", lower = lowers, upper = uppers)  

# replace GP 
w4$gps[[4]] <- loglik_model
w4$loglik_thresholds[4] <- loglik_threshold_new
w4$loglik_optim <- op
saveRDS(w4, "/home/ivan/Insync/Fire spread modelling/fire_spread/files/pseudolikelihood_estimation/wave_4_all_fires_gp-fitted-overlap-scale.rds")
p4 <- loglik_plot(w4, varying = "all", color_point = "in_prop",
                  response = "overlap")
w4$loglik_thresholds
# esto de arriba quedo en R. si el plot es decente, correr la ola 5.
# quizas desde ahora haga falta correr el gp con todos los datos.
# pero es probable que haya algun error de codigo, y que el gp
# se este ajustando mal.
# verificar esos detalles.
# sera que el GP no convirgio?

# ahi anduvo perfecto!


w5 <- loglik_update_cumulative(
  previous_wave = w4,
  loglik_tolerance = 0.03,
  loglik_lower = 0.12,
  n_pw = 800,
  n_fit = "all",
  prop_best = 0.3,
  gp_optimize = TRUE,
  use_in_bounds = FALSE,
  n_sim = n_sim,
  response = "overlap"
)
p5 <- loglik_plot(w5, varying = "all", color_point = "in_prop",
                  response = "overlap")


w6 <- loglik_update_cumulative(
  previous_wave = w5,
  loglik_tolerance = 0.03,
  loglik_lower = 0.13,
  n_pw = 800,
  n_fit = "all",
  prop_best = 0.3,
  gp_optimize = TRUE,
  use_in_bounds = FALSE,
  n_sim = n_sim,
  response = "overlap"
)
p6 <- loglik_plot(w6, varying = "all", color_point = "in_prop",
                  response = "overlap")


w7 <- loglik_update_cumulative(
  previous_wave = w6,
  loglik_tolerance = 0.03,
  loglik_lower = 0.14,
  n_pw = 800,
  n_fit = "all",
  prop_best = 0.3,
  gp_optimize = TRUE,
  use_in_bounds = FALSE,
  n_sim = n_sim,
  response = "overlap"
)
p7 <- loglik_plot(w7, varying = "all", color_point = "in_prop",
                  response = "overlap")
# corriendo w7


w8 <- loglik_update_cumulative(
  previous_wave = w7,
  loglik_tolerance = 0.02,
  loglik_lower = 0.14,
  n_pw = 800,
  n_fit = "all",
  prop_best = 0.3,
  gp_optimize = TRUE,
  use_in_bounds = FALSE,
  n_sim = n_sim,
  response = "overlap"
)
p8 <- loglik_plot(w8, varying = "all", color_point = "in_prop",
                  response = "overlap")
# corriendo w8
# then fit the GP using all the particles, and that's the final
# likelihood emulator
# 
# de correr una ola más, setear antes el umbral, para que sea restrictivo.


# Fit full GP at log-sum scale? --------------------------------------------

data_path <- file.path("files", "pseudolikelihood_estimation")
ff <- list.files(data_path)
ff2 <- ff[grep("gp-fitted", ff, invert = T)]

tabs <- lapply(ff2, function(x) {
  fname <- file.path(data_path, x)
  readRDS(fname)
})
lapply(tabs, dim)

aa <- abind::abind(tabs, along = 1)[, "overlap", ]
dim(aa)
logsum <- apply(log(aa), 1, sum)

# plot(logsum)
range(logsum)
range(exp(logsum))
logsum_shift <- logsum - max(logsum) # shift in log to get numbers 


# get latest wave
w8 <- readRDS(file.path(data_path, "wave_8_all_fires_gp-fitted.rds"))


data8 <- w8$particles_data
data8$logsum <- logsum
data8$logsumexp <- exp(logsum)

plot(data8$logsum ~ data8$par_values[, "intercept"])
abline(h = -150, col = 2, lwd = 2, lty = 2)

plot(data8$logsumexp ~ data8$par_values[, "intercept"])

plot(data8$logsum ~ data8$par_values[, "vfi"])
plot(data8$logsum ~ data8$par_values[, "tfi"])
plot(data8$logsum ~ data8$par_values[, "slope"])
plot(data8$logsum ~ data8$par_values[, "wind"])
plot(data8$logsum ~ data8$par_values[, "fwi"])

# plot(data8$logsumexp ~ data8$par_values[, "intercept"])

# Fit gp at logsum scale, using the 800 best
data8 <- data8[order(data8$logsum, decreasing = TRUE), ]
data8 <- data8[order(data8$overlap, decreasing = TRUE), ]
# View(data8)

plot(data8$overlap ~ data8$logsum, col = rgb(0, 0, 0, 0.1),
     pch = 19)
plot(data8$overlap ~ data8$logsumexp, col = rgb(0, 0, 0, 0.1),
     pch = 19)

gp8 <- gp_fit(data8[1:800, ], response = "logsum")
logsum[order(logsum, decreasing = T)]
plot(density(logsum)); abline(v = max(logsum), lty = 2)



# comparar la distribución de overlap across fires 
# para la partícula que maximice el logsum y el overlap
best_logsum <- which.max(logsum)
best_om <- which.max(rowMeans(aa))

hist(aa[best_logsum, ], breaks = 20)
hist(aa[best_om, ], breaks = 20)
plot(density(aa[best_logsum, ], from = 0, to = 1))
lines(density(aa[best_om, ], from = 0, to = 1), col = 2)
abline(v = 0.1, lty = 2)
# mirando la mejor partícula, logsum parece tener menor media y 
# menor varianza.
plot(ecdf(aa[best_logsum, ]))
plot(ecdf(aa[best_om, ]))


summary(aa[best_logsum, ])
summary(aa[best_om, ])
plot(aa[best_om, ] ~ aa[best_logsum, ],
     ylim = c(0, 0.8), xlim = c(0, 0.8)); abline(0, 1)


# y si agarro las mejores 100?
ovmean <- rowMeans(aa)
ov_ord_mean <- aa[order(ovmean, decreasing = T), ]
ov_ord_log <- aa[order(logsum, decreasing = T), ]

best <- 200
ovs_mean <- ov_ord_mean[1:best, ] %>% as.vector()
ovs_log <- ov_ord_log[1:best, ] %>% as.vector()

plot(ovs_mean ~ ovs_log, ylim = c(0, 0.8), xlim = c(0, 0.8),
     pch = 19, col = rgb(0, 0, 0, 0.2))
abline(0, 1)

plot(density(ovs_log, from = 0, to = 1))
lines(density(ovs_mean, from = 0, to = 1), col = 2)
abline(v = 0.1, lty = 2)


# There is no clear reason to prefer the overlap_logsum over the
# overlap_mean. 
# The latter has a wider distribution, with mean and median a bit 
# higher.

# Proceed to sample the posterior for the fixed effects model, using 
# a GP fitted with lots of data (3000?) using the overlap_mean.


# Final GP (w8) -----------------------------------------------------------

w8 <- readRDS(file.path(data_path, "wave_8_all_fires_gp-fitted.rds"))

particles_fit <- subset_particles(w8$particles_data,
                                  n_fit = 3000, prop_best = 0.05)
rows_fit <- which(w8$particles_data$particle_id %in% particles_fit)
plot(w8$particles_data$overlap ~ w8$particles_data$par_values[, "intercept"])
points(w8$particles_data$overlap[rows_fit] ~ w8$particles_data$par_values[rows_fit, "intercept"],
       pch = 19, col = rgb(1, 0, 0, 1))

data_fit <- w8$particles_data[rows_fit, ]

gp8 <- gp_fit(data_fit) # expected 1735 s = 28 min
gp8_wave <- gp2plot(gp8)
p8 <- loglik_plot(gp8_wave, varying = "all", color_point = "in_prop",
                  response = "overlap")
# biutiful

saveRDS(gp8_wave,
        file.path(data_path, "wave_8_all_fires_gp-fitted-with-more-data.rds"))
# check how long it takes to predict new values as a function of the training 
# set size. If it's too slow, train with fewer points (2000)

# fist try to sample 1000 iter from the posterior.
