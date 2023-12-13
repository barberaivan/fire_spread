# Code to choose the similarity metric to approximate the likelihood, based on
# simulated fires. 
# We will use the fire 2015_53, "Alerces 2015", which is small enough to reduce 
# the computational burden. 

# The distribution (%) of vegetation type in this landscape is

# shrubland   subalpine   wet          dry
# 41.02752    23.69743    13.09468     22.18038

# and pixel counts are
# 103823      59968       33137        56129

# n_rep parameter vectors will be simulated from the prior (n_rep = 15?).
# For each one n_sim fire events are going to be simulated, called "observed" 
# hereafter (n_sim = 10), using 4 fixed ignition points. By fixing the ignition 
# points, the simulated fires will be comparable with all the n_rep * n_sim 
# observed ones.

# For every observed fire and similarity metric (n_met = 7?) a Gaussian Process 
# will be fitted in n_waves waves to emulate the likelihood function. 
# In this model we will exclude the FWI effect and use a fixed-intercept. 
# After the last wave, we will obtain the MLE using optim (n_rep * n_met MLEs).

# For each parameter separately we will compute the squared difference between
# the real value and the posterior mean. We will choose the metric that
# makes the smaller difference in most parameters, probably with a preference
# for the simplest one (spatial overlap) if differences are small. 

# Overall, n_rep * n_met estimations will be carried out. In the first 
# wave all instances will start with the same particles, so the simulated fires
# at particles_01 will not be saved (only the similarity metrics will be). 
# After the first wave, the following particles to be evaluated will depend on
# the first GP, so the particles used will probably differ across estimations. 

# From the second wave on, simulated fires by particle will be written to disk,
# and when the particle is required, the fires will be loaded instead of 
# simulated to compute the discrepancy metrics.

# TEST: is simulating 10 fires slower than writing and reading them?


# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)

library(Rcpp)
library(terra)
library(abind)         # manage arrays

library(randtoolbox)   # sobol sequences
library(DiceKriging)   # Fit Gaussian Processes
library(mgcv)          # fit spline to choose non-bounded particles
library(microbenchmark)

library(foreach)       # parallelization
library(doMC)          # doRNG had problems with cpp functions

sourceCpp("spread_functions.cpp")
sourceCpp("similarity_functions.cpp")


# Multicore settings -----------------------------------------------------

n_cores <- parallel::detectCores()
registerDoMC(n_cores)

# Data and constants -----------------------------------------------------

# landscape to run the simulations
land_full <- readRDS(file.path("data", "focal fires data",
                               "landscapes_ig-known_non-steppe", "2015_53.rds"))
land <- land_full$landscape

dnames <- dimnames(land)[[3]]

# terra-raster of the same landscape (with elevation) for plotting purposes.
# only the raster structure is needed, not its data.
land_raster <- rast(file.path("data", "focal fires data",
                              "wind ninja files", "2015_53.tif"))

# constants for fire spread simulation

distances <- rep(30, 8) # sides
distances[c(1, 3, 6, 8)] <- 30 * sqrt(2)

upper_limit <- 0.5

# ignition points
ig_rowcol <- matrix(c(round(ncol(land) * 0.4), round(nrow(land) * 0.3))) - 1

# number of particles to choose by wave
n_pw <- 1000 

# number of real parameters to simulate
n_rep <- 20

# number of fires to simulate by parameter
n_sim <- 10

# total number of real fires
n_obs <- n_rep * n_sim

# number of similarity metrics to compare
n_met <- 7

metric_names <- c("overlap_sp",
                  "sp_norm_5050", "sp_norm_7525",
                  "sp_expquad_5050", "sp_expquad_7525",
                  "sp_quad_5050", "sp_quad_7525")

par_names <- c("intercept", "subalpine", "wet", "dry", "fwi",
               "aspect", "wind", "elevation", "slope") 

par_names_sub <- par_names[-(which(par_names == "fwi"))]

# Functions ---------------------------------------------------------------

# prior simulator for graphical prior predictive checks
prior_sim <- function(mu_int = 0, sd_int = 20, sd_veg = 5,
                      r_01 = 0.05, r_z = 0.15) {
  
  betas <- c(
    "intercept" = rnorm(1, mu_int, sd_int),   # shrubland logit (reference class)
    "subalpine" = rnorm(1, 0, sd_veg),        # veg coefficients
    "wet" = rnorm(1, 0, sd_int),
    "dry" = rnorm(1, 0, sd_int),
    "fwi" = 0,                                # ZERO HERE
    "aspect" = rexp(1, r_01),                 # positive (northing)
    "wind" = rexp(1, r_01),                   # positive
    "elevation" = (-1) * rexp(1, r_z),        # negative
    "slope" = rexp(1, r_01)                   # positive
  )
  
  return(betas)
}

# Prior simulator from cumulative probabilities obtained from Sobol sequences 
# in [0, 1] ^ d
# p is a matrix of n_particles * d sobol samples of cumulative probability.
prior_q <- function(p, mu_int = 0, sd_int = 20, sd_veg = 5,
                    r_01 = 0.05, r_z = 0.15) {
  
  # add null column for fwi in p
  p <- cbind(p[, 1:4], rep(0, nrow(p)), p[, 5:8])
  
  q <- matrix(NA, nrow(p), ncol(p))
  colnames(q) <- c("intercept", "subalpine", "wet", "dry", "fwi",
                   "aspect", "wind", "elevation", "slope") 
  colnames(p) <- colnames(q)
  
  # fill matrix with parameters samples
  q[, 1] <- qnorm(p[, 1], mean = mu_int, sd = sd_int)
  
  for (i in 2:4) {
    q[, i] <- qnorm(p[, i], mean = 0, sd = sd_veg)
  }
  
  names_01 <- which(colnames(q) %in% c("aspect", "wind", "slope"))
  for(i in names_01) { # [0, 1] predictors
    q[, i] <- qexp(p[, i], rate = r_01)
  }
  
  # negative effect for elevation
  q[, "elevation"] <- (-1) * qexp(p[, "elevation"], rate = r_z)
  
  # fill null fwi parameter
  q[, "fwi"] <- 0
  
  return(q)
}

# function to make particles from a sobol sequence
particles_sim <- function(N) {
  prior_q(sobol(N, dim = 8, seed = 123, init = TRUE), # without FWI
          mu_int = -1.0, sd_int = 4, sd_veg = 4, # changed sd from 0.05
          r_z = 1.3, r_01 = 1.0) 
}
# this prior is widened relative to the one used to simulate fires


# Function to turn burned matrix into SpatRaster (for plotting)
rast_from_mat <- function(m, fill_raster) { # fill_raster is a SpatRaster from terra
  mt <- t(m)
  for(i in 1:nrow(m)) mt[, i] <- m[i, ]
  r <- fill_raster[[1]]
  values(r) <- as.numeric(mt)
  return(r)
}

# Function to simulate and plot a few fires to choose parameters that make sense
# in the focal landscape.
fire_prior_sim <- function(prior = NULL) {
  
  sizes <- numeric(3)
  par(mfrow = c(1, 3))
  
  for(i in 1:3) {
    fire_prior <- simulate_fire_cpp(
      landscape = land[, , 1:7],
      burnable = land[, , "burnable"],
      ignition_cells = ig_rowcol,
      coef = prior,
      wind_layer = which(dnames == "wind") - 1,
      elev_layer = which(dnames == "elev") - 1,
      distances = distances,
      upper_limit = upper_limit
    )
    # plot
    burnable_rast <- rast_from_mat(land[, , "burnable"], land_raster)
    burned_rast <- rast_from_mat(fire_prior, land_raster)
    values(burnable_rast)[values(burned_rast) == 1] <- 2
    plot(burnable_rast, col = c("black", "green", "red"),
         main = paste("intercept =", round(prior["intercept"], 3)))
    
    sizes[i] <- sum(fire_prior)
  }
  
  par(mfrow = c(1, 1))
  
  return(sizes)
}


# Simulate similarity for a given particle on the set of observed fires.
# Returns an array [n_obs, n_sim, n_met] with the similarity metrics
# by observed fire and replicate.
# Add size difference to rule out particles making extremely large fires
similarity_simulate_particle <- function(particle, n_sim = 10) {
  
  metrics <- array(NA, dim = c(n_obs, n_sim, n_met + 1),
                   dimnames = list(
                     fire_id = 1:n_obs,
                     simulation = 1:n_sim,
                     metric = c(metric_names, "size_diff")
                   ))
  
  for(i in 1:n_sim) {
    
    fire_sim <- simulate_fire_compare(
      landscape = land[, , 1:7],
      burnable = land[, , "burnable"],
      ignition_cells = ig_rowcol,
      coef = particle, 
      wind_layer = which(dnames == "wind") - 1,
      elev_layer = which(dnames == "elev") - 1,
      distances = distances,
      upper_limit = upper_limit
    )
    
    for(o in 1:n_obs) {
      comp <- compare_fires_try(fire_sim, fires_real[[o]])
      
      # subset metrics to evaluate
      m_keep <- which(names(comp) %in% metric_names) 
      metrics[o, i, 1:n_met] <- comp[m_keep]
                                # remove non-spatial raw metrics (2 to 5)
      # add size_diff
      metrics[o, i, n_met + 1] <- sum(fire_sim$counts_veg) - 
                                  sum(fires_real[[o]]$counts_veg)                   
    }
  }
  
  return(metrics)
}

# Emulate the loglik over a list of particles in parallel. Returns the loglik
# arrays binded in a 4D array, where the last dimension is the particle
similarity_simulate_parallel <- function(particle_ids = 1:100) {
  
  # turn particle matrix into list for parallel evaluation
  particles_mat_local <- particles_all[particle_ids, ]
  particles_list <- lapply(1:nrow(particles_mat_local), 
                           function(x) particles_mat_local[x, ])
  
  # simulate in parallel
  result <- foreach(pp = particles_list) %dopar% {
    similarity_simulate_particle(pp)
  }
  
  # inherit previous dimnames
  dn_single <- dimnames(result[[1]])
  
  # turn into array
  res <- abind(result, along = length(dn_single) + 1)
  
  # set dimnames
  dimnames(res) = c(
    dn_single,
    list(particle = particle_ids)
  )

  return(res)
}

# function to compute mean and variance of a desired metric by particle
summarize_by_particle <- function(metrics_matrix, name = "ll") {
  
  ### TEST
  # metrics_matrix = metrics_array[obs_id, , "size_diff", ]; name = "size_diff"
  ### TEST
  
  summ_fun <- function(x) {
    a <- c(mean(x), var(x))
    names(a) <- paste(name, c("mean", "var"), sep = "_")
    return(a)
  }
  
  summaries <- apply(
    metrics_matrix, 
    which(names(dimnames(metrics_matrix)) %in% c("particle")),
    FUN = summ_fun
  ) %>% t %>% as.data.frame()
  
  return(summaries)
}

# extract_focal_loglik: function to extract values from the large loglik array, 
# subsetting the parameters replicate, the "observed" fire and the desired
# metrics. It summarizes the results from n_sim simulated fires with mean and var,
# returning a df where each row is a particle, and column names are 
# ids: particle id,
# fire_id: id of the observed fire. in data_id the corresponding replicate can 
#   be obtained,
# ll_mean: mean of the metric across the n_sim simulations,
# ll_var: var of the metric across the n_sim simulations,
# metric: name of the metric used.

extract_focal_loglik <- function(metrics_array, metric = "overlap_sp",
                                 obs_id = 1, get_size = TRUE) {
  
  ### TEST
  # metrics_array = metrics_array; metric = metric; obs_id = obs_id; get_size = get_size
  ### TEST
  
  if(!(metric %in% dimnames(metrics_array)[["metric"]])) {
    stop("Metric does not exist.")
  }
  
  # subset array for focal fire and metric:
  ll <- qlogis(metrics_array[obs_id, , metric, ])  ## USING LOGIT-LIKELIHOOD
  
  # get mean and var by particle
  ll_summ <- summarize_by_particle(ll)
  
  if(get_size) {
    sizediff_summ <- summarize_by_particle(
      metrics_matrix = metrics_array[obs_id, , "size_diff", ], 
      name = "size_diff"
  ) } else sizediff_summ <- NULL
  
  # particles dimension
  particle_dim <- which(names(dimnames(ll)) == "particle")
  
  # metadata
  metadata <- data.frame(ids = dimnames(ll)[[particle_dim]] %>% as.numeric,
                         obs_id = obs_id,
                         metric = metric,
                         wave = NA)
  
  # parameters values 
  par_values <- as.data.frame(particles_all[metadata$ids, par_names_sub])
  
  # merge all
  df <- cbind(metadata, par_values, ll_summ, sizediff_summ)
  
  return(df)
}


# Function to filter particles based on the previously fitted gps
gp_filter <- function(gp_list = NULL, threshold = 5, 
                      loglik_max = NULL, particle_ids = NULL) {
  
  # compute predictions for all particles at all GPs
  preds <- lapply(gp_list, function(x) {
    predict(
      x, type = "UK",
      newdata = particles_all[particle_ids, par_names_sub, drop = F]
    )
  })
  
  # extract mean and sd
  pred_means <- do.call("cbind", lapply(1:length(preds), 
                                        function(x) preds[[x]]$mean))
  pred_sds <- do.call("cbind", lapply(1:length(preds), 
                                      function(x) preds[[x]]$sd))
  
  # evaluate plausibility in all previous GPs, using each ones' loglik_max
  keep_wide <- do.call("cbind", lapply(1:ncol(pred_means), function(x) {
    keep <- (pred_means[, x] + 3 * pred_sds[, x]) >= (loglik_max[x] - threshold)
    return(as.numeric(keep))
  }))
  
  keep_these <- which(rowSums(keep_wide) == length(preds))
  
  return(particle_ids[keep_these])
}

# builder of new particles (filtering with the previous GPs). It returns the 
# particle ids, from particles_all, that should be included in the next wave.
more_particles <- function(gp_list = NULL, 
                           loglik_max = NULL,
                           last_particle = NULL,
                           threshold = 5,
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
      
      new_ids <- c(new_ids, 
                   gp_filter(gp_list = gp_list, loglik_max = loglik_max, 
                             threshold = threshold, particle_ids = try_ids))
      
      added <- length(new_ids)
      last_particle <- try_ids[length(try_ids)]
    }
    
    new_particles <- new_ids[1:n_pw]
  } else {
    new_particles <- 1:n_pw
  }
  
  return(new_particles)
}


# Function to get more data passing given filters. 
#   first, it predicts which particles are good based on the filtering, based
#   on simulations according to the probability of the filters. Once n_pw 
#   particles are selected, fires are simulated. 
#   Once the data is available, it checks how many simulated particles have 
#   actually passed the filters, and continues to simulate more particles until
#   the good realized particles is over n_pw.
 
# filter: character vector indicating the filters to pass. if NULL, no filter 
# is applied. "in_bounds" and/or "high_loglik" may be specified. 
# for "in_bounds", the in_model is required, and the previous data is used
# to find the maximum probability found at the moment in previous particles.
# for "high_loglik", a gp model and a loglik threshold are required. 
 
# returns a data.frame with the data of all simulated particles, 
more_data <- function(filters = NULL, 
                      gp = NULL, 
                      loglik_threshold = NULL,
                      in_model = NULL,
                      in_prob_max = NULL, 
                      last_particle = NULL,
                      n_pw = 1000) {
  
  
  if(!is.null(filters)) {
    
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
      
      new_ids <- c(new_ids, 
                   gp_filter(gp_list = gp_list, loglik_max = loglik_max, 
                             threshold = threshold, particle_ids = try_ids))
      
      added <- length(new_ids)
      last_particle <- try_ids[length(try_ids)]
    }
    
    new_particles <- new_ids[1:n_pw]
  } else {
    new_particles <- 1:n_pw
  }
  
  return(new_particles)
}

# fit the in_bounds model, specifying the data and a formula (default is 
# s(intercept, k = 10))
fit_in_bounds <- function(particles_data, 
                          form = formula(in_bounds ~ s(intercept, k = 10))) {
  
  m <- gam(form, data = particles_data, family = binomial(link = "logit"), 
           method = "REML")
  
  return(m)
}

# function to simulate which particles are to be evaluated, based on the relative
# probability of being in range.
simulate_in_bounds <- function(in_model, particle_ids) {
  
  # get maximum fitted probability to relativize predictions.
  pmax <- predict(in_model, type = "response") %>% max
  
  # predict in_range_prob
  prob <- predict(in_model, as.data.frame(particles_all[particle_ids, ]),
                  type = "response") / pmax
  use_bin <- rbinom(length(prob), size = 1, prob = prob)
  
  return(particle_ids[use_bin == 1])
}

# function to simulate which particles are to be evaluated, based on the the 
# probability of the loglik to be above a certain threshold.

simulate_high_loglik <- function(gp, particle_ids, loglik_threshold) {
  
  # predict
  preds <- predict(gp, particles_all[particle_ids, ], se = TRUE)
  
  # get probability for x > loglik_threshold
  prob <- 1 - pnorm(loglik_threshold, preds$mean, preds$sd)
  
  # predict in_range_prob
  use_bin <- rbinom(length(prob), size = 1, prob = prob)
  
  return(particle_ids[use_bin == 1])
}

# Add the in_bounds filter to get new particles.
more_particles2 <- function(gp_list = NULL, 
                           loglik_max = NULL,
                           last_particle = NULL,
                           threshold = 5,
                           n_pw = 1000,
                           in_model) {
  
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
      
      in_range_ids <- simulate_in_bounds(in_model, try_ids)
      
      new_ids <- c(new_ids, 
                   gp_filter(gp_list = gp_list, loglik_max = loglik_max, 
                             threshold = threshold, particle_ids = in_range_ids))
      
      added <- length(new_ids)
      last_particle <- try_ids[length(try_ids)]
    }
    
    new_particles <- new_ids[1:n_pw]
  } else {
    new_particles <- 1:n_pw
  }
  
  return(new_particles)
}



# remove_bounded: Filter particles according to their proximity to the bounds.
# it returns the indexes in data that should remain after the particles near the
# bounds are removed.
remove_bounded <- function(data, obs_id, 
                           prop = 0.96, var_threshold = 0.15) {
  
  # obs_id = 1
  # prop = 0.96
  # var_threshold = 0.15
  # data = w1_1$particles_data
  size_diff_max <- similarity_bounds["size_diff", obs_id, "largest"]
  size_diff_min <- similarity_bounds["size_diff", obs_id, "smallest"]
  
  size_diff_lower <- prop * size_diff_min
  size_diff_upper <- prop * size_diff_max
  
  keep <- as.numeric(
    (
      (data$size_diff_mean >= size_diff_lower) &  # size in bounds
      (data$size_diff_mean <= size_diff_upper)
    ) & (
      data$ll_var >= var_threshold                # non-zero variance
    )
  )
  
  # returns binary vector to fit bernoulli model
  return(list(keep_bin = keep, keep_ids = data$ids[keep == 1]))
} 



# In each wave, keep track of the following objects:

# gp: Gaussian process fitted in the previous step.
# particles_data: data.frame with the following columns
#   ids: particle id,
#   wave: latest wave in which the particle was used to fit a GP,
#   mean: mean of the simulated likelihoods in each particle,
#   var: variance of the simulated likelihoods in each particle,
#   [parameter_values]: columns with the parameter values.

loglik_update <- function(previous_wave = NULL,
                          # list with gp and particles_data. The last is a df
                          # with ids, wave, mean, var and param values
                          threshold = 5,
                          n_pw = 1000, 
                          metric = "overlap_sp", 
                          obs_id = 1,
                          get_size = TRUE,
                          trend_intercept = 0) {
  
  # # TEST --
  # previous_wave = NULL
  # metric = "overlap_sp"; obs_id = 1; get_size = TRUE
  
  # previous_wave <- w1; threshold = 5
  # # TEST --
  
  
  
  # when a previous wave has been run
  if(!is.null(previous_wave)) {
    
    gp_list <- previous_wave$gps
    wave_last <- length(gp_list)    # number of previous waves
    gp_last <- gp_list[[wave_last]]
    loglik_max <- previous_wave$loglik_max 
    # max loglik according to each of the previous GPs (except the last) on 
    # the particles used to fit each one.
    
    particles_data <- previous_wave$particles_data
    
    # if there is a previous wave, override the following arguments:
    metric <- particles_data$metric[1]
    obs_id <- particles_data$obs_id[1]
    get_size <- length(grep("size_diff", names(particles_data))) > 0
    
    # OVERRIDE COEF_TREND!
    # coef_trend <- 
    
    # filter evaluated particles according to last gp
    print("Filtering old particles")
    ids_eval <- particles_data[particles_data$wave == wave_last, "ids"]
    ids_keep <- gp_filter(list(gp_last), particle_ids = ids_eval,
                          threshold = threshold, 
                          loglik_max = loglik_max[wave_last])
    
    if(length(ids_keep) > 0) {
      ids_keep_index <- which(particles_data$ids %in% ids_keep)
      particles_data$wave[ids_keep_index] <- wave_last + 1  
    }
    
    # get new particles meeting the same criterion
    print("Getting more particles")
    new_particles <- more_particles(gp_list = gp_list, 
                                    loglik_max = loglik_max,
                                    last_particle = max(particles_data$ids),
                                    threshold = threshold,
                                    n_pw = n_pw)
    
    # The new particles could be evaluated for previous simulation here, 
    # because a previous estimation (for another metric) might have been 
    # performed.
    
  } else {
    print("Getting more particles")
    new_particles <- more_particles(n_pw = n_pw)
  }
  
  # Simulate loglik on new particles
  print("Simulating fires")
  metrics_array <- similarity_simulate_parallel(particle_ids = new_particles) 
  
  # tidy array as df, previously summarizing over simulations
  particles_data_new <- extract_focal_loglik(
    metrics_array, metric = metric, obs_id = obs_id, get_size = get_size
  )
  
  # get wave order
  this_wave <- ifelse(is.null(previous_wave), 
                      1, 
                      wave_last + 1)
  
  particles_data_new$wave <- this_wave
  
  # merge old and new particles data
  if (is.null(previous_wave)) {
    particles_data_join <- particles_data_new
  } else {
    particles_data_join <- rbind(particles_data, particles_data_new)
  }
  
  # filter for current gp
  this <- particles_data_join$wave == max(particles_data_join$wave)
  
  # fit gp
  print("Fitting GP")
  gp_new <- km(design = particles_data_join[this, par_names_sub], 
               response = particles_data_join[this, "ll_mean"], 
               noise.var = particles_data_join[this, "ll_var"])
  ## ADD FIXED INTERCEPT LATER!
  
  # evaluate which was the max loglik at the training particles, to use as
  # filter in the next wave
  fitted <- predict(gp_new, type = "UK",
                    newdata = particles_data_join[this, par_names_sub])$mean
  loglik_max_new <- max(fitted)
  
  # put all together
  
  # (make NULL the absent objects)
  if(is.null(previous_wave)) {
    gp_list <- NULL
    loglik_max <- NULL
  }
  
  result <- list(gps = c(gp_list, gp_new),
                 loglik_max = c(loglik_max, loglik_max_new),
                 particles_data = particles_data_join)
  
  return(result)
}  


# this function includes a particle filter to use only particles in_bounds, 
# and is more astringent regarding the similarity threshold.

loglik_update_2 <- function(previous_wave = NULL,
                          # list with gp and particles_data. The last is a df
                          # with ids, wave, mean, var and param values
                          threshold = 5,
                          n_pw = 1000, 
                          metric = "overlap_sp", 
                          obs_id = 1,
                          get_size = TRUE,
                          trend_intercept = 0) {
  
  # # TEST --
  # previous_wave = NULL
  # metric = "overlap_sp"; obs_id = 1; get_size = TRUE
  
  # previous_wave <- w1; threshold = 5
  # # TEST --
  
  
  
  # when a previous wave has been run
  if(!is.null(previous_wave)) {
    
    gp_list <- previous_wave$gps
    wave_last <- length(gp_list)    # number of previous waves
    gp_last <- gp_list[[wave_last]]
    loglik_max <- previous_wave$loglik_max 
    # max loglik according to each of the previous GPs (except the last) on 
    # the particles used to fit each one.
    
    particles_data <- previous_wave$particles_data
    
    # if there is a previous wave, override the following arguments:
    metric <- particles_data$metric[1]
    obs_id <- particles_data$obs_id[1]
    get_size <- length(grep("size_diff", names(particles_data))) > 0
    
    # OVERRIDE COEF_TREND!
    # coef_trend <- 
    
    # filter evaluated particles according to last gp
    print("Filtering old particles")
    ids_eval <- particles_data[particles_data$wave == wave_last, "ids"]
    ids_keep <- gp_filter(list(gp_last), particle_ids = ids_eval,
                          threshold = threshold, 
                          loglik_max = loglik_max[wave_last])
    
    if(length(ids_keep) > 0) {
      ids_keep_index <- which(particles_data$ids %in% ids_keep)
      particles_data$wave[ids_keep_index] <- wave_last + 1  
    }
    
    # get new particles meeting the same criterion
    print("Getting more particles")
    new_particles <- more_particles(gp_list = gp_list, 
                                    loglik_max = loglik_max,
                                    last_particle = max(particles_data$ids),
                                    threshold = threshold,
                                    n_pw = n_pw)
    
    # The new particles could be evaluated for previous simulation here, 
    # because a previous estimation (for another metric) might have been 
    # performed.
    
  } else {
    print("Getting more particles")
    new_particles <- more_particles(n_pw = n_pw)
  }
  
  # Simulate loglik on new particles
  print("Simulating fires")
  metrics_array <- similarity_simulate_parallel(particle_ids = new_particles) 
  
  # tidy array as df, previously summarizing over simulations
  particles_data_new <- extract_focal_loglik(
    metrics_array, metric = metric, obs_id = obs_id, get_size = get_size
  )
  
  # get wave order
  this_wave <- ifelse(is.null(previous_wave), 
                      1, 
                      wave_last + 1)
  
  particles_data_new$wave <- this_wave
  
  # merge old and new particles data
  if (is.null(previous_wave)) {
    particles_data_join <- particles_data_new
  } else {
    particles_data_join <- rbind(particles_data, particles_data_new)
  }
  
  # filter for current gp
  this <- particles_data_join$wave == max(particles_data_join$wave)
  
  # fit gp
  print("Fitting GP")
  gp_new <- km(design = particles_data_join[this, par_names_sub], 
               response = particles_data_join[this, "ll_mean"], 
               noise.var = particles_data_join[this, "ll_var"])
  ## ADD FIXED INTERCEPT LATER!
  
  # evaluate which was the max loglik at the training particles, to use as
  # filter in the next wave
  fitted <- predict(gp_new, type = "UK",
                    newdata = particles_data_join[this, par_names_sub])$mean
  loglik_max_new <- max(fitted)
  
  # put all together
  
  # (make NULL the absent objects)
  if(is.null(previous_wave)) {
    gp_list <- NULL
    loglik_max <- NULL
  }
  
  result <- list(gps = c(gp_list, gp_new),
                 loglik_max = c(loglik_max, loglik_max_new),
                 particles_data = particles_data_join)
  
  return(result)
}  



# function to make new data varying only one predictor.
# the mle, if provided, must be named.
make_newdata <- function(varying = "intercept", data, mle = NULL) {
  
  values_list_mean <- lapply(par_names_sub, function(v) {
    if(v != varying) {
      res <- mean(data[, v])
    } else {
      res <- seq(min(data[, v]), max(data[, v]), length.out = 150)
    }
    return(res)
  })
  
  values_list_mle <- lapply(par_names_sub, function(v) {
    if(v != varying) {
      res <- mle[v]
    } else {
      res <- seq(min(data[, v]), max(data[, v]), length.out = 150)
    }
    return(res)
  })
  
  names(values_list_mean) <- names(values_list_mle) <- par_names_sub
  
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
partial_predictions <- function(varying = "intercept", data, gp, mle) {
  
  ### TEST
  # varying = "all"; data = w1$particles_data; gp = w1$gps[[length(w1$gps)]]
  # mle = mle
  ### 
  
  if(varying != "all") {
    new_data <- make_newdata(varying = varying, data = data, mle = mle)
    pred <- predict(gp, newdata = new_data[, par_names_sub], type = "UK")
    
    new_data$mle <- pred$mean
    new_data$upper <- pred$upper95
    new_data$lower <- pred$lower95
  } else {
    new_data <- do.call("rbind", lapply(par_names_sub, function(x) {
      make_newdata(varying = x, data = data, mle = mle)
    }))
    
    pred <- predict(gp, newdata = new_data[, par_names_sub], type = "UK")
    new_data$mle <- pred$mean
    new_data$upper <- pred$upper95
    new_data$lower <- pred$lower95
  }
  
  return(new_data)
}

# gp_plot: function to explore the advance of the gp.
gp_plot <- function(fitting_wave, 
                    varying = "intercept", # parameter to vary in the 1d plot 
                    fixed_at = "mean", # fix the not-varying parameters at mean or mle
                    local_range = FALSE,
                    filter_bounds = TRUE) { 
  
  #### TEST
  # fitting_wave <- w1; varying = "intercept"; fixed_at = "mean"
  # varying = "all"; local_range = FALSE
  #### 
  
  gp <- fitting_wave$gps[[length(fitting_wave$gps)]]
  data <- fitting_wave$particles_data
  data$in_bounds <- factor(as.character(data$in_bounds), levels = c("0", "1"))
  
  # get current wave and add to data
  current <- which(data$wave == max(data$wave))
  old <- which(data$wave < max(data$wave))
  
  data$wave_plot <- "previous"
  data$wave_plot[current] <- "current"
  
  # subset data for partial predictions at local range
  if(local_range) {
    data_pred <- data[current, ]
  } else {
    data_pred <- data
  }
  
  
  # if(fixed_at == "mean") mle <- NULL
  # else {
  #   ll <- function(betas) {
  #     nd <- matrix(betas, nrow = 1)
  #     colnames(nd) = names(betas)
  #     -predict(gp, newdata = nd, type = "UK")$mean
  #   }
  #   opt <- optim(colMeans(as.matrix(data_pred[, par_names_sub])), fn = ll)
  #   mle <- opt$par
  #   names(mle) <- par_names_sub
  # }
  
  # get mle for partial predictions
  
  ll <- function(betas) {
    nd <- matrix(betas, nrow = 1)
    colnames(nd) = names(betas)
    -predict(gp, newdata = nd, type = "UK")$mean
  }
  opt <- optim(colMeans(as.matrix(data_pred[, par_names_sub])), fn = ll)
  mle <- opt$par
  
  # compute partial predictions 
  preds <- partial_predictions(varying = varying, 
                               data = data_pred, 
                               gp = gp, mle = mle)
  
  color_point <- ifelse(filter_bounds, "in_bounds", "wave_plot")
  
  if(varying != "all") {
    
    p <-
      ggplot() +
      
      # data
      geom_point(data = data, 
                 mapping = aes_string(x = varying, y = "ll_mean",
                                      color = color_point, shape = color_point),
                 size = 2, alpha = 0.5) + 
      scale_shape_manual(values = c(16, 17)) +
      scale_color_viridis(end = 0.7, discrete = TRUE) +      
      
      # gp predictions 
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
      ylab("loglik") + 
      xlab("parameter value")
    
    print(p)
    return(list(plot = p, data_points = data, data_preds = preds))    
  } else {
    
    # when all parameters are varied, data must be replicated to be plotted 
    # against every parameter.
    
    data_expand <- do.call("rbind", lapply(par_names_sub, function(n) {
      data$varying_val <- data[, n]
      data$varying_var <- n
      return(data)
    }))
    
    data_expand$varying_var <- factor(data_expand$varying_var,
                                      levels = par_names_sub)
    
    preds$varying_var <- factor(preds$varying_var,
                                levels = par_names_sub)

    p <-
    ggplot() +

      # data
      geom_point(data = data_expand, 
                 mapping = aes_string(x = "varying_val", y = "ll_mean",
                                      color = color_point, shape = color_point),
                 size = 2, alpha = 0.5) + 
      scale_shape_manual(values = c(16, 17)) +
      scale_color_viridis(end = 0.7, discrete = TRUE) +      
      
      # gp predictions 
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
      ylab("loglik") + 
      xlab("parameter value")
    
    print(p)
    return(list(plot = p, data_points = data_expand, data_preds = preds))
  }
}


# Parameters simulation --------------------------------------------------

# Choose parameter vectors by hand, visually inspecting a few simulations
# from every one to check that fire sizes are OK with this landscape.
# (Finding a prior that would simulate proper-sized fires was too hard.)

# pp <- prior_sim(mu_int = -1.0, sd_int = 0.02, sd_veg = 0.02,
#                 r_z = 1.5, r_01 = 1.1) 
# pp
# fire_prior_sim(pp)
# 
# # param_mat <- matrix(NA, n_rep, 9)
# param_mat[30, ] <- pp
# saveRDS(param_mat, "files/param_mat.rds")
param_mat <- readRDS("files/param_mat.rds")

# Simulate real fires ------------------------------------------------------

# n_rep (parameters) * n_sim (fires) will be simulated in the same landscape. 
seeds <- matrix(100:(100 + n_rep * n_sim - 1),
                nrow = n_sim, ncol = n_rep)

data_id <- expand.grid(sim = 1:n_sim, 
                       rep = 1:n_rep)

fires_real <- vector(mode = "list", length = n_obs)

for(i in 1:n_obs) {
  # i = 1  
  r = data_id$rep[i]
  s = data_id$sim[i]
    
  set.seed(seeds[s, r])
    
  fires_real[[i]] <- simulate_fire_compare(
    landscape = land[, , 1:7],
    burnable = land[, , "burnable"],
    ignition_cells = ig_rowcol,
    coef = param_mat[r, ],
    wind_layer = which(dnames == "wind") - 1,
    elev_layer = which(dnames == "elev") - 1,
    distances = distances,
    upper_limit = upper_limit
  )
}


# Get bounds on similarity metrics ----------------------------------------

# Given the landscape resolution and size, similarity metrics between observed
# and simulated fires are bounded from below, hindering a proper approximation 
# of the unknown likelihood. This causes the similarity to be overestimated
# when simulated fire sizes are more extreme than what the landscape resolution
# and size allow to notice. For example, if a given particle would simulate a 
# fire that burns all the landscape and more, the similarity metric would
# correspond to a simulated fire that burned all the landscape exactly, not more.

# A simple way to overcome this problem is to reject particles close to the 
# lower bounds of the similarity metrics. In most cases, the lower bound for 
# large simulated fires is higher than the lower bound for small fires
# (considering overlap_sp: large_fires ~ 0.05; small_fires ~ 1 / 111 = 0.009).

# To remove regions of the parameter space below the highest of the lower bounds
# (similarity_min), we define the similarity_threshold for every metric as
#   similarity_min = max(lower_bound_large, lower_bound_small)
#   similarity_threshold = similarity_min + delta
# Adding delta is necessary because we cannot know whether particles near 
# similarity_min would have been lower or not if the landscape would have 
# allowed it. Another option is
#   similarity_threshold = similarity_min + p * (1 - similarity_min) 
# with p ~ 0.25
p <- 0.25
similarity_min <- 0.07
(similarity_threshold = similarity_min + p * (1 - similarity_min))

# With this threshold we can emulate the likelihood using a GP, and filter out
# particles < similarity_threshold. After this first astringent wave, we can 
# proceed with the method proposed by Wilkinson (2014).

# Ideally, the first wave should be estimated with many particles, so mistakes 
# are not performed there. Then, GPs could use data in a cumulative manner, so 
# the last one is a good predictor of the likelihood above similarity_threshold.
# (This allows to compute the likelihood instead of just rejecting the particle).

# simulate largest and smallest fires
fire_largest <- simulate_fire_compare(
  landscape = land[, , 1:7],
  burnable = land[, , "burnable"],
  ignition_cells = ig_rowcol,
  coef = c(1e6, rep(0, 8)),
  wind_layer = which(dnames == "wind") - 1,
  elev_layer = which(dnames == "elev") - 1,
  distances = distances,
  upper_limit = upper_limit
)

fire_smallest <- simulate_fire_compare(
  landscape = land[, , 1:7],
  burnable = land[, , "burnable"],
  ignition_cells = ig_rowcol,
  coef = c(-1e6, rep(0, 8)),
  wind_layer = which(dnames == "wind") - 1,
  elev_layer = which(dnames == "elev") - 1,
  distances = distances,
  upper_limit = upper_limit
)

similarity_bounds <- array(NA, dim = c(n_met + 1, n_obs, 2),
                          dimnames = list(
                            metric = c(metric_names, "size_diff"),
                            fire_id = 1:n_obs,
                            extreme = c("smallest", "largest")
                          ))
for(i in 1:n_obs) {
  print(i)# = 1
  similarity_bounds[1:n_met, i, 1] <- compare_fires_try(
    fire_smallest, fires_real[[i]]
  )[metric_names]
  
  similarity_bounds[1:n_met, i, 2] <- compare_fires_try(
    fire_largest, fires_real[[i]]
  )[metric_names]
  
  # add size difference
  similarity_bounds[n_met + 1, i, 1] <- sum(fire_smallest$counts_veg) - 
                                        sum(fires_real[[i]]$counts_veg)
  
  similarity_bounds[n_met + 1, i, 2] <- sum(fire_largest$counts_veg) - 
                                        sum(fires_real[[i]]$counts_veg)
}

simil_lower_large <- similarity_simulate_particle(c(1e6, rep(0, 8)), 
                                                  n_sim = 1)[, 1, ]
simil_lower_small <- similarity_simulate_particle(c(-1e6, rep(0, 8)), 
                                                  n_sim = 1)[, 1, ]

similarity_min <- pmax(simil_lower_large, simil_lower_small)

# compare them with observed ones.


# Long sobol sequence -----------------------------------------------------

# for all estimations to run over the same particles (when possible), a very long
# sobol sequence will be produced, expecting not to need further "simulations".
# A larger sequence will be created from scratch if more particles are needed,
# to ensure the same particles are always used.
n_large <- 300000
particles_all <- particles_sim(n_large)

# Exploring 1 wave --------------------------------------------------------

w1 <- loglik_update(n_pw = 1500)
p1 <- gp_plot(w1, varying = "all")
p1 +  ylim(-4, 4)

w1_1 <- w1
w1_1$particles_data$in_bounds <- remove_bounded(w1_1$particles_data)

# fit model to in_bounds using 1500 particles
in_model <- gam(in_bounds ~ s(intercept, k = 10), method = "REML", 
                data = w1_1$particles_data, family = binomial(link = "logit")) 
# in_model <- gam(in_bounds ~ intercept + I(intercept ^ 2), method = "REML", 
#                 data = w1_1$particles_data, family = "binomial") 

d <- data.frame(intercept = seq(-10, 10, length.out = 150))
d$prob <- predict(in_model, d, type = "response")
plot(prob ~ intercept, d, type = "l")

# If only a few points are judged to be in range, the average in_probability 
# could be low. Then, we relativize predictions to the maximum fitted 
# probability, so the not-in-range values are not simulated so often:
in_prob_max <- max(predict(in_model, type = "response")) # DO NOT USE FITTED()!

d$prob_scaled <- d$prob / in_prob_max
plot(prob_scaled ~ intercept, d, type = "l")

# simulate more in_range particles:
try_these <- max(w1$particles_data$ids) : (max(w1$particles_data$ids) * 3)
ss <- simulate_in_bounds(in_model, try_these)
hist(particles_all[ss, "intercept"])
hist(particles_all[try_these, "intercept"], add = TRUE)
plot(density(particles_all[ss, "intercept"]))
lines(density(particles_all[try_these, "intercept"]), col = 2)
# Bien! Este filtro funciona hermoso



# loglik_update por partes ------------------------------------------------

previous_wave <- w1; threshold = 3 # menor!!!


gp_list <- previous_wave$gps
wave_last <- length(gp_list)    # number of previous waves
gp_last <- gp_list[[wave_last]]
loglik_max <- previous_wave$loglik_max 
# max loglik according to each of the previous GPs (except the last) on 
# the particles used to fit each one.

particles_data <- previous_wave$particles_data

# if there is a previous wave, override the following arguments:
metric <- particles_data$metric[1]
obs_id <- particles_data$obs_id[1]
get_size <- length(grep("size_diff", names(particles_data))) > 0

# OVERRIDE COEF_TREND!
# coef_trend <- 

# filter evaluated particles according to last gp
print("Filtering old particles")

##### NEW : remove the out of bounds, but keep the others
ids_eval <- particles_data[particles_data$wave == wave_last, "ids"]

# first, remove those near the bounds
particles_data$in_bounds <- remove_bounded(particles_data, 
                                           obs_id = obs_id)$keep_bin

# ids_keep <- gp_filter(list(gp_last), particle_ids = ids_in_bounds,
#                       threshold = threshold, 
#                       loglik_max = loglik_max[wave_last])

# keep track of the wave irrespective of bounded or not.

# particles_data$wave[particles_data$in_bounds == 1] <- wave_last + 1  


# get new particles meeting the same criterion
print("Getting more particles")
new_particles <- more_particles2(gp_list = gp_list, 
                                loglik_max = loglik_max,
                                last_particle = max(particles_data$ids),
                                threshold = threshold,
                                n_pw = 3000,
                                in_model = in_model) ### NEW


# Simulate loglik on new particles
print("Simulating fires")
metrics_array <- similarity_simulate_parallel(particle_ids = new_particles) 

# tidy array as df, previously summarizing over simulations
particles_data_new <- extract_focal_loglik(
  metrics_array, metric = metric, obs_id = obs_id, get_size = get_size
)

## remove particles near the bounds
particles_data_new$in_bounds <- remove_bounded(particles_data_new, 
                                               obs_id = obs_id)$keep_bin
sum(particles_data_new$in_bounds) # cool!
## Acá podría hacer un while para simular más fuegos si los nuevos son muy pocos
# mantener los out_of_bounds que fueron simulados ayuda a actualizar el
# in_model

# get wave order
this_wave <- ifelse(is.null(previous_wave), 
                    1, 
                    wave_last + 1)

particles_data_new$wave <- this_wave

# merge old and new particles data
if (is.null(previous_wave)) {
  particles_data_join <- particles_data_new
} else {
  particles_data_join <- rbind(particles_data, particles_data_new)
}


# TRY a harder constraint on size limits:
particles_data_join$in_bounds_shrink <- remove_bounded(
  particles_data_join, obs_id, 
  prop = 0.75, var_threshold = 0.5
  )$keep_bin


# filter for current gp (USE ALL IN_BOUNDS PARTICLES)
# this <- particles_data_join$wave == max(particles_data_join$wave)
this <- particles_data_join$in_bounds == 1
this2 <- particles_data_join$in_bounds_shrink == 1

# fit gp
print("Fitting GP")
gp_new <- km(design = particles_data_join[this, par_names_sub], 
             response = particles_data_join[this, "ll_mean"], 
             noise.var = particles_data_join[this, "ll_var"])

gp_new2 <- km(design = particles_data_join[this2, par_names_sub], 
             response = particles_data_join[this2, "ll_mean"], 
             noise.var = particles_data_join[this2, "ll_var"])

gp_new3 <- km(design = particles_data_join[this2, par_names_sub], 
              response = particles_data_join[this2, "ll_mean"], 
              noise.var = particles_data_join[this2, "ll_var"],
              # formula = ~ intercept + I(intercept ^ 2) + 
              #   subalpine + I(subalpine ^ 2) +
              #   wet + I(wet ^ 2) +
              #   dry + I(dry ^ 2) +
              #   aspect + I(aspect ^ 2) +
              #   wind + I(wind ^ 2) +
              #   elevation + I(elevation ^ 2) +
              #   aspect + I(aspect ^ 2)                 ### esto no anduvo
              
              formula = ~ I((intercept - (-1)) ^ 2) #+
                # I(subalpine ^ 2) +
                # I(wet ^ 2) +
                # I(dry ^ 2) +
                # I(aspect ^ 2) +
                # I(wind ^ 2) +
                # I(elevation ^ 2) +
                # I(aspect ^ 2)                 ### ANDA, pero sigue estimando un intercept
              
              )

gp_new4 <- km(design = particles_data_join[this2, par_names_sub], 
              response = particles_data_join[this2, "ll_mean"], 
              noise.var = particles_data_join[this2, "ll_var"],
              formula = ~ 1, coef.trend = -10)

## ADD FIXED INTERCEPT LATER!

# evaluate which was the max loglik at the training particles, to use as
# filter in the next wave
fitted <- predict(gp_new, type = "UK",
                  newdata = particles_data_join[this, par_names_sub])$mean
loglik_max_new <- max(fitted)

# put all together

# (make NULL the absent objects)
if(is.null(previous_wave)) {
  gp_list <- NULL
  loglik_max <- NULL
}

result <- list(gps = c(gp_list, gp_new),
               loglik_max = c(loglik_max, loglik_max_new),
               particles_data = particles_data_join)
w2 <- result

w2_shrink <- list(gps = c(gp_list, gp_new2),
               loglik_max = c(loglik_max, loglik_max_new),
               particles_data = particles_data_join)

w2_quad <- list(gps = c(gp_list, gp_new3),
                  loglik_max = c(loglik_max, loglik_max_new),
                  particles_data = particles_data_join)

w2_int <- list(gps = c(gp_list, gp_new4),
                loglik_max = c(loglik_max, loglik_max_new),
                particles_data = particles_data_join)

p2 <- gp_plot(w2, varying = "intercept")
p2_shrink <- gp_plot(w2_shrink, varying = "intercept")
p2_quad <- gp_plot(w2_quad, varying = "intercept")
p2_int <- gp_plot(w2_int, varying = "intercept")

p2$plot + geom_hline(yintercept = qlogis(0.2), linetype = "dashed")
p2$plot + geom_hline(yintercept = qlogis(0.07), linetype = "dashed")

# el cuadrático así peladito suaviza demasiado. restringirle para que el optimo
# se quede en el MLE


# No están andando bien, a pesar de shrinkear más la definición del bound, 
# pareciera que siguen quedando partículas en el bound a la hora de hacer el ajuste.

# AUNQUE PUEDE SER ALGO GRÁFICO, QUE EN REALIDAD HAYA QUE PLOTEAR CON LA MEDIA DE LO
# USADO PARA EL MODELO!!! ALTERAR LA GP_PLOT PARA QUE HAGA ESO, QUE USE IN_BOUNDS.
# PORQUE EL MLE NO ESTÁ TAAAN MAL, SOLO QUE TIENE MUCHA INCERTIDUMBRE.

# TAMBIÉN ESTOY TENIENDO MUY POCOS DATOS PARA ESTIMAR EL GP


# OTRO APPROACH SERÍA PRIMERO HACER UN CLASIFICADOR DE LA P DE ESTAR ABAJO DEL
# UMBRAL MÍNIMO PARA LA SIMILARITY USADA, Y QUE SOLO SI ESTÁ POR ARRIBA LA 
# PARTÍCULA SEA USADA.

# AAAAH Y LA MEDIA SUBE HACIA INTERCEPT ALTOS PORQUE NO HAY INTERCEPTS BAJOS
# PARA QUE LOS VEA Y LOS AJUSTE. INCLUSO SI FILTRARA AFUERA A TODOS, NO PODRÍA
# ESTIMARLO BIEN.

# cosas de antes ----------------------------------------------------------



good_ones <- which(w1$particles_data$ll_mean > -2)#log(0.2))
length(good_ones)

w1_good <- w1
w1_good$particles_data <- w1$particles_data[good_ones, ]

w1_good$gps[[1]] <- km(
  design = w1_good$particles_data[, par_names_sub], 
  response = w1_good$particles_data[, "ll_mean"], 
  noise.var = w1_good$particles_data[, "ll_var"],
  formula = ~ intercept + I(intercept ^ 2) +
            subalpine + I(subalpine ^ 2) +
            wet + I(wet ^ 2) +
            dry + I(dry ^ 2) +
            aspect + I(aspect ^ 2) +
            wind + I(wind ^ 2) +
            elevation + I(elevation ^ 2) +
            aspect + I(aspect ^ 2)
)
p1_good <- gp_plot(w1_good, varying = "all")
# esto pasa porque el mle se nos va de rango. habría que resgringirlo para 
# que fuera local.


w2 <- loglik_update(previous_wave = w1, n_pw = 1000)
p2 <- gp_plot(w2, varying = "all")

good_ones <- which(w2$particles_data$ll_mean > -1.6)#log(0.2))
length(good_ones)

w2_good <- w2
w2_good$particles_data <- w2$particles_data[good_ones, ]

w2_good$gps[[1]] <- km(
  design = w2_good$particles_data[, par_names_sub], 
  response = w2_good$particles_data[, "ll_mean"], 
  noise.var = w2_good$particles_data[, "ll_var"],
  # formula = ~ intercept + I(intercept ^ 2) +
  #   subalpine + I(subalpine ^ 2) +
  #   wet + I(wet ^ 2) +
  #   dry + I(dry ^ 2) +
  #   aspect + I(aspect ^ 2) +
  #   wind + I(wind ^ 2) +
  #   elevation + I(elevation ^ 2) +
  #   aspect + I(aspect ^ 2)
)
p2_good <- gp_plot(w2_good, varying = "all")


w3 <- loglik_update(previous_wave = w2, n_pw = 1000)
p3 <- gp_plot(w3, varying = "all")

# alterar GP de la wave 3 para ver qué pasa con la unconditional mean.

# filtrar partículas
d3 <- w3$particles_data

# size_diff distribution and limits
size_obs <- sum(fires_real[[1]]$counts_veg)
size_diff_min <- 1 - size_obs
size_diff_max <- size_max - size_obs
size_diff_lower <- 0.95 * size_diff_min
size_diff_upper <- 0.95 * size_diff_max
ll_var_low <- 0.5

filtrador <- function(smean, svar) {
  filt <- ((smean < size_diff_lower) | (smean > size_diff_upper)) & 
          (svar < ll_var_low)
}

d3$inside <- as.numeric(!filtrador(d3$size_diff_mean, d3$ll_var))
d3$wave[(d3$inside == 1) & (d3$wave == 3)] <- 4
table(d3$wave)

use <- sample(which(d3$wave == 4), size = 500, replace = F)

# add fake data where I want
more <- 100
dextra <- d3[1:more, ]
dextra$ll_mean <- log(min_overlap)
dextra$intercept <- runif(more, 5, 10)

dadd <- rbind(d3[use, ], dextra)

gp4 <- km(design = dadd[, par_names_sub], 
          response = dadd[, "ll_mean"], 
          noise.var = dadd[, "ll_var"])#,
          #formula = ~ 1, coef.trend = log(min_overlap) * 2)

w4 <- list(gps = c(w3$gps, gp4),
           loglik_max = c(w3$loglik_max, 0),
           particles_data = dadd)

p4 <- gp_plot(w4, varying = "all", local_range = F)
p4 + ylim(-10, 5)

# forzar una unconditional mean cualquiera corrompe el ajuste de los GP.
# quizás haya que asumir una función determinada a ojo cuando el GP se sale 
# de rango y chau.


# explore simulation variance as a function of intercept
d3 <- w3$particles_data

ggplot(d3, aes(x = intercept, y = ll_var, color = size_diff_mean)) +
  geom_point() +
  geom_hline(yintercept = 0.75)

ggplot(d3[d3$size_diff_mean < 0 , ], 
       aes(x = size_diff_mean, y = ll_var)) +
  geom_point() + 
  geom_smooth()


# size_diff distribution and limits
size_obs <- sum(fires_real[[1]]$counts_veg)
size_diff_min <- 1 - size_obs
size_diff_max <- size_max - size_obs
size_diff_lower <- 0.98 * size_diff_min
size_diff_upper <- 0.98 * size_diff_max

par(mfrow = c(2, 1))
hist(d3$size_diff_mean[d3$size_diff_mean > 0], breaks = 100)
abline(v = c(size_diff_upper), col = 2)

hist(d3$size_diff_mean[d3$size_diff_mean < 0], breaks = 100)
abline(v = c(size_diff_lower), col = 2)
par(mfrow = c(1, 1))

# size limits to impose better constraints


range(d3$size_diff_mean)


# use loglik variance, not size variance. it has a better scale
# threshold for variance: 
# pero en realidad podemos saber qué es mucho y qué es poco en la varianza
# de log-metricas que van en [0, 1]
bbb <- 0.1
curve(dbeta(x, bbb, bbb), ylim = c(0, 10)); abline(h = 0)
var(log(rbeta(1e5, bbb, bbb))) # 75 is a lot

bbb <- 0.5
curve(dbeta(x, bbb, bbb), ylim = c(0, 10)); abline(h = 0)
var(log(rbeta(1e5, bbb, bbb))) # 3.32 tambien is a lot

bbb <- 1
curve(dbeta(x, bbb, bbb), ylim = c(0, 10)); abline(h = 0)
var(log(rbeta(1e5, bbb, bbb))) # 1 es la var de una uniforme!! o sea, un montón

bbb <- 10
curve(dbeta(x, bbb, bbb), ylim = c(0, 10)); abline(h = 0)
var(log(rbeta(1e5, bbb, bbb))) # 1 es la var de una uniforme!! o sea, un montón


filtrador <- function(smean, svar) {
  size_diff_lower <- 0.98 * size_diff_min
  size_diff_upper <- 0.98 * size_diff_max
  ll_var_low <- 0.1# quantile(d3$ll_var, probs = c(0.1))
  # esto podría ser algo como 0.05 * quantile(ll_var, 0.95),
  # para que sea invariable a la métrica.

  
  filt <- ((smean < size_diff_lower) | (smean > size_diff_upper)) & 
          (svar < ll_var_low)
}

bad <- filtrador(d3$size_diff_mean, d3$ll_var)
d3$keep_or_not <- "good"
d3$keep_or_not[bad] <- "bad"

dl <- pivot_longer(d3, cols = all_of(which(names(d3) %in% par_names_sub)),
                   names_to = "parameter", values_to = "par_value")
dl$parameter <- factor(dl$parameter, levels = par_names_sub)

ggplot(dl, aes(x = par_value, y = ll_mean, color = keep_or_not)) +
  geom_point(alpha = 0.36) +
  facet_wrap(vars(parameter), scales = "free_x") + 
  scale_color_viridis(option = "A", end = 0.7, discrete = TRUE)
# Esto parece bueno

GGally::ggpairs(d3, aes(color = keep_or_not), 
                columns = which(names(d3) %in% par_names_sub))

# model to predict limitness
d3$inside <- 1
d3$inside[d3$keep_or_not == "bad"] <- 0
mm <- mgcv::gam(inside ~ s(intercept, k = 10), data = d3,
                family = "binomial")
plot(mm)
ndat <- data.frame(intercept = seq(min(d3$intercept), max(d3$intercept), length.out = 150),
                   slope = 1)
pp <- predict(mm, ndat, type = "response")
plot(pp ~ ndat$intercept, type = "l")
# parece una buena idea. con esto se pueden filtrar a priori algunas partículas, 
# aquellas con prob < 0.4
# luego se vuelven a filtrar en base a la badness observada, y recién ahí 
# se ajusta el GP.
# quizás, en vez de filtrar duramente convenga simular la inclusión o no, para
# que cada tanto caiga alguna partícula por ahí.


# para llevar los umbrales de tamaño a una escala relativa aplicable a cada fuego,
# puedo calcular el diffsize mínimo y máximo para cada uno, usando partículas
# que no quemen nada o que quemen todo. 
# luego, los umbrales permitidos serían:


# diff_size_upper <- 0.9 * diff_size_max
# diff_size_lower <- 0.9 * diff_size_min

# en realidad lo mismo podría hacer con una combinación de overlap y signo del 
# diffsize. solo debería conocer los límites del overlap para cada fuego.

# queda la pregunta de cómo definir una varianza entre partículas que sea
# razonable e invariable entre fuegos.

# será más rápido usar la proporción del paisaje quemada por los fuegos simulados?
# quemada considerando solo lo quemable.
# Eso sería independiente del fuego observado, pero variable entre paisajes.
# mejor una medida que contemple el paisaje y el fuego observado.





# BASURERO ------------------------------------------------------------------

library(randtoolbox)
a <- sobol(5, seed = 1)
b <- sobol(5, seed = 1, init = FALSE)

a <- sobol(10, seed = 1, init = TRUE)
b <- sobol(20, seed = 1, init = TRUE)

all.equal(a, b[1:10])

cc <- sobol(10, seed = 1)
all.equal(c(a,b), cc)

# get the first n_pw particles
parts_01 <- particles_all[1:20, ]
ll_array <- similarity_simulate_parallel(particles_all[1:n_pw, ]) # just 1 min


# fit gaussian process to particles

# subset rep and metric
rep_id <- 1
obs_ids <- data_id$sim[data_id$rep == rep_id]
metric <- "overlap_sp"

# ll_local <- log(ll_array[obs_ids, , metric, ])
ll_local <- log(ll_array[1, , metric, ])

# get loglik sums
# ll_sums <- apply(ll_local, which(names(dimnames(ll_local)) %in% c("simulation", "particle")), sum)
ll_sums <- ll_local

# get mean and var by particle
mv <- apply(ll_sums, which(names(dimnames(ll_sums)) %in% c("particle")),
            FUN = function(x) c("mean" = mean(x), "var" = var(x))) %>% t %>% 
      as.data.frame()

# join with particle data
part_local <- particles_all_raw[as.numeric(rownames(mv)), 
                                -which(colnames(particles_all_raw) == "fwi")] %>% 
              as.data.frame()
names(part_local) <- par_names_sub

# plots

for(i in 1:ncol(part_local)) {
  par(mfrow = c(1, 2))
  plot(mv$mean ~ part_local[, i], main = names(part_local)[i])
  plot(mv$var ~ part_local[, i], main = names(part_local)[i])
  par(mfrow = c(1, 1))
}

mv$expmean <- exp(mv$mean)
for(i in 1:ncol(part_local)) {
  par(mfrow = c(1, 2))
  plot(mv$expmean ~ part_local[, i], main = names(part_local)[i])
  plot(mv$var ~ part_local[, i], main = names(part_local)[i])
  par(mfrow = c(1, 1))
}


# With narrow priors, the loglik does not seem to be bounded, so estimation can
# proceed over the whole particles (use prior_sd = 2 to see how the loglik is 
# bounded)

plot(mean ~ var, mv) # tanto muchos que ajustan muy bien como muchos que ajustan
                     # muy mal tienen varianza baja


# GP fit

# colnames(p) <- colnames(w0)
# d0 <- as.data.frame(p)
# d0$y <- rnorm(nrow(d0), 
#               (10) * d0$intercept + (-10) * d0$intercept ^ 2,
#               sd = 1)
# d0$noise <- exp((-3) * d0$intercept + (3) * d0$intercept ^ 2)
# 
# # plot(noise ~ intercept, d0)
# 
# des <- d0[, "intercept", drop = F]

form <- formula(~ intercept + I(intercept ^ 2) + 
                  subalpine + I(subalpine ^ 2) + 
                  wet + I(wet ^ 2) + 
                  dry + I(dry ^ 2) +
                  aspect + I(aspect ^ 2) +
                  wind + I(wind ^ 2) +
                  elevation + I(elevation ^ 2) +
                  slope + I(slope ^ 2)) 

gp_01 <- km(formula = form, design = part_local, 
            response = mv$mean, noise.var = mv$var)

# set the coef.trend at a very low value (log(1e-50))

gp_02 <- km(formula = ~1, coef.trend = log(1e-50),
            design = part_local, 
            response = mv$mean, noise.var = mv$var)

gp_02 <- km(formula = ~1, coef.trend = log(1e-50),
            design = part_local, 
            response = mv$mean, noise.var = mv$var)

gp_03 <- km(formula = ~1, coef.trend = -1000, #min(mv$mean) * 2,
            design = part_local, 
            response = mv$mean, noise.var = mv$var)


gp <- gp_03

# Get fitted values
fitted <- predict(gp, newdata = part_local, type = "UK")
# do not use the trend, it doesn't make sense. It seems that adding a quadratic
# effect is not helping.

pred_fool <- part_local
pred_fool[, "intercept"] <- seq(-4, 4, length.out = nrow(part_local))
for(i in 2:8) pred_fool[, i] <- 1

fitted_fool <- predict(gp, newdata = pred_fool, type = "UK")

plot(mv$mean ~ part_local$intercept, col = rgb(0, 0, 0, 0.1), pch = 19)
lines(fitted_fool$mean ~ pred_fool$intercept, type = "l")
lines(fitted_fool$upper95 ~ pred_fool$intercept, type = "l", col = "blue")
lines(fitted_fool$lower95 ~ pred_fool$intercept, type = "l", col = "blue")

plot(exp(mv$mean) ~ part_local$intercept, col = rgb(0, 0, 0, 0.1), pch = 19)
lines(exp(fitted_fool$mean) ~ pred_fool$intercept, type = "l")
lines(exp(fitted_fool$upper95) ~ pred_fool$intercept, type = "l", col = "blue")
lines(exp(fitted_fool$lower95) ~ pred_fool$intercept, type = "l", col = "blue")

# optimize likelihood
ll_fun <- function(betas) {
  m <- matrix(betas, nrow = 1)
  colnames(m) <- par_names_sub
  m <- as.data.frame(m)
  neg_loglik <- -predict(gp, newdata = m, type = "UK")$mean
  return(neg_loglik)
}

opt_ll_01 <- optim(par = colMeans(part_local), fn = ll_fun)
opt_ll_01

# slice for intercept at MLE
pred_mle <- do.call("rbind", lapply(1:200, function(x) opt_ll_01$par)) %>% as.data.frame
pred_mle$intercept <- seq(min(part_local$intercept), max(part_local$intercept), 
                          length.out = nrow(pred_mle))

fitted_mle <- predict(gp, newdata = pred_mle, type = "UK")

plot(mv$mean ~ part_local$intercept, col = rgb(0, 0, 0, 0.1), pch = 19,
     ylim = c(-20, 1))
lines(fitted_mle$mean ~ pred_mle$intercept, type = "l")
lines(fitted_mle$upper95 ~ pred_mle$intercept, type = "l", col = "blue")
lines(fitted_mle$lower95 ~ pred_mle$intercept, type = "l", col = "blue")

plot(exp(mv$mean) ~ part_local$intercept, col = rgb(0, 0, 0, 0.1), pch = 19,
     ylim = c(0, 2))
lines(exp(fitted_mle$mean) ~ pred_mle$intercept, type = "l")
lines(exp(fitted_mle$upper95) ~ pred_mle$intercept, type = "l", col = "blue")
lines(exp(fitted_mle$lower95) ~ pred_mle$intercept, type = "l", col = "blue")


# Filter new particles

fitted$mean %>% summary()
mv$mean %>% summary()
thres <- 4 # the larger, the more inclusive

fitted_uppers <- fitted$mean + qnorm(0.995) * fitted$sd
fitted_max <- max(fitted$mean)
keep <- which((fitted_max - fitted_uppers) <= thres)

plot(mv$mean ~ part_local$intercept, col = rgb(0, 0, 0, 0.1), pch = 19,
     ylim = c(-20, 1))

points(mv$mean[keep] ~ part_local$intercept[keep], col = rgb(1, 0, 0, 0.5), 
       pch = 19)

abline(v = param_mat[1, 1])

new_parts <- particles_all_raw[]
# hacer un while que consiga más partículas hasta llegar a 1000, siempre 
# aceptando el criterio de plausibilidad

1/10000




aaa <- list(a = 1, b = 8:10)
bbb <- list(c = "nada")
c(aaa, bbb)


x <- matrix(NA, 2, 1)
do.call("cbind",lapply(1:ncol(x), function(x) x + 1))

# KESIGUE -----------------------------------------------------------------

# Implementar filtrado de partículas antes de simular prediciendo 
# estocásticamente si van a ser razonables o no (far from the bounds).

# Aplicar filtrado duro post-simulación, y pedir nuevas simulaciones solo si 
# las nuevas son muy pocas.

# Evaluar qué fixed lower limit sería razonable forzar en los GPs (unconditional
# expectation). Una opción sería un poquito menos que el mínimo posible.
# eso es un problema porque cualquier cosa forzada que le metamos 
# (incluso datos inventados) corrompe el ajuste de la zona interesante.

# Pensaba por qué Pájaro no tuvo este problema, y puede ser que su criterio
# de exclusión haya excuido a las partículas que estaban muy por afuera de
# la zona "evaluable". O sea, una forma bruta de decir que esas partículas
# tienen likelihood 0.

# Quizás sea más sensato ajustar 2 modelos en paralelo:
# * un clasificador (estocástico) que diga si estamos en una región plausible,
# * un GP que modele la likelihood en esa región.

# El MCMC sortearía primero la plausibilidad y dado eso procedería a calcular la 
# loglik.

#### OPCION SUPERADORA ----------------------------------------

# que el primer GP sea nuestro enfocador... Ese será el primer gran filtro.
# pero su filtro es "eliminá a las partículas con likelihood ligeramente mejor
# a la del fuego que quema todo".
# osea, definimos a la likelihood del fuego que quema todo como 
# loglik_huge
# entonces, la condición de plausibilidad que imponemos en el primer paso es
# loglik_hat > loglik_huge + delta
# donde loglik_hat se estima a partir del primer GP.

# La loglik_huge se puede calcular para todas las métricas... hermoso.

# El problema teórico quizás sería que de esta forma estaremos achicando la 
# posterior, pero en la práctica, lo que buscamos, es obtener una posterior que
# nos simule fuegos similares a los observados.

# Con los siguientes GPs, el filtro puede seguir de la misma manera que dice
# Wilkinson 2014.

# Como hay que confiarle mucho a ese GP1, idealmente correrlo con muchas partículas
# (3000), total en los siguientes podemos usar menos.

# OTRA idea: en el GP, para eviar errarle cuando se va de rango, 
# se puede sumar el mínimo posible a los valores de loglik antes de ajustar el modelo.
# esto haría que en las regiones con pocos datos la "loglik" tienda a cero, 
# que al transformarlo, sería el mínimo. 
# [PERO ESO NO ES EXACTAMENTE LO QUE HICE ANTES CUANDO LE PUSE QUE TENDIERA
# AL MÍNIMO?]

loglik_huge_arr <- similarity_simulate_particle(c(1e6, rep(0, 8)))
loglik_huge_arr_fire1 <- loglik_huge_arr[1, , ]

options(scipen = 999)
view(colMeans(loglik_huge_arr_fire1))

# pensá que ese overlap o similarity sería la que genera cierta partícula
# EN PROMEDIO. o sea, no podemos permitir que sea muy mala.

# tarea -------------------------------------------------------------------

# hacer secuencial la evaluación de los GPs para incluir más partículas.
# así no calcula predicciones inutilmente.

# TAREAS ------------------------------------------------------------------

# algunos fuegos simulados reales no propagaron. 
# quitar esos o usar un metodo de estimación que evalue el ajuste a varias
# replicas, ajustan la un GP al promedio de las metricas.


# Problema con rejection --------------------------------------------------

# Un problema del approach de wilkinson es que cuando una partícula es 
# supuestamente muy mala, la rechaza por completo. Eso tiene sentido si
# aproximamos la likelihood conjunta. Pero cuando queremos usar likelihoods 
# separadas, sería más sensato obtener el valor, aunque sea muy bajo, antes de
# decidir si la partícula es rechazada o no. Esto es porque una partícual puede
# ser muy mala para un fuego y muy buena para otros.

# Todo sería más fácil si no tuviéramos el maldito problema de los lower bounds
# en la similitud simulada. 
# Lo que podemos hacer es estimar un modelo cuadrático simple que nos prediga
# cómo decae la likelihood cuando estamos fuera del rango simulable. 

# Quizás una forma grosera sería aproximar una forma normal a la superficie de
# likelihood y estimar los parámetros de esta curva cuadrática desde ahí. Eso
# se haría de la misma forma en que los frecuentistas calculan la vcov.

# Otra es definir un modelo de regresión cuadrático a mano y estimarlo a lo 
# bruto: centro a todos los params en su MLE, cosa de que midan "distancia 
# desde el óptimo", y fuerzo a todos los parámetros cuadráticos a ser negativos.
  
# Leer desde el Bolker o desde McElreath esto de la cuadrature.

# Quizás lo mejor sea ajustar una función con la misma forma de la log-PDF de
# una normal multivariada:
# https://en.wikipedia.org/wiki/Multivariate_normal_distribution

# Otra es volver a intentar ponerle una función cuadrática a la media del GP, 
# pero solo usando los datos con loglik arriba del threshold.

# la cuadrática es menos elegante que la mvn_pdf, pero a la vez sería una única
# función que usamos para juzgar all the particles. 
# chequear qué onda la varianza de ese GP. 
# hay que ajustarlo con muchas muestras, dejarle que tenga muestras malas para
# que la curva sepa bajar. por ahí se pueden filtrar con distinto criterio 
# las partículas de intercept positivo y las de intercept negativo.
# 
# Eso debería ayudar a ajustar bien la curva que baje desde arriba en el lado
# que tiene pocos datos.


# mvn curve to fit --------------------------------------------------------

# x and mu are vectors, Sigma is the vcov matrix
mvn_pdf <- function(x, mu, Sigma) {
  
  k <- length(x)
  
  # density
  prob <- (2 * pi) ^ (-k / 2) * det(Sigma) ^ (-1 / 2) * 
          exp(-(1 / 2) * t(x - mu) %*% solve(Sigma) %*% (x - mu))
  
  # unnormalized density
  prob_unnorm <- exp(-(1 / 2) * t(x - mu) %*% solve(Sigma) %*% (x - mu))
  
  return(c("prob" = prob, "prob_unnorm" = prob_unnorm, 
           "log_prob_unnorm" = log(prob_unnorm)))
}

# Example
mus <- c(1, 2)
k <- length(mus)
cor <- -0.9
rho <- matrix(c(1, cor,
                cor, 1), ncol = k, byrow = TRUE)
sds <- c(0.9, 0.8)
Sigma <- rho
for(i in 1:k) {
  for(j in 1:k) {
    Sigma[i, j] <- sds[i] * sds[j] * rho[i, j]
  }
}

solve(Sigma)

# data for predicting the density
gg <- expand.grid(x1 = seq(-2, 5, length.out = 100),
                  x2 = seq(-2, 5, length.out = 100))
gg$prob <- NA
gg$prob_un <- NA
gg$log_prob <- NA

for(i in 1:nrow(gg)) {
  gg[i, 3:5] <- mvn_pdf(x = as.numeric(gg[i, 1:2]), mu = mus, Sigma = Sigma)
}

ggplot(gg) +
  #geom_tile(aes(x = x1, y = x2, fill = prob)) +
  #geom_contour(aes(x = x1, y = x2, z = prob)) +
  geom_contour_filled(aes(x = x1, y = x2, z = prob)) +
  # scale_fill_gradientn(colors = c("blue", "red")) + 
  ggtitle("Likelihood")

# fit a normal model to the loglik, where the mean follows a log-mvn-pdf 
# function. rho is parameterized at the logit scale, then transformed
# as inv_logit(rho_logit) * 2 - 1. sds are parameterized at the log scale.
# all starting points can be set by obtaining the MLE of the GP and then
# getting the vcov.
# remember to add an intercept to the model.


# más fácil: usar 
?mgcv::rmvn()
logden_real <- mgcv::dmvn(t(as.matrix(gg[, 1:2])), mu = mus, V = Sigma)
logden_max <- mgcv::dmvn(mus, mu = mus, V = Sigma)
max(logden_real)

# fix the maximum logdensity at zero, so the intercept is identifiable
logden_zero <- logden_real - logden_max
max(logden_zero) # OK

# entonces, la func a ajustar es
intercept + zero_top_mvn(x)

zero_top_mvn <- function(x, mu, V) {
  logden <- mgcv::dmvn(x, mu, V)
  logden_max <- mgcv::dmvn(mu, mu, V)
  return(logden - logden_max)
}

# mvn_curve predicts the log-density of a multivariate normal given
# the x value and the parameters. It is not normalized and another intercept
# can be used.
# mu and sd are vectors, and corr is the correlation matrix
 
mvn_curve <- function(x, intercept, mu, sd, corr) {
  
}




# just playing: list with dimensions --------------------------------------



l <- vector(mode = "list", length = 2 * 3)
dim(l) <- c(2, 3)
l[, 3, drop = F]


l <- list(1, "a", "hola", rep(1, 5), NULL , rnorm(10))
dim(l) <- c(2, 3)
l[, 3, drop = F]
l[, 3]
l[2, 3]
l[[2, 3]]
l[[6]]


lcol <- l[, 3, drop = F]
str(lcol)
lcol[[2]]


dimnames(l) <- list(filas = 1:2, columnas = letters[1:3])

str(l)

l[2, 1]

# luego aplicar esto a la lista de fires_real así queda más organizada.