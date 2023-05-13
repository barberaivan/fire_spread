# Code to choose the similarity metric to approximate the likelihood, based on
# simulated fires. 
# We will use the fire 2015_53, "Alerces 2015", which is small enough to reduce 
# the computational burden. 

# The distribution (%) of vegetation type in this landscape is

# shrubland   subalpine   wet          dry
# 41.02752    23.69743    13.09468     22.18038

# and pixel counts are
# 103823      59968       33137        56129

# dry forest here is cypress, and there are no more vegetation types
# for now.

# n_rep parameter vectors will be simulated from the prior (n_rep = 15?).
# For each one n_sim fire events are going to be simulated, called "observed" 
# hereafter (n_sim = 10), using 1 fixed ignition point. By fixing the ignition 
# point, the simulated fires will be comparable with all the n_rep * n_sim 
# observed ones.

# For every parameter vector (n_rep) and similarity metric (n_met = 7?) a 
# Gaussian Process will be fitted in n_waves waves to emulate the likelihood 
# function. 
# In this model we will exclude the FWI effect and use a fixed-intercept. 
# After the last wave, we will obtain the MLE using optim (n_rep * n_met MLEs).
# As for every parameter replicate (n_rep) there will be n_sim fires, the
# similarity function used will be the average between a each simulated fire
# and reference fire.

# For each parameter separately we will compute the squared difference between
# the real value and the posterior mean. We will choose the metric that
# makes the smaller difference in most parameters, probably with a preference
# for the simplest one (spatial overlap) if differences are small. 

# Overall, n_rep * n_met estimations will be carried out. In the first 
# wave all instances will start with the same particles, so the simulated fires
# at particles_01 will not be saved (only the similarity metrics will be). 
# After the first wave, the following particles to be evaluated will depend on
# the first GP, so the particles used will probably differ across estimations. 

# In this code the simulation variance in similarty metrics is not considered, 
# i.e., I just take the mean across all similarities between simulated fires
# and observed fires corresponding to a {parameter_vector, particle, metric} 
# set.

# TRY TO IMPLEMENT MEMOISED FUNCTION similarity_simulate_parallel(). If possible,
# ignore the following:

# From the second wave on, simulated fires by particle will be written to disk,
# and when the particle is required, the fires will be loaded instead of 
# simulated to compute the discrepancy metrics.

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
library(memoise)       # memoization of similarity_simulate_parallel

sourceCpp("spread_functions.cpp")
sourceCpp("similarity_functions.cpp")


# Multicore settings -----------------------------------------------------

n_cores <- parallel::detectCores()
registerDoMC(n_cores)

# Data and constants -----------------------------------------------------

# landscape to run the simulations
land_full <- readRDS(file.path("data", "focal fires data",
                               "landscapes_ig-known_non-steppe", "2015_53.rds"))
terrain <- land_full$terrain   # cambiar "land" por terrain
vegetation <- land_full$vegetation

dnames <- dimnames(terrain)[[3]]

# terra-raster of the same landscape (with elevation) for plotting purposes.
# only the raster structure is needed, not its data.
land_raster <- rast(file.path("data", "focal fires data",
                              "wind ninja files", "2015_53.tif"))

# constants for fire spread simulation

upper_limit <- 0.5

# ignition points
ig_rowcol <- matrix(c(round(ncol(terrain) * 0.4), round(nrow(terrain) * 0.3))) - 1

# number of particles to choose by wave
n_pw <- 1000 

# number of real parameters to simulate
n_rep <- 10

# number of fires to simulate by parameter
n_sim <- 10

# total number of real fires
n_obs <- n_rep * n_sim

rep_obs_data <- data.frame(rep = rep(1:10, each = n_sim),
                           obs = 1:n_obs)

# parameters stuff
n_veg_types <- 4
vegetation_names <- c("shrubland", "subalpine", "wet", "dry")
terrain_names <- c("northing", "elev", "wind", "slope")
par_names <- c(vegetation_names, terrain_names)
d <- length(par_names)

# number of similarity metrics to compare
n_met <- 7

metric_names <- c("overlap_sp",
                  "sp_norm_5050", "sp_norm_7525",
                  "sp_expquad_5050", "sp_expquad_7525",
                  "sp_quad_5050", "sp_quad_7525")

# empty array to store simulations (manual memoization)
simulations_vault <- array(
  NA, dim = c(n_met + 1, n_rep, 2, 1),
  dimnames = list(
    metric = c(metric_names, "size_diff"),
    param_rep = 1:n_rep, 
    mean_var = c("mean", "var"),
    particle = 0
  )
)
# initialize the number of particles with which simulations_vault was written
# to disk
saved_with <- 0

# Functions ---------------------------------------------------------------

# prior simulator for graphical prior predictive checks
prior_sim <- function(mu_int = 0, sd_int = 20, r_01 = 0.05, r_z = 0.15,
                      n_veg_types = 4) {
  
  intercepts <- rnorm(n_veg_types, 0, sd_int)
  names(intercepts) <- vegetation_names
  
  slopes <- c(
    "northing" = rexp(1, r_01),                 # positive 
    "wind" = rexp(1, r_01),                   # positive
    "elev" = (-1) * rexp(1, r_z),        # negative
    "slope" = rexp(1, r_01)                   # positive
  )
  
  return(c(intercepts, slopes))
}

# Prior distribution inverse cumulative distribution function, to transform sobol 
# sequences in [0, 1] ^ d to the prior scale
# p is a matrix of n_particles * d sobol samples of cumulative probability.
prior_icdf <- function(p, n_veg_types = 4, 
                       mu_int = -1.0, sd_int = 4, r_z = 1.3, r_01 = 1.0) {
  
  q <- matrix(NA, nrow(p), ncol(p))
  colnames(q) <- par_names
  colnames(p) <- colnames(q)
  
  # fill matrix with parameters samples
  q[, 1:n_veg_types] <- qnorm(p[, 1:n_veg_types], mean = mu_int, sd = sd_int)
  
  # transform [-1, 1] variables
  names_01 <- which(colnames(q) %in% c("northing", "wind", "slope"))
  q[, names_01] <- qexp(p[, names_01], rate = r_01)
  
  # negative effect for elevation
  q[, "elev"] <- (-1) * qexp(p[, "elev"], rate = r_z)

  return(q)
}

# function to make particles from a sobol sequence
particles_sim <- function(N, d = 8) {
  prior_icdf(sobol(N, dim = d, seed = 123, init = TRUE)) 
}

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
      terrain = terrain, 
      vegetation = vegetation,
      ignition_cells = ig_rowcol,
      coef = prior,
      n_veg_types = n_veg_types,
      upper_limit = upper_limit
    )
    # plot
    mat <- vegetation
    mat[,] <- 0
    mat[vegetation < 99] <- 1
    mat[fire_prior == 1] <- 2
    plot_rast <- rast_from_mat(burnable_mat, land_raster)
    plot(plot_rast, col = c("black", "forestgreen", "orange"),
         type = "classes", levels = c("non-burnable", "unburned", "burned"),
         main = paste("mean intercept =", 
                      round(mean(prior[vegetation_names]), 3)))
    
    sizes[i] <- sum(fire_prior)
  }
  
  par(mfrow = c(1, 1))
  
  return(sizes)
}

# summarize with mean and variance
summ_mean_var <- function(x) c("mean" = mean(x), "var" = var(x))

# Simulate similarity for a given particle on the set of observed fires.
# Returns an array [n_rep, n_met, c("mean", "var")] with the similarity metrics 
# by real parameter vector, after summarizing siilarities across simulated and 
# observed replicates.The variance is used to detect bounded particles.
# It also includes size difference to rule out particles making extremely large 
# fires.
similarity_simulate_particle <- function(particle, n_sim = 10) {
  
  ## testo
  # particle <- particles_sim(N = 1)
  ## end testo
  
  metrics_raw <- array(
    NA, dim = c(n_rep, n_sim, n_sim, n_met + 1),
    dimnames = list(
      param_rep = 1:n_rep, # p index
      fire_obs = 1:n_sim, # o index
      simulation = 1:n_sim,# i index
      metric = c(metric_names, "size_diff")
    )
  )

  for(i in 1:n_sim) { # simulated fires by particle
    fire_sim <- simulate_fire_compare(
      terrain = terrain, 
      vegetation = vegetation,
      ignition_cells = ig_rowcol,
      coef = particle,
      n_veg_types = n_veg_types,
      upper_limit = upper_limit
    )
    
    for(p in 1:n_rep) {  # loop over "real" parameter vectors
      for(o in 1:n_sim) { # loop over "observed" fires
        comp <- compare_fires_try(fire_sim, fires_real[[p]][[o]])
        
        # subset metrics to evaluate
        m_keep <- which(names(comp) %in% metric_names) 
        metrics_raw[p, o, i, 1:n_met] <- comp[m_keep]
                                  # remove non-spatial raw metrics (2 to 5)
        # add size_diff
        metrics_raw[p, o, i, n_met + 1] <- sum(fire_sim$counts_veg) - 
                                           sum(fires_real[[p]][[o]]$counts_veg)                   
      }
    }
  }
  
  # compute mean and variance across simulations within every observed fire.
  met_meanvar <- apply(metrics_raw, c("param_rep", "fire_obs", "metric"), 
                       summ_mean_var)
  # name the new dimension (mean_var), put as first
  names(dimnames(met_meanvar))[1] <- "mean_var"
  
  # summarize mean and variance across observed fires for every parameter
  # replicate
  metrics <- apply(met_meanvar, c("metric", "param_rep", "mean_var"), mean)

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

# memoized version (by hand) of simulate_particles_parallel.
# It stores everything in an array called simulations_vault, and it's written
# to disk every time the number of non-written simulated particles 
# is > save_every.
# It although it saves all particles, it returns an array containing only the 
# requested ones.
similarity_simulate_memo <- function(
    particle_ids = 1:100, save_every = 3000
  ) {

  stored <- dimnames(simulations_vault)[["particle"]] %>% as.numeric()
  compute_this <- particle_ids[!(particle_ids %in% stored)]
  
  # compute new particles (if there are any)
  if(length(compute_this) > 0) {
    new_ones <- similarity_simulate_parallel(compute_this)
  } else new_ones <- NULL
  
  # merge with previously computed ones
  simulations_vault <<- abind(simulations_vault, new_ones, along = 4,
                              use.dnns = TRUE, use.first.dimnames = TRUE)
  
  result <- simulations_vault[, , , as.character(particle_ids)]
  
  # write simulations_vault to disk every 3000 simulations
  n_computed <- dim(simulations_vault)[4]

  if(n_computed > (saved_with + save_every)) {
    saved_with <<- n_computed
    pathh <- file.path("data", "simulations", 
                       "similarity-metric-selection_simulations-vault.rds")
    saveRDS(simulations_vault, pathh)
  }
  
  return(result)
}

# extract_focal_loglik: function to extract similarity values from the 
# simulations_vault array, subseting the parameters vector (rep_id) and 
# the desired metric(s). It merges this with the particles values, allowing
# to fit a model to the similarity.
# The variance and the size_diff are used only to filter out bounded particles.
# Returns a df where each row is a particle, and with the following columns:
#   ids: particle id,
#   rep_id: id of the parameter vector that simulated the "observed" fires. 
#   ll_mean: mean of the metric across the n_sim simulations (at logit scale!),
#   var: variance of the metric across the n_sim simulations, at its original 
#     scale,
#   metric: name of the metric used.
#   size_diff_mean, size_diff_var: mean and variance of the size difference.
extract_focal_loglik <- function(metrics_array, metric = "overlap_sp",
                                 rep_id = 1) {
  
  ### TEST
  # metrics_array <- similarity_simulate_memo(1:10)
  # metric = "overlap_sp"; rep_id = 1; 
  ### TEST
  
  if(!(metric %in% dimnames(metrics_array)[["metric"]])) {
    stop("Metric does not exist.")
  }
  
  # extract values
  ll_summ <- data.frame(
    ll = qlogis(metrics_array[metric, rep_id, "mean", ]),  ## USING LOGIT-LIKELIHOOD
    var = metrics_array[metric, rep_id, "var", ],
    size_diff = metrics_array["size_diff", rep_id, "var", ]
  )
  
  # particles dimension
  particles_dim <- which(names(dimnames(metrics_array)) == "particle")
  
  # metadata
  metadata <- data.frame(
    particle_id = dimnames(metrics_array)[[particles_dim]] %>% as.numeric,
    rep_id = rep_id,
    metric = metric,
    wave = NA
  )
  
  # merge all
  df <- cbind(metadata, ll_summ)
  df$par_values <- particles_all[metadata$particle_id, par_names]

  return(df)
}

# in_bounds: evaluate whether the particles are in_bounds, returning a binary
# vector to fit in_bounds model.
in_bounds <- function(data, rep_id, 
                      prop = 0.95, var_threshold = 0.15) {
  
  # # test
  # particles_data <- particles_data_new
  # size_diff_limits <- data.frame(size_diff_min = rep(-190000, n_rep),
  #                                size_diff_max = rep(190000, n_rep))
  # ll_var_low = 0.1; prop = 0.95
  # # end testo
  
  size_diff_max <- similarity_bounds["size_diff", rep_id, "largest"]
  size_diff_min <- similarity_bounds["size_diff", rep_id, "smallest"]
    
  size_diff_lower <- prop * size_diff_min
  size_diff_upper <- prop * size_diff_max
    
  keep <- as.numeric(
    (
      (data$size_diff_mean >= size_diff_lower) &  # size in bounds
      (data$size_diff_mean <= size_diff_upper)
    ) & (
      data$var >= var_threshold                # non-zero variance
    )
  )
  
  return(as.numeric(filt))  
} 





# fit the in_bounds model, specifying the data and a formula (default is 
# s(intercept, k = 10))
fit_in_bounds <- function(
    particles_data, 
    form = formula(in_bounds ~ s(shrubland, k = 8) + 
                     s(subalpine, k = 8) + 
                     s(wet, k = 8) +     
                     s(dry, k = 8))
    
) {
  
  data <- cbind(in_bounds = particles_data$in_bounds, 
                particles_data$par_values) %>% as.data.frame
  m <- gam(form, data = data, family = binomial(link = "logit"), method = "REML")
  
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
  prob[prob > 1] <- 1
  
  use_bin <- rbinom(length(prob), size = 1, prob = prob)
  
  return(use_bin)
}

# Simulate which particles are to be evaluated, based on the the 
# probability of the loglik to be above a certain threshold, as predicted by the
# gp model.
simulate_high_loglik <- function(gp, particle_ids, 
                                 loglik_threshold = logit(0.15)) {
  
  # predict
  preds <- predict(gp, particles_all[particle_ids, ], se.compute = TRUE,
                   type = "UK", light.return = TRUE)
  
  # get probability for x > loglik_threshold
  prob <- 1 - pnorm(loglik_threshold, preds$mean, preds$sd)
  
  # predict in_range_prob
  use_bin <- rbinom(length(prob), size = 1, prob = prob)
  
  return(use_bin)
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
 
# returns a list, containing the new_data data.frame with the data of all 
# simulated particles, and the id of the last tried particle. It may be larger
# than the largest simulated particle.


# SEGUIR POR ACÁ ----------------------------------------------------------
# testear more_data y loglik_update. Por las dudas, también revisar los plots

more_data <- function(previous_wave = NULL,
                      loglik_tolerance = NULL,
                      loglik_lower = NULL,
                      n_new = 500,
                      metric = "overlap_sp",
                      rep_id = 1) {
  
  # is this the first wave?
  first_wave <- is.null(previous_wave)
  
  if(first_wave) {
    metrics_array <- similarity_simulate_memo(particle_ids = 1:n_new) 
    
    # tidy array as df, previously summarizing over simulations
    new_data <- extract_focal_loglik(
      metrics_array, metric = metric, rep_id = rep_id
    )
    
    return(list(new_data = new_data, last_particle = n_new))
    
  # if there are previous waves, it's not so easy...
  } else {
    
    # check we have a likelihood criterion
    if(is.null(loglik_lower) & is.null(loglik_tolerance)) {
      stop("No likelihood criterion was provided.")
    }
    
    # extract models to filter particles
    gp <- previous_wave$gp
    in_model <- previous_wave$in_model
    
    # define likelihood_threshold
    if(!is.null(loglik_tolerance)) {
      model_lower <- previous_wave$loglik_max_fitted - loglik_tolerance
    } else {
      model_lower <- NULL
    }
    
    loglik_threshold <- max(model_lower, loglik_lower)
    
    # initialize number of added particles and data.
    added <- 0       # new useful particles counter
    new_data <- NULL # new data with all the new simulated particles
                     # (will be larger than n_new because n_new counts the 
                     # useful ones after they are simulated)
    
    while(added < n_new) {
      
      try_ids <- (last_particle + 1) : (last_particle + n_new * 10)
      
      # If you run out of particles, make more
      if(max(try_ids) > nrow(particles_all)) {
        particles_all <<- particles_sim(nrow(particles_all) + 1e5)
        # this is inefficient because recreates the sequence from scratch,
        # but it ensures the particles used are always the same even if
        # the computation is run over separate R sessions.
      }
      
      # simulate filters pass
      use_ll <- simulate_high_loglik(gp, try_ids, loglik_threshold)
      use_in <- simulate_in_bounds(in_model, particle_ids)
      use <- use_ll * use_in
      
      # simulate fire in particles that passed
      if(sum(use) > 0) {
        
        try_sim <- try_ids[use == 1]
        metrics_array <- similarity_simulate_memo(particle_ids = try_sim) 
        
        # tidy array as df, previously summarizing over simulations
        new_data_forward <- extract_focal_loglik(
          metrics_array, metric = metric, rep_id = rep_id
        )
        
        # count how many particles actually passed the filters
        high_ll <- as.numeric(new_data_forward$ll >= loglik_threshold)
        inside <- in_bounds(new_data_forward)
        added <- added + sum(high_ll * inside)
        
        # merge the new simulations with the previous ones in this search
        new_data <- rbind(new_data, new_data_forward)
      }
      
      last_particle <- max(try_ids)
    } # end while
  
    return(list(new_data = new_data, last_particle = last_particle))
  }
}


# loglik_update runs a wave of gp-fitting.

# Arguments

# previous_wave: list generated in the previous wave (default is NULL, which 
#   behaves accordingly for the first wave),
# loglik_tolerance: new particles will be judged as good if they have high
#   probability of being > (loglik_max_fitted - loglik_tolerance)
# loglik_lower: alternatively, the particles are jugded as good if they are
#   > loglik_lower. If present, the loglik_threshold is the maximum between
#   (loglik_max_fitted - loglik_tolerance) and loglik_lower,
# n_pw: number of new good particles to include in the new wave,
# metric: similarity metric,
# rep_id: parameter vector id of the fires evaluated.

# returns

# gp: Gaussian process fitted in the previous step.
# particles_data: data.frame with the following columns
#   ids: particle id,
#   rep_id: parameter vector id that originated the "observed" fires analysed,
#   wave: latest wave in which the particle was used to fit a GP,
#   ll: logit of the mean of the simulated likelihoods in each particle,
#   var: variance of the simulated likelihoods in each particle,
#   parameter_values: matrix with the parameter values.
#   in_bounds: integer [0, 1] indicating if the particles are in_bounds,
loglik_update <- function(previous_wave = NULL,
                          # list from the previous wave
                          loglik_tolerance = 3,
                          loglik_lower = NULL,
                          n_pw = 1500, 
                          metric = "overlap_sp", 
                          rep_id = 1) {
  
  # Simulate new particles
  print("Getting more data (simulating fires)")
  fresh_ones <- more_data(
    previous_wave = previous_wave,
    loglik_tolerance = 3,
    loglik_threshold = loglik_lower,
    n_new = n_pw,
    metric = metric,
    rep_id = rep_id
  )
  
  particles_data_new <- fresh_ones$new_data
  
  # get wave order
  this_wave <- ifelse(is.null(previous_wave), 
                      1, 
                      max(previous_wave$particles_data$wave) + 1)
  
  particles_data_new$wave <- this_wave
  
  # define in_bounds particles (binary)
  particles_data_new$in_bounds <- in_bounds(particles_data_new)
  
  # merge old and new particles data
  if (is.null(previous_wave)) {
    particles_data_join <- particles_data_new
  } else {
    particles_data_join <- rbind(particles_data, particles_data_new)
  }

  # fit or update in_bounds_model using all the data
  in_model <- fit_in_bounds(particles_data_join)
  
  # fit gp using only in_bounds observations
  fit_this <- which(particles_data_join$in_bounds == 1)
  
  print("Fitting GP")
  gp_new <- km(design = particles_data_join$par_values[fit_this, ], 
               response = particles_data_join$ll[fit_this], 
               nugget.estim = TRUE,
               multistart = n_cores,
               control = list(maxiter = 5e4))

  # compute loglik_max_fitted
  fitted <- predict(gp_new, type = "UK",
                    newdata = particles_data_join$par_values[fit_this, ],
                    se.compute = TRUE, light.return = TRUE)

  loglik_max_fitted <- max(fitted$mean)
  
  # get fitted mle
  id_fitted <- particles_data_join$particle_id[fit_this]
  id_max <- id_fitted[which.max(fitted$mean)]
  id_filter <- particles_data_join$particle_id == id_max
  par_max_fitted <- particles_data_join$par_values[id_filter, ]
  
  # optimize the GP
  op <- optim(par_max_fitted, fn = function(x) {
    predict(gp_new, type = "UK", se.compute = F, light.return = TRUE,
            newdata = matrix(x, nrow = 1))$mean
    }, control = list(fn_scale = -1))
  
  # put all together
  result <- list(gp = gp_new,
                 loglik_max_fitted = loglik_max_fitted,
                 loglik_max_optim = op$value,
                 mle = op$par,
                 particles_data = particles_data_join,
                 last_particle = fresh_ones$last_particle)
  
  return(result)
}  





# function to make new data varying only one predictor.
# the mle, if provided, must be named.
make_newdata <- function(varying = "shrubland", data, mle = NULL) {
  
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
partial_predictions <- function(varying = "shrubland", data, gp, mle) {
  
  ### TEST
  # varying = "all"; data = w1$particles_data; gp = w1$gps[[length(w1$gps)]]
  # mle = mle
  ### 
  
  if(varying != "all") {
    new_data <- make_newdata(varying = varying, data = data, mle = mle)
    pred <- predict(gp, newdata = new_data[, par_names], type = "UK",
                    light.return = TRUE, se.compute = TRUE)
    
    new_data$mle <- pred$mean
    new_data$upper <- pred$upper95
    new_data$lower <- pred$lower95
  } else {
    new_data <- do.call("rbind", lapply(par_names, function(x) {
      make_newdata(varying = x, data = data, mle = mle)
    }))
    
    pred <- predict(gp, newdata = new_data[, par_names], type = "UK",
                    light.return = TRUE, se.compute = TRUE)
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
                    filter_bounds = TRUE) { 
  
  #### TEST
  # fitting_wave <- w1; varying = "intercept"; fixed_at = "mean"
  # varying = "all"; local_range = FALSE
  #### 
  
  gp <- fitting_wave$gp
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
  #   opt <- optim(colMeans(as.matrix(data_pred[, par_names])), fn = ll)
  #   mle <- opt$par
  #   names(mle) <- par_names
  # }
  
  # get mle for partial predictions
  mle <- fitting_wave$mle
  
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
                 mapping = aes_string(x = varying, y = "ll",
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
                 mapping = aes_string(x = "varying_val", y = "ll",
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
# remove fwi
param_mat <- param_mat[, -5] 
# remove reference category in intercept
param_mat[, 2:3] <- param_mat[, 2:3] + param_mat[, 1]

# Simulate real fires ------------------------------------------------------

# n_rep (parameters) * n_sim (fires) will be simulated in the same landscape. 
seeds <- matrix(100:(100 + n_rep * n_sim - 1),
                nrow = n_sim, ncol = n_rep)

fires_real <- vector(mode = "list", length = n_rep)
names(fires_real) <- paste("rep", 1:n_rep, sep = "_")

for(p in 1:n_rep) {
  fires_real_internal <- vector(mode = "list", length = n_sim)
  names(fires_real_internal) <- paste("sim_obs", 1:n_sim, sep = "_")
  
  for(i in 1:n_sim) {
    
    set.seed(seeds[i, p])
      
    fires_real_internal[[i]] <- simulate_fire_compare(
      terrain = terrain, 
      vegetation = vegetation,
      ignition_cells = ig_rowcol,
      coef = param_mat[p, ],
      n_veg_types = n_veg_types,
      upper_limit = upper_limit
    )
  }
  
  fires_real[[p]] <- fires_real_internal
}

# EDITAR ESTO: -------------------------------------------------------------
# (solo un poco)
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



# lo que sigue ------------------------------------------------------------

# meter loglik_update en un while() o en un for, para hacer una función que 
# estime lo loglik de un saque y devuelva el resultado de la última ronda.

# Loopear sobre replicas y sobre metricas


# notas -------------------------------------------------------------------

# El código de ahora es bastante distinto al de antes. No ajusto muchos gps 
# sucesivos, sino que hago uno solo acumulativo usando siempre las partículas
# in_bounds (excepto en la wave 1).

# la función more_data por un lado predice cuáles van a ser buenas, las simula,
# y luego evalua si fueron realmente buenas, y termina recién cuando el n
# de partículas buenas nuevas requeridas es igualado o superado. 
 
# Par high loglik se pueden usar 2 criterios: 
# loglik > loglik_crit
# loglik_crit puede ser provided (llamada loglik_lower), o se puede dar una 
# tolerance y un loglik_max, para calcularla así:
# model_lower = loglik_max - tolerance.
# Then, loglik_crit = max(loglik_lower, model_lower)

# probar cuándo esto puede tener sentido. quizás usando siempre el criterio de
# la tolerancia funciona mejor.
