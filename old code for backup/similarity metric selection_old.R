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
# function (n_waves = 5?). 
# In this model we will exclude the FWI effect and use a fixed-intercept. 
# After the last wave, we will obtain the MLE using optim (n_rep * n_met MLEs).
# As for every parameter replicate (n_rep) there will be n_sim fires, the
# similarity function used will be the average between a each simulated fire
# and reference fire.

# For each parameter separately we will compute the squared difference between
# the real value and the MLE. We will choose the metric that
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

# As simulations performed on certain particles will probably be required for
# other parameter vectors and metrics, I save the results from every particle 
# in a simulations_vault, which is written to disk as it becomes large. I
# tried to use the memoise package, but it was not easy to save the result of
# a parallelized function.

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)

library(Rcpp)
library(terra)
library(abind)         # manage arrays

library(randtoolbox)   # sobol sequences
library(DiceKriging)   # Fit Gaussian Processes // probably not used
library(GauPro)
library(mgcv)          # fit spline to choose non-bounded particles
library(microbenchmark)

library(foreach)       # parallelization
library(doMC)          # doRNG had problems with cpp functions

library(FireSpread)    # spread and similarity functions

# Multicore settings -----------------------------------------------------

n_cores <- parallel::detectCores()
registerDoMC(n_cores)

# Data and constants -----------------------------------------------------

# landscape to run the simulations
land_full <- readRDS(file.path("data", "focal fires data",
                               "landscapes_ig-known_non-steppe", "2015_53.rds"))
terrain <- land_full$terrain   # cambiar "land" por terrain
vegetation <- land_full$vegetation
rows <- nrow(vegetation)
cols <- ncol(vegetation)

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
  NA, dim = c(n_met + 2, n_rep, 2, 1), # metrics + size + edge
  dimnames = list(
    metric = c(metric_names, "size_diff", "edge"),
    param_rep = 1:n_rep, 
    mean_var = c("mean", "var"),
    particle = 0
  )
)
# initialize the number of particles with which simulations_vault was written
# to disk
saved_with <- 0

# formula for in_bounds model (gam)
bounds_model_formula <- formula(
  in_bounds ~ 
    s(shrubland, k = 8, bs = "cr") + 
    s(subalpine, k = 8, bs = "cr") + 
    s(wet, k = 8, bs = "cr") +     
    s(dry, k = 8, bs = "cr") + 
    s(northing, k = 8, bs = "cr") +
    s(elev, k = 8, bs = "cr") + 
    s(wind, k = 8, bs = "cr") + 
    s(slope, k = 8, bs = "cr")
  )

loglik_formula <- formula(
  ll ~ 
    # additive effects
    s(shrubland, k = 7, bs = "cr") + 
    s(subalpine, k = 7, bs = "cr") + 
    s(wet, k = 7, bs = "cr") +     
    s(dry, k = 7, bs = "cr") + 
    s(northing, k = 7, bs = "cr") +
    s(elev, k = 7, bs = "cr") + 
    s(wind, k = 7, bs = "cr") + 
    s(slope, k = 7, bs = "cr") + # this estimates 6 parameters (k-1)
    
    # interactions
    ti(shrubland, subalpine, k = 3) + # this estimates 4 params (k-1)^2
    ti(shrubland, wet, k = 3) +
    ti(shrubland, dry, k = 3) +
    ti(shrubland, northing, k = 3) +
    ti(shrubland, elev, k = 3) +
    ti(shrubland, wind, k = 3) +
    ti(shrubland, slope, k = 3) +
    
    ti(subalpine, wet, k = 3) +
    ti(subalpine, dry, k = 3) +
    ti(subalpine, northing, k = 3) +
    ti(subalpine, elev, k = 3) +
    ti(subalpine, wind, k = 3) +
    ti(subalpine, slope, k = 3) +
    
    ti(wet, dry, k = 3) +
    ti(wet, northing, k = 3) +
    ti(wet, elev, k = 3) +
    ti(wet, wind, k = 3) +
    ti(wet, slope, k = 3) +
    
    ti(dry, northing, k = 3) +
    ti(dry, elev, k = 3) +
    ti(dry, wind, k = 3) +
    ti(dry, slope, k = 3) +
    
    ti(northing, elev, k = 3) +
    ti(northing, wind, k = 3) +
    ti(northing, slope, k = 3) +
    
    ti(elev, wind, k = 3) +
    ti(elev, slope, k = 3) +
    
    ti(elev, wind, k = 3) +
    ti(wind, slope, k = 3)
)

loglik_formula_additive <- formula(
  ll ~ 
    # additive effects
    s(shrubland, k = 11, bs = "cr") + 
    s(subalpine, k = 11, bs = "cr") + 
    s(wet, k = 11, bs = "cr") +     
    s(dry, k = 11, bs = "cr") + 
    s(northing, k = 11, bs = "cr") +
    s(elev, k = 11, bs = "cr") + 
    s(wind, k = 11, bs = "cr") + 
    s(slope, k = 11, bs = "cr") # this estimates 6 parameters (k-1)
)


# Escala y centralidad de las predictoras
terr_vec0 <- apply(terrain, 3, as.vector)
veg_vec <- as.vector(vegetation)
terrain_vec <- terr_vec0[veg_vec < 99, ]
summary(terrain_vec)
# wind y veg claramente no están centradas en cero.


# Functions ---------------------------------------------------------------

# prior simulator for graphical prior predictive checks
prior_sim <- function(mu_int = 0, sd_int = 20, r_01 = 0.05, r_z = 0.15,
                      n_veg_types = 4) {
  
  intercepts <- rnorm(n_veg_types, mu_int, sd_int)
  names(intercepts) <- vegetation_names
  
  slopes <- c(
    "northing" = rexp(1, r_01),               # positive 
    "elev" = (-1) * rexp(1, r_z),             # negative
    "wind" = rexp(1, r_01),                   # positive
    "slope" = rexp(1, r_01)                   # positive
  )
  
  return(c(intercepts, slopes))
}

# Prior distribution inverse cumulative distribution function, to transform sobol 
# sequences in [0, 1] ^ d to the prior scale
# p is a matrix of n_particles * d sobol samples of cumulative probability.
prior_icdf <- function(p, n_veg_types = 4, 
                       mu_int = 0, sd_int = 4, r_z = 1, r_01 = 0.75) {
  
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

# plot fires
fire_plot <- function(burned_layer) {
  mat <- vegetation
  mat[,] <- 0
  mat[vegetation < 99] <- 1
  mat[burned_layer == 1] <- 2
  plot_rast <- rast_from_mat(mat, land_raster)
  plot(plot_rast, col = c("black", "forestgreen", "orange"),
       type = "classes", levels = c("non-burnable", "unburned", "burned"))
}

# Function to simulate and plot a few fires to choose parameters that make sense
# in the focal landscape.
fire_prior_sim <- function(prior = NULL, maxcell = 500000) {
  
  sizes <- numeric(3)
  par(mfrow = c(1, 3))
  
  for(i in 1:3) {
    fire_prior <- simulate_fire(
      terrain = terrain, 
      vegetation = vegetation,
      ignition_cells = ig_rowcol,#land_full$ig_rowcol,#
      coef = prior,
      n_veg_types = n_veg_types,
      upper_limit = upper_limit
    )
    # plot
    mat <- vegetation
    mat[,] <- 0
    mat[vegetation < 99] <- 1
    mat[fire_prior == 1] <- 2
    plot_rast <- rast_from_mat(mat, land_raster)
    plot(plot_rast, maxcell = maxcell,
         col = c("black", "forestgreen", "orange"),
         type = "classes", levels = c("non-burnable", "unburned", "burned"),
         main = paste("mean intercept =", 
                      round(mean(prior[vegetation_names]), 3)))
    
    sizes[i] <- sum(fire_prior)
  }
  
  par(mfrow = c(1, 1))
  
  return(sizes)
}

# summarize with mean and variance
summ_mean_var <- function(x) c("mean" = mean(x, na.rm = TRUE), 
                               "var" = var(x, na.rm = TRUE))

# evaluando los bordes
edge_count <- function(burned_layer, r = rows, c = cols) {
  rows_edge <- sum(burned_layer[c(1, r), ])
  cols_edge <- sum(burned_layer[-c(1, r), c(1, c)])
  return(rows_edge + cols_edge)
}

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
    NA, dim = c(n_rep, n_sim, n_sim, length(c(metric_names, "size_diff", "edge"))),
    dimnames = list(
      param_rep = 1:n_rep, # p index
      fire_obs = 1:n_sim, # o index
      simulation = 1:n_sim,# i index
      metric = c(metric_names, "size_diff", "edge")
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
      for(o in which(fires_use[, p])) { # loop over "observed" fires with decent size
        comp <- compare_fires_try(fire_sim, fires_real[[p]][[o]])
        
        # subset metrics to evaluate
        m_keep <- which(names(comp) %in% metric_names) 
        metrics_raw[p, o, i, 1:n_met] <- comp[m_keep]
                                  # remove non-spatial raw metrics (2 to 5)
        # add size_diff
        metrics_raw[p, o, i, "size_diff"] <- sum(fire_sim$counts_veg) - 
                                             sum(fires_real[[p]][[o]]$counts_veg)
        
        # count number of pixels burned in the edge of the landscape
        # to detect escaping fires
        metrics_raw[p, o, i, "edge"] <- edge_count(fire_sim$burned_layer)
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
  metrics <- apply(met_meanvar, c("metric", "param_rep", "mean_var"), mean, 
                   na.rm = TRUE)

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
# simulations_vault array, subsetting the parameters vector (rep_id) and 
# the desired metric(s). It merges this with the particles values, allowing
# to fit a model to the similarity.
# The variance and the size_diff are used only to filter out bounded particles.
# Returns a df where each row is a particle, and with the following columns:
#   ids: particle id,
#   rep_id: id of the parameter vector that simulated the "observed" fires. 
#   ll: mean of the metric across the n_sim simulations (at logit scale!),
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
    size_diff = metrics_array["size_diff", rep_id, "mean", ],
    edge =  metrics_array["edge", rep_id, "mean", ]
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
in_bounds <- function(data, rep_id = 1, metric = "overlap_sp",
                      prop = 0.85, var_factor = 100,
                      edge_prop_upper = 0.05,
                      edge_abs_upper = NULL) {
  
  # # test
  # data <- particles_data_new
  # var_factor = 10
  # # end testo
  
  size_diff_max <- similarity_bounds["size_diff", rep_id, "mean", "largest"]
  size_diff_min <- similarity_bounds["size_diff", rep_id, "mean", "smallest"]
    
  size_diff_lower <- prop * size_diff_min
  size_diff_upper <- prop * size_diff_max
  
  # get variance in metric
  low_var <- c(
    similarity_bounds[metric, rep_id, "var", "largest"],
    similarity_bounds[metric, rep_id, "var", "smallest"]
  ) %>% max
    
  # keep <- as.numeric(
  #   (
  #     (data$size_diff >= size_diff_lower) &  # size in bounds
  #     (data$size_diff <= size_diff_upper)
  #   ) & (
  #     data$var >= low_var * var_factor # non-zero variance
  #   )
  # )
  lowest_var <- low_var * var_factor #0.001 # 
  
  edge_upper <- ifelse(
    is.null(edge_abs_upper),
    edge_prop_upper * edge_large,
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

# Simulate which particles are to be evaluated, based on the the 
# probability of the loglik to be above a certain threshold, as predicted by the
# loglik_model.
simulate_high_loglik <- function(loglik_model, particle_ids, 
                                 loglik_threshold = qlogis(0.15)) {
  
  ## TESTO
  # particle_ids = 1:2000
  # loglik_threshold = qlogis(0.03)
  
  if(any(class(loglik_model) == "km")) {
    preds <- predict(loglik_model, particles_all[particle_ids, ], se.compute = TRUE,
                     type = "UK", light.return = TRUE)
    
    # get probability for x > loglik_threshold
    prob <- 1 - pnorm(loglik_threshold, preds$mean, preds$sd)
  }
  
  if(any(class(loglik_model) == "gam")) {
    nd <- as.data.frame(particles_all[particle_ids, ])
    names(nd) <- par_names
    mu <- predict(loglik_model, nd, type = "response")
    ss <- loglik_model$scale %>% sqrt
    # get probability for x > loglik_threshold
    prob <- 1 - pnorm(loglik_threshold, mu, ss)
  }
  
  # predict in_range_prob
  use_bin <- rbinom(length(prob), size = 1, prob = prob)
  
  return(use_bin)
}

# Function to get more data passing bound and optionally loglik filters. 
#   first, it predicts which particles are good based on the filtering, based
#   on simulations (bernoulli) according to the probability of the filtering
#   models (bounds and loglik). Once n_pw 
#   particles are selected, fires are simulated. 
#   Once the data is available, it checks how many simulated particles have 
#   actually passed the filters, and continues to simulate more particles until
#   the good realized particles is above n_pw.

# loglik_filter determines whether the loglik model (gp) is used to filter new
# particles in advance, by predicting their probability of high loglik before
# simulating. This is useful to disable its use at the first waves, when many
# particles are out of bounds and the gp is unstable, fitted with little data.
# However, the search for new particles may count as useful only those above
# a certain threshold. In that case, the GP is not used to predict good particles,
# but only particles with high loglik are deemed as useful. Note that if 
# gp_filter is TRUE, loglik_tolerance or loglik_lower (or both) must be provided.

# returns a list, containing the new_data data.frame with the data of all 
# simulated particles, and the id of the last tried particle. It may be larger
# than the largest simulated particle.

# Arguments
# previous wave: result from loglik_update()
# loglik_filter: whether to use or not a loglik model to choose new particles,
# loglik_high: is so, a high loglik value to take as reference in accordance 
#   to the loglik_tolerance.
# loglik_lower: lower limit for good particles. It is usually set at a value
#   that is above out-of-bounds particles.
#   The likelihood criterion for particles to count as useful will be
#   new_loglik >= loglik_threshold,
#   loglik_threshold = max(loglik_high - loglik_tolerance, loglik_lower)
# n_new: number of new useful particles to provide (checked after fires are
#   simulated)
# metric.
# rep_id.
# strict: simulate until n_new particles meet the filters strictly, i.e., post-
#   simulation. This requires a lot of computation, and is not recommended.

# EDITAR ACÁ 
# Evitar que simule hasta encontrar todas cosas válidas.
# ser más laxo, como antes. buscar n_pw partículas buenas a priori, simularlas,
# y devolver, aunque no todas sean efectivamente buenas.

# LISTO, ya comenté lo que lo hacía estricto, pero la documentación lo cuenta
# como si fuera estricto.

more_data <- function(previous_wave = NULL,
                      loglik_filter = TRUE,
                      loglik_high = NULL,
                      loglik_tolerance = NULL,
                      loglik_lower = NULL,
                      n_new = 1000,
                      metric = "overlap_sp",
                      rep_id = 1,
                      strict = FALSE) {
  
  # override metric and rep_id if there is a previous_wave
  if(!is.null(previous_wave)) {
    metric <- previous_wave$metric
    rep_id <- previous_wave$rep_id
  }
  
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
    
    # define likelihood-based threshold for counting particles
    relative_lower <- if(is.null(loglik_tolerance) | is.null(loglik_high)) NULL else {
      loglik_high - loglik_tolerance
    } 
      
    loglik_threshold <- if(is.null(relative_lower) & is.null(loglik_lower)) NULL else {
      max(relative_lower, loglik_lower)
    }
    
    if(loglik_filter & is.null(loglik_threshold)) {
      stop("A loglik criterion is needed to pre-filter particles with a loglik model.")
    }
    
    # extract models to filter particles
    loglik_model <- previous_wave$loglik_model
    in_model <- previous_wave$in_model
    
    # initialize number of added particles and data.
    added <- 0       # new useful particles counter
    new_data <- NULL # new data with all the new simulated particles
                     # (may be larger than n_new because n_new counts the 
                     # useful ones after they are simulated)
    last_particle <- previous_wave$last_particle
    
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
      use_ll <- ifelse(loglik_filter, 
                       simulate_high_loglik(loglik_model, try_ids, loglik_threshold),
                       rep(1, length(try_ids)))
      use_in <- simulate_in_bounds(in_model, try_ids)
      use <- use_ll * use_in
      
      # simulate fire in particles that passed the pre-simulation filters
      if(sum(use) > 0) {
        
        try_sim <- try_ids[use == 1]
        metrics_array <- similarity_simulate_memo(particle_ids = try_sim) 
        
        # tidy array as df, previously summarizing over simulations
        new_data_forward <- extract_focal_loglik(
          metrics_array, metric = metric, rep_id = rep_id
        )
        
        ## NON-STRICT FILTERING
        
        # # count how many particles actually passed the filters
        # high_ll <- ifelse( # determines whether loglik is used to judge simulated particles
        #   !is.null(loglik_threshold),
        #   as.numeric(new_data_forward$ll >= loglik_threshold),
        #   rep(1, length(try_ids))
        # )
        # inside <- in_bounds(new_data_forward)
        # added <- added + sum(high_ll * inside)
        
        # merge the new simulations with the previous ones in this search
        new_data <- rbind(new_data, new_data_forward)
      }
      
      added <- added + sum(use)
      last_particle <- max(try_ids)
    } # end while
  
    return(list(new_data = new_data, last_particle = last_particle))
  }
}

# Functions used to find the maximum of a GAM or a GP
gam_predict <- function(x, loglik_model) {
  nd <- as.data.frame(matrix(x, nrow = 1))
  names(nd) <- par_names
  return(predict(loglik_model, nd, type = "response"))
}

# for DiceKrigging GP
gp_predict_dice <- function(x, loglik_model) {
  nd <- matrix(x, nrow = 1)
  colnames(nd) <- par_names
  predict(loglik_model, type = "UK", se.compute = F, light.return = TRUE,
          newdata = nd)$mean
}

# for GauPro gp
gp_predict <- function(x, loglik_model) {
  nd <- matrix(x, nrow = 1)
  return(loglik_model$pred(nd, se.fit = FALSE))
}


# Function to subset particles to fit a GP.
# It takes the best 200 plus 800 with an even-cumulative probability below 
# them. It returns the particle_ids to use for GP fitting
subset_particles <- function(data, n_best = 200, n_others = 800) {
  data <- data[data$in_bounds == 1, ] # use only in_bounds
  
  if(nrow(data) < (n_best + n_others)) {
    warning("Not so much data, we'll give you what we can.")
  }
  
  if(n_best > nrow(data)) {
    n_best <- nrow(data) 
    warning("n_best was set to nrow(data), considering only in_bounds.")
  }
  
  data <- data[order(data$ll, decreasing = TRUE), ]
  ids_best <- data$particle_ids[1:n_best]
  
  if(nrow(data) > (n_best + n_others)) {
    if((nrow(data) - nbest) < n_others) {
      n_others <- nrow(data) - nbest
      warning(paste("n_others was set to", n_others))
    }
    
    ids_others <- sample(data$particle_ids[(n_best+1):nrow(data)], 
                         size = n_others, replace = F)
  } else {ids_others <- NA}
  
  res <- as.numeric(na.omit(c(ids_best, ids_others)))
  return(res)
}


# loglik_update runs a wave of likelihood emulator fitting.

# returns
# loglik_model: likelihood model fitted in the previous step.
# in_model: gam model predicting the probability of particles laying in_bounds.
# loglik_max_obs, _fitted, and _optim: high loglik values to use as benchmarks
#   for selecting more particles in the following wave.
# params_max_optim: parameter values reaching loglik_max_optim.
# particles_data: data.frame with the following columns
#   ids: particle id,
#   rep_id: parameter vector id that originated the "observed" fires analysed,
#   wave: latest wave in which the particle was used to fit a GP,
#   ll: logit of the mean of the simulated likelihoods in each particle,
#   var: variance of the simulated likelihoods in each particle,
#   parameter_values: matrix with the parameter values.
#   in_bounds: integer [0, 1] indicating if the particles are in_bounds.
# metric: metric to emulate likelihood being evaluated.
# rep_id: id for the real parameter vector to be estimated.

# arguments
# previous_wave: list generated in the previous wave (default is NULL, which 
#   starts the estimation process),
# model_filter: see more_data() (old gp_filter)
# loglik_high: a high value of loglik to take as reference to judge particles.
#   It might be a number or "max_fitted", "max_observed", or "max_estimated".
# loglik_tolerance: new particles will be judged as good if they have high
#   probability of being > (loglik_high - loglik_tolerance)
# loglik_lower: alternatively, the particles are jugded as good if they are
#   > loglik_lower. If present, the loglik_threshold is the maximum between
#   (loglik_high - loglik_tolerance) and loglik_lower,
# n_pw: number of new good particles to include in the new wave,
# fit_ll_model: whether to fit a loglik model or not. It is convenient to wait 
#   until many many simulations in bounds are available (maybe fit it in the 
#   second wave?).
# fit_ll_bounded: use only particles in_bounds to fit the loglik model,
# ll_model = c("gam", "gp")
# ll_model_optimize = should the fitted loglik model be optimized? Not recommended
#   if it's first fitted with many particles out of bounds.
# # metric: similarity metric,
# rep_id: parameter vector id of the fires evaluated.

loglik_update <- function(previous_wave = NULL, 
                          loglik_filter = TRUE,
                          loglik_high = NULL,
                          loglik_tolerance = 3,
                          loglik_lower = NULL,
                          n_pw = 1000, 
                          fit_ll_model = TRUE,
                          fit_ll_bounded = F,
                          use_best_particles = T,
                          n_best = 500,
                          ll_model_optimize = T,
                          metric = "overlap_sp", 
                          rep_id = 1) {
  
  ### TESTO
  # previous_wave = NULL#w1#result#NULL
  # loglik_filter = F
  # loglik_high = NULL# w1$loglik_max_obs
  # loglik_tolerance = NULL
  # loglik_lower = NULL
  # n_pw = 2000
  # fit_ll_model = F#TRUE
  # fit_ll_bounded = T
  # ll_model = "gp"
  # ll_model_optimize = T
  # metric = "overlap_sp"
  # rep_id = 1
   
   
  # previous_wave = w1
  # n_pw = 1000
  # loglik_filter = F
  # fit_ll_model = T
  # fit_ll_bounded = T
  # use_best_particles = T
  # n_best = 800
  # loglik_model = "gp"
  ### TESTO ENDO
  
  # override metric and rep_id if there is a previous_wave
  if(!is.null(previous_wave)) {
    metric <- previous_wave$metric
    rep_id <- previous_wave$rep_id
  }
  
  # Simulate new particles
  print("Getting more data (simulating fires)")
  fresh_ones <- more_data(
    previous_wave = previous_wave,
    loglik_filter = loglik_filter,
    loglik_high = loglik_high,
    loglik_tolerance = loglik_tolerance,
    loglik_lower = loglik_lower,
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
  particles_data_new$in_bounds <- in_bounds(particles_data_new,
                                            rep_id = rep_id,
                                            metric = metric,
                                            prop = 0.8,
                                            var_factor = 100,
                                            edge_prop_upper = 0.01)
  
  ### Chequeo bounds
  # data_plot <- do.call(data.frame, particles_data_new)
  # cols_change <- grep("par_values", names(data_plot))
  # names(data_plot)[cols_change] <- par_names
  # ggplot(data_plot, aes(shrubland, ll, colour = subalpine)) + 
  #   geom_point() +
  #   facet_wrap(vars(in_bounds))
  ### Chequeo bounds
  
  # merge old and new particles data
  if (is.null(previous_wave)) {
    particles_data_join <- particles_data_new
  } else {
    particles_data_join <- rbind(previous_wave$particles_data, 
                                 particles_data_new)
  }
  
  # fit or update in_bounds_model using all the data
  in_model <- fit_in_bounds(particles_data_join)
  
  if(fit_ll_model) {
    # filter (or not) in_bounds particles to fit model
    fit_this_0 <- if(fit_ll_bounded) {
      which(particles_data_join$in_bounds == 1)
    } else 1:nrow(particles_data_join)
    
    # use only the n_best particles to fit the loglik_model

# agregar subset_particles() ----------------------------------------------

    
    # 
    if(use_best_particles) {
      n_in_bounds <- length(fit_this_0)
      cand_sort <- fit_this_0[order(particles_data_join$ll[fit_this_0], 
                                    decreasing = TRUE)]
      fit_this <- if(n_in_bounds >= n_best) cand_sort[1:n_best] else cand_sort
    } else {
      fit_this <- fit_this_0
    }
    
    print(paste("Fitting", ll_model))
    
    # with DiceKrigging    
    # loglik_model <- km(design = particles_data_join$par_values[fit_this, ], 
    #                    response = particles_data_join$ll[fit_this], 
    #                    nugget.estim = TRUE,
    #                    # multistart = ceiling(n_cores / 2),
    #                    control = list(maxiter = 5e4))
    # # max fitted value
    # fitted_ll <- predict(loglik_model, type = "UK",
    #                      newdata = particles_data_join$par_values[fit_this, ],
    #                      se.compute = TRUE, light.return = TRUE)$mean
    
    # with GauPro
    loglik_model <- gpkm(
      X = particles_data_join$par_values[fit_this, ], 
      Z = particles_data_join$ll[fit_this], 
      parallel = F, useC = TRUE, nug.max = 100,
      kernel = "matern52"
    )
    
    # max fitted loglik value
    fitted_ll <- loglik_model$pred(particles_data_join$par_values[fit_this, ], 
                                   se.fit = F)
    loglik_max_fitted <- max(fitted_ll)
    
    # max fitted parameters value
    id_fitted <- particles_data_join$particle_id[fit_this]
    id_max <- id_fitted[which.max(fitted_ll)]
    id_filter <- which(particles_data_join$particle_id == id_max)
    params_max_fitted <- particles_data_join$par_values[id_filter, ]
    
    if(ll_model_optimize) {
      lowers <- apply(particles_data_join$par_values[fit_this, ], 2, min)
      uppers <- apply(particles_data_join$par_values[fit_this, ], 2, max) 
      
      op <- optim(params_max_fitted, gp_predict, loglik_model = loglik_model,
                  control = list(fnscale = -1, maxit = 1e5),
                  method = "L-BFGS-B", lower = lowers, upper = uppers)  
    }
  }
  
  # put all together
  result <- list(
    loglik_model = if(fit_ll_model) loglik_model else NULL,
    in_model = in_model,
    loglik_max_obs = max(particles_data_join$ll),
    loglik_max_fitted = if(fit_ll_model) loglik_max_fitted else NULL,
    loglik_max_optim = if(fit_ll_model & ll_model_optimize) op$value else NULL,
    loglik_optim = if(fit_ll_model & ll_model_optimize) op else NULL,
    params_max_optim = if(fit_ll_model & ll_model_optimize) op$par else NULL,
    particles_data = particles_data_join,
    last_particle = fresh_ones$last_particle,
    metric = metric,
    rep_id = rep_id
  )
  
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
partial_predictions <- function(varying = "shrubland", data, loglik_model, mle) {
  
  ### TEST
  # varying = "shrubland" 
  # data = data_pred 
  # loglik_model = loglik_model; mle = fitting_wave$loglik_max_optim
  ### 
  
  if(varying != "all") {
    new_data <- make_newdata(varying = varying, data = data, mle = mle)
    
    if(any(class(loglik_model) == "km")) {
      pred <- predict(loglik_model, newdata = new_data[, par_names], type = "UK",
                      light.return = TRUE, se.compute = TRUE)
      
      new_data$mle <- pred$mean
      new_data$upper <- pred$upper95
      new_data$lower <- pred$lower95
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
    
    if(any(class(loglik_model) == "km")) {
      pred <- predict(loglik_model, newdata = new_data[, par_names], type = "UK",
                      light.return = TRUE, se.compute = TRUE)
      
      new_data$mle <- pred$mean
      new_data$upper <- pred$upper95
      new_data$lower <- pred$lower95
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
                        varying = "shrubland", # parameter to vary in the 1d plot 
                        color_point = c("wave_plot", "in_bounds"),
                        latest_wave = FALSE) { 
  
  #### TEST
  # fitting_wave <- w1; varying = "all"
  # color_point <- "wave_plot"
  # latest_wave <- F
  #### 
  
  loglik_model <- fitting_wave$loglik_model
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
  mle <- fitting_wave$params_max_optim
  
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
                 mapping = aes_string(x = varying, y = "ll",
                                      color = color_point, shape = color_point),
                 size = 2, alpha = 0.5) + 
      scale_shape_manual(values = c(16, 17)) +
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
      ylab("loglik") + 
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
                 mapping = aes_string(x = "varying_val", y = "ll",
                                      color = color_point, shape = color_point),
                 size = 2, alpha = 0.5) + 
      scale_shape_manual(values = c(16, 17)) +
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

# most elevation values are below 0, so the elevation effect will be correlated
# with the intercepts
# elev_vec <- terrain[, , "elev"] %>% as.vector()
# veg_vec <- vegetation %>% as.vector()
# hist(elev_vec)
# hist(elev_vec[veg_vec < 99])
# mean(elev_vec[veg_vec < 99])
# sum(elev_vec[veg_vec < 99] < 0) / length(elev_vec[veg_vec < 99])
# 0.8287471
# como el cero está en el percentil 0.82, sería esperable que haya alta corr
# entre los intercepts y la elevation. Es recomendable ajustar luego la MVN
# en la escala centrada.

# param_mat <- matrix(NA, n_rep, d)
# colnames(param_mat) <- par_names
# pp <- prior_sim(mu_int = 1, sd_int = 0.05, r_z = 2, r_01 = 1.5)
# pp
# fire_prior_sim(pp, maxcell = 90000)
# print("cuidado")

# param_mat[10, ] <- pp
# saveRDS(param_mat, "files/param_mat.rds")
# param_mat <- readRDS("files/param_mat.rds")
# parameterization without FWI and no vegetation class as reference

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


# identify very small fires
fire_sizes <- do.call("cbind", lapply(fires_real, function(par) {
  do.call("c", lapply(par, function(f) {
    f$counts_veg %>% sum
  }))
})) %>% as.data.frame

# get size limits for estimable fires. Before simulate large and small fires.

coef_burn_all <- c(rep(1e6, n_veg_types), rep(0, length(terrain_names)))
coef_burn_none <- c(rep(-1e6, n_veg_types), rep(0, length(terrain_names)))

small_fire <- simulate_fire(
  terrain = terrain,
  vegetation = vegetation,
  ignition_cells = ig_rowcol,
  coef = coef_burn_none,
  n_veg_types = n_veg_types,
  upper_limit = upper_limit
)

large_fire <- simulate_fire(
  terrain = terrain,
  vegetation = vegetation,
  ignition_cells = ig_rowcol,
  coef = coef_burn_all,
  n_veg_types = n_veg_types,
  upper_limit = upper_limit
)

# define limits

decent_large <- 0.2 * as.vector(large_fire) %>% sum
decent_small <- 111 # smallest in our data set (10 ha)

table(as.numeric(as.matrix(fire_sizes)) >= decent_small)
table(as.numeric(as.matrix(fire_sizes)) <= decent_large)

# matrix indicating which fires have a proper size to be estimated.
# Fires too large or too small do not discard parameters that make fires
# out of bounds (not spreading or burning everything).
# This introduces bias in the estimation, but the point here is to find
# the metric that takes us closest to the true value.
# This matrix is used in similarity_simulate_particle to ignore too large and
# too small fires.
fires_use <- fire_sizes
for(i in 1:n_rep) {
  fires_use[, i] <- fire_sizes[, i] >= decent_small &
                    fire_sizes[, i] <= decent_large
}
# colSums(fires_use)

# Get bounds on similarity metrics ----------------------------------------

# (quizás tenga algunas ideas desactualizadas)

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
p <- 0.10
similarity_min <- 0.07
(similarity_threshold = similarity_min + p * (1 - similarity_min))

# With this threshold we can emulate the likelihood using a GP, and filter out
# particles < similarity_threshold. After this first astringent wave, we can 
# proceed with the method proposed by Wilkinson (2014).

# Ideally, the first wave should be estimated with many particles, so mistakes 
# are not performed there. Then, GPs could use data in a cumulative manner, so 
# the last one is a good predictor of the likelihood above similarity_threshold.
# (This allows to compute the likelihood instead of just rejecting the particle).

simil_lower_large <- similarity_simulate_particle(coef_burn_all, 
                                                  n_sim = 10)
simil_lower_small <- similarity_simulate_particle(coef_burn_none, 
                                                  n_sim = 10)

# junto large y small en un array
similarity_bounds <- abind::abind(simil_lower_large, simil_lower_small,
                                  along = length(dim(simil_lower_large)) + 1)
dimnames(similarity_bounds)[[length(dim(similarity_bounds))]] <- c("largest", "smallest")
names(dimnames(similarity_bounds)) <- c(
  "metric", "rep_id", "mean_var", "extreme"
)



similarity_bounds["overlap_sp", 1, "mean", "largest"] %>% qlogis()

# evaluando los bordes
edge_large <- edge_count(large_fire)


# Long sobol sequence -----------------------------------------------------

# for all estimations to run over the same particles (when possible), a very long
# sobol sequence will be produced, expecting not to need further "simulations".
# A larger sequence will be created from scratch if more particles are needed,
# to ensure the same particles are always used.
n_large <- 300000
particles_all <- particles_sim(n_large)

# Exploring waves --------------------------------------------------------

w1 <- loglik_update(n_pw = 1500, use_best_particles = F, fit_ll_bounded = F,
                    fit_ll_model = F,
                    rep_id = 8)
# loglik_plot(w1, varying = "all", "in_bounds")

w2 <- loglik_update(w1, n_pw = 1000, loglik_filter = F,
                    fit_ll_model = T, fit_ll_bounded = T,
                    use_best_particles = T,
                    n_best = 800,
                    rep_id = 8)
# p2 <- loglik_plot(w2, varying = "all", "in_bounds")

w3 <- loglik_update(w2, n_pw = 1000, loglik_filter = T,
                    loglik_high = w2$loglik_max_obs,
                    loglik_tolerance = 1,
                    fit_ll_model = T, fit_ll_bounded = T,
                    use_best_particles = T,
                    n_best = 800,
                    rep_id = 8)
# p3 <- loglik_plot(w3, varying = "all", "in_bounds")
# p3$plot + ylim(-3, 1)

w4 <- loglik_update(w3, n_pw = 1000, loglik_filter = T,
                    loglik_high = w3$loglik_max_obs,
                    loglik_tolerance = 1,
                    fit_ll_model = T, fit_ll_bounded = T,
                    use_best_particles = T,
                    n_best = 800,
                    rep_id = 8)
p4 <- loglik_plot(w4, varying = "all", "in_bounds")

# OK, en 4 olas estima muy bien, al menos el máximo.

# lo que sigue ------------------------------------------------------------

# meter loglik_update en un while() o en un for, para hacer una función que 
# estime la loglik de un saque y devuelva el resultado de la última ronda.

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

# quizás seleccionar partículas que sean >= max_obs - 0.25 * (max_obs - min_obs) 
# más abajo del máximo.
# Por otro lado, usar solo las in_bounds para ajustar el GAM lo hace mucho más 
# económico.
# VAMOS CON GAM Y NO CON GP para buscar partículas. 
# A LO SUMO SE PUEDE VER SI VALE LA PENA AL FINAL USAR
# GP EN VEZ DE BAYESIAN S-MVN, pero será en otro momento.
# También hay que probar estimar el map del S-MVN en vez de muestrear la posterior.



# probar cuándo esto puede tener sentido. quizás usando siempre el criterio de
# la tolerancia funciona mejor.


# El more_data tarda muchísimo si le pedimos que todo lo que devuelva sean 
# partículas efectivamente buenas. Además, devuelve mucha porquería, muchos datos.
# Mejor no pedirle tanto.

# updating the km takes too long, it's better to fit from scratch

# El gam es super rápido, pero se va al upite fuera de rango. Donde no hay datos
# no se aplana sino que se va a la mierda.

# CORRELACION en las likelihood -------------------------------------------

# conviene ajustarlas en escala centrada, o sea, centrando las predictoras a ni
# vel paisaje. Así las mvn tendrán menos corr entre intercept y el resto. 




# TAREA -------------------------------------------------------------------

# Comparar GauPro, que parece estar escrito en C++ una parte.
# Quizás corra mejor y sea más preciso.
# función gpkm()
dd <- w4$particles_data
dd <- dd[order(dd$ll, decreasing = TRUE), ]
dd <- dd[1:2000, ]

X <- dd$par_values
Z <- dd$ll

system.time(
gp1 <- gpkm(X, Z, parallel = F, useC = TRUE)
)
# 28 s con n 800
# 402.382 s con n 2000
# quizás valga la pena

# gp1 <- gpkm(X, Z, parallel = T, useC = TRUE, parallel_cores = 8)
# en paralelo tira error

system.time(
gp2 <- km(design = X, 
          response = Z, 
          nugget.estim = TRUE,
          control = list(maxiter = 5e4))
)
# 43 s con n 800
# 494.626 con n 2000
# lo feo de este es que verbosea y que no sabemos si converge o no.


# GauPro es ligeramente más rápido, no verbosea, y aparentemente no se detiene
# en las 101 iterations.
 

# Cambiar todo por GauPro -------------------------------------------------
# y eliminar las opciones de gam para estimar la loglik.
# para el optimizador, usar gauPro$predict_C
