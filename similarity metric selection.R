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

library(tidyverse)
library(Rcpp)
library(terra)

library(randtoolbox)   # sobol sequences
library(DiceKriging)   # Fit Gaussian Processes
library(microbenchmark)

library(foreach)       # parallelization
library(doParallel)
library(doRNG)         # reproducible parallelization

sourceCpp("spread_functions.cpp")
sourceCpp("similarity_functions.cpp")

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
  
  q[, "elevation"] <- qexp(p[, "elevation"], rate = r_z)
  
  # fill null fwi parameter
  q[, "fwi"] <- 0
  
  return(q)
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
fire_prior_sim <- function(prior = prior_sim()) {
  
  sizes <- numeric(3)
  
  par(mfrow = c(1, 3))
  
  for(i in 1:3) {
    fire_prior <- simulate_fire_cpp(
      landscape = land[, , 1:7],
      burnable = land[, , "burnable"],
      ignition_cells = ig_rowcol,
      coef = pp,
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
         main = paste("intercept =", round(pp["intercept"], 3)))
    
    sizes[i] <- sum(fire_prior)
  }
  
  par(mfrow = c(1, 1))
  
  return(sizes)
}


# Fire simulator as a function of the particle (fixed landscape)
# just a loop over for simulate_fire_compare()
emulate_loglik_particle <- function(particle, n_sim = 10) {
  
  metrics <- array(NA, dim = c(n_obs, n_sim, n_met))
  
  for(o in 1:n_obs) {
    for(i in 1:n_sim) {
      fire_sim <- simulate_fire_compare(
        landscape = land[, , 1:7],
        burnable = land[, , "burnable"],
        ignition_cells = ig_rowcol,
        coef = particle, # without id
        wind_layer = which(dnames == "wind") - 1,
        elev_layer = which(dnames == "elev") - 1,
        distances = distances,
        upper_limit = upper_limit
      )
      
      metrics[i, o, ] <- compare_fires_try(fire_sim, fires_real[[o]])
    }
  }
  
  
  #### SEGUIR ACÁ
  
  return(sims)
}

# Fire simulator to simulate over a list of particles. It returns a list as long
# as the evaluated particles, where the names of the list have the id of the
# particle in the long sobol sequence. This is used to avoid recomputing 
# particles.

cl <- makeCluster(parallel::detectCores(), type = "FORK")
registerDoParallel(cl)
registerDoRNG(seed = 123)

emulate_loglik_parallel <- function(particle_mat) {
  
  # particle_mat <- particles_all[1:20, ]
  # prepare particles and ids
  ids <- particle_mat[, 1]
  particles_list <- lapply(1:nrow(particle_mat), function(x) particle_mat[x, -1])
  
  # simulate in parallel
  result <- foreach(pp = particles_list) %dopar% {
    simulate_fire_compare_particle(pp)  ## USAR EMULATE_LOGLIK_PARTICLE
  }
  
  # set ids
  names(result) <- ids
  
  return(result)
}


# Data and constants -----------------------------------------------------

# landscape to run the simulations
land_full <- readRDS(file.path("..", "fire_spread_data", "focal fires data",
                              "landscapes_ig-known_non-steppe", "2015_53.rds"))
land <- land_full$landscape

dnames <- dimnames(land)[[3]]

# terra-raster of the same landscape (with elevation) for plotting purposes.
# only the raster structure is needed, not its data.
land_raster <- rast(file.path("..", "fire_spread_data", "focal fires data",
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


# Simulate fires ----------------------------------------------------------

# n_rep (parameters) * n_sim (fires) will be simulated in the same landscape. 
seeds <- matrix(100:(100 + n_rep * n_sim - 1),
                nrow = n_sim, ncol = n_rep)

data_id <- expand.grid(sim = 1:n_sim, 
                       rep = 1:n_rep)
n_obs <- nrow(data_id)

fires_real <- vector(mode = "list", length = n_obs)

for(i in 1:n_obs) {
    
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


# Long sobol sequence -----------------------------------------------------

# for all estimations to run over the same particles (when possible), a very long
# sobol sequence will be produced, expecting not to need further "simulations".
n_large <- 300000

particles_all_raw <- prior_q(sobol(n_large, dim = 8, seed = 123), # without FWI
                         mu_int = -1.0, sd_int = 0.05, sd_veg = 0.05,
                         r_z = 1.3, r_01 = 1.0)
# this prior is widened relative to the one used to simulate fires

# add particle ID to avoid recomputing fires
particles_all <- cbind(part_id = 1:nrow(particles_all_raw), particles_all_raw)

# Wave 1 ------------------------------------------------------------------

# get the first n_pw particles
parts_01 <- particles_all[1:20, ]

part_list <- lapply(1:nrow(parts_01), function(x) parts_01[x, -1])

fires_01 <- simulate_fires_parallel(particles_all[1:n_pw, ])

# compute 200 * 1000 DISCREPANCIES
# 
# 