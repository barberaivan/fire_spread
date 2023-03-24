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
  
  q <- matrix(NA, nrow(p), ncol(p))
  colnames(q) <- c("intercept", "subalpine", "wet", "dry",
                   "aspect", "wind", "elevation", "slope") # without fwi
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
ig_cumprob <- expand.grid(x = ppoints(2), y = ppoints(2))
# in this way, the third ignition point (following terra ordering) falls in non-
# burnable. edit it.
ig_cumprob$x[3] <- 0.15

ig_rowcol <- rbind(
  row = quantile(0:(nrow(land) - 1), probs = ig_cumprob$y) %>% round,
  col = quantile(0:(ncol(land) - 1), probs = ig_cumprob$x) %>% round
)
ig_cells <- cellFromRowCol(land_raster, 
                           row = ig_rowcol["row", ] + 1,
                           col = ig_rowcol["col", ] + 1)

# plot
# check_rast <- rast_from_mat(land[, , "burnable"], land_raster)
# values(check_rast)[ig_cells] <- 2
#plot(check_rast)

# check_rast2 <- rast_from_mat(land[, , "burnable"], land_raster)
# values(check_rast2)[ig_cells] # OK

# Prior predictive check --------------------------------------------------

# Obtain a prior distribution that avoids fires that burn the whole landscape
# or more. The intercept is the most important parameter in this regard.

ig_rowcol <- matrix(c(round(ncol(land) * 0.4), round(nrow(land) * 0.3))) - 1

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

# Choose parameter vectors by hand, visually inspecting a few simulations
# from every one to check that fire sizes are OK with this
# landscape.

# pp <- prior_sim(mu_int = -1.0, sd_int = 0.02, sd_veg = 0.02,
#                 r_z = 1.5, r_01 = 1.1) 
# pp
# fire_prior_sim(pp)
# 
# n_rep <- 30
# # param_mat <- matrix(NA, n_rep, 9)
# param_mat[30, ] <- pp
# saveRDS(param_mat, "files/param_mat.rds")
param_mat <- readRDS("files/param_mat.rds")

# Simulation variance is small when at least a few covariates show large effects.
# Consider, for example, elevation. A large negative effect implies very high
# burn probability in all the lowlads (below average elevation).

curve(dexp(x, 0.5), to = 20); abline(h = 0, col = "red")
