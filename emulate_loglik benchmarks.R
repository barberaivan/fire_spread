# Compare the timing of looping across particles in R to compute the loglik 
# vs looping in C++. The second implies passing R objects to C++ only once, 
# as the fires are simulated in the same landscape.

# Import fires data and use a large one (1999_27j_S).

library(Rcpp)
library(terra)
library(tidyverse)
library(microbenchmark)

# load cpp functions:
sourceCpp("similarity_functions_many_particles.cpp")

# fires data
fires <- readRDS(file.path("..", "fire_spread_data", "landscapes_ig-known_non-steppe.rds"))
fire <- fires[["1999_27j_S"]]

# inputs for fire spread
distances <- rep(30, 8) # sides
distances[c(1, 3, 6, 8)] <- 30 * sqrt(2)

# function to simulate from the prior
prior_sim <- function(mu_int = 0, sd_int = 20, r_01 = 0.05, r_z = 0.15) {
  betas <- c(
    "intercept" = rnorm(1, mu_int, sd_int),   # shrubland logit (reference class)
    "subalpine" = rnorm(1, 0, sd_int),        # veg coefficients
    "wet" = rnorm(1, 0, sd_int),
    "dry" = rnorm(1, 0, sd_int),
    "fwi" = rexp(1, r_z),                     # positive
    "aspect" = rexp(1, r_01),                 # positive (northing)
    "wind" = rexp(1, r_01),                   # positive
    "elevation" = (-1) * rexp(1, r_z),        # negative
    "slope" = rexp(1, r_01)                   # positive
  )
  
  return(betas)
}

# make particles
n_part <- 10
parts <- do.call("rbind", lapply(1:n_part, function(i) prior_sim()))
parts[, 1] <- 1e5 # make all burnable

# Loop in cpp
emu_cpp <- function(parts) {
  loglik_try_many <- emulate_loglik_try_par(
    landscape = fire$landscape[, , 1:7], 
    burnable = fire$landscape[, , "burnable"],
    ignition_cells = fire$ig_rowcol,
    particles = parts,
    wind_layer = which(dimnames(fire$landscape)[[3]] == "wind") - 1,
    elev_layer = which(dimnames(fire$landscape)[[3]] == "elev") - 1,
    distances = distances,
    upper_limit = 1.0,
    
    fire_ref = fire[c("burned_layer", "burned_ids", "counts_veg")],
    n_replicates = 10
  )
  return(loglik_try_many)
}

# loop in r
emu_r <- function(parts) {
  
  # array to fill
  simil_array <- array(NA, dim = c(10, 11, nrow(parts)))
  
  for(j in 1:dim(simil_array)[3]) {
    
    # j = 1
    simil_array[, , j] <- emulate_loglik_try(
      landscape = fire$landscape[, , 1:7],
      ignition_cells = fire$ig_rowcol,
      burnable = fire$landscape[, , "burnable"],
      coef = parts[j, ],
      wind_layer = which(dimnames(fire$landscape)[[3]] == "wind") - 1,
      elev_layer = which(dimnames(fire$landscape)[[3]] == "elev") - 1,
      distances = distances,
      upper_limit = 1.0,
      
      fire_ref = fire[c("burned_layer", "burned_ids", "counts_veg")],
      n_replicates = 10
    )
  }
  
  return(simil_array)
}

# benchmark
mbm <- microbenchmark(
  r_loop = {rr <- emu_r(parts = parts)},
  cpp_loop = {cpp <- emu_cpp(parts = parts)},
  times = 1
)
mbm

identical(rr, cpp)
all.equal(rr, cpp)

# muy raro, más rapido loopear en r que en cpp