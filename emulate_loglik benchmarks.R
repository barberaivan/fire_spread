# Compare the timing of looping across particles in R to compute the loglik 
# vs looping in C++. The second implies passing R objects to C++ only once, 
# as the fires are simulated in the same landscape.

# Import fires data and use a large one (1999_27j_S).

library(Rcpp)
library(terra)
library(tidyverse)
library(microbenchmark)
library(foreach)
library(doParallel)
library(doRNG)

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
  times = 3
)
mbm
# muy raro, ligeramente más rapido loopear en r que en cpp.

identical(rr, cpp)
all.equal(rr, cpp)

# compare RAM
pr <- peakRAM::peakRAM(
  {rr <- emu_r(parts = parts)},
  {cpp <- emu_cpp(parts = parts)}
)
pr # very similar


# Compare functions when paralelized ----------------------------------------

# wrapper functions for evaluate_loglik, taking one particle or many.

# one particle
emulate_loglik_one <- function(fire, particle) {
  loglik <- emulate_loglik_try(
    landscape = fire$landscape[, , 1:7], 
    burnable = fire$landscape[, , "burnable"],
    ignition_cells = fire$ig_rowcol,
    coef = particle,
    wind_layer = which(dimnames(fire$landscape)[[3]] == "wind") - 1,
    elev_layer = which(dimnames(fire$landscape)[[3]] == "elev") - 1,
    distances = distances,
    upper_limit = 1.0,
    
    fire_ref = fire[c("burned_layer", "burned_ids", "counts_veg")],
    n_replicates = 10
  )
  return(loglik)
}

# many
emulate_loglik_many <- function(fire, particles_matrix) {
  loglik <- emulate_loglik_try_par(
    landscape = fire$landscape[, , 1:7], 
    burnable = fire$landscape[, , "burnable"],
    ignition_cells = fire$ig_rowcol,
    particles = particles_matrix,
    wind_layer = which(dimnames(fire$landscape)[[3]] == "wind") - 1,
    elev_layer = which(dimnames(fire$landscape)[[3]] == "elev") - 1,
    distances = distances,
    upper_limit = 1.0,
    
    fire_ref = fire[c("burned_layer", "burned_ids", "counts_veg")],
    n_replicates = 10
  )
  return(loglik)
}

# Loop over every particle by core. for every particle the data should be copied
# to every worker.
loglik_eval_par1 <- function(fire, particles_list) {
  # particles must be a list!
  result <- foreach(pp = particles_list) %dopar% {
    emulate_loglik_one(fire, pp)
  }
  return(result)
}

# every core takes a matrix of particles. In this way the data is copied only 
# once to each worker.The particles_mats is a list containing the matrix of 
# particles for every core. Its length should be equal to n_cores
loglik_eval_par2 <- function(fire, particles_mats) {
  # particles must be a list!
  result <- foreach(pp = particles_mats) %dopar% {
    emulate_loglik_many(fire, pp)
  }
  return(result)
}


# Set up cluster
n_cores <- 10
cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)
registerDoRNG(seed = 123)

# make particles (all with high intercept, so all fires are the same)
n_part <- n_cores * 10
p_base <- c(1e5, rep(0, 8))
p_list_single <- lapply(1:n_part, function(x) p_base)
n_part_core <- n_part / n_cores
p_list_mat <- lapply(1:n_cores, function(x) {
  matrix(p_base, nrow = n_part_core, ncol = length(p_base), byrow = TRUE)
})

mbm_par <- microbenchmark(
  one = {test1 <- loglik_eval_par1(fire, p_list_single)},
  many = {test2 <- loglik_eval_par2(fire, p_list_mat)},
  times = 5
); mbm_par

# Ram start at 5 Gb
# one: peaks at ~75 % (with 10 cores) and remains quite steady.
# dan lo mismo.