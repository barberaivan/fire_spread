library(FireSpread)
library(terra)
library(tidyverse)
library(testthat)

# r functions to compare
similarity_dir <- file.path("..", "FireSpread", "tests", "testthat", "R_similarity_functions.R")
spread_dir <- file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R")
source(similarity_dir)
source(spread_dir)

# constants ---------------------------------------------------------------

n_coef <- 5

# fire spread parameters (coef)
coefs <- rnorm(n_coef, 0, 3)

# landscape raster
size <- 30
n_rows <- size
n_cols <- size
res <- 30
n_cells <- n_rows * n_cols

landscape <- rast(
  ncol = n_cols, nrow = n_rows, res = res, crs = "EPSG:5343",
  nlyrs = n_coef, # no intercept, but wind uses two lawers
  xmin = 0, xmax = res * n_cols, ymin = 0, ymax = res * n_rows,
  names = c("vfi", "tfi", "elev", "wdir", "wspeed")
)

# fill data
landscape$vfi <- rnorm(ncell(landscape))
landscape$tfi <- rnorm(ncell(landscape))
landscape$elev <- runif(ncell(landscape), 0, 2200)
landscape$wdir <- runif(ncell(landscape), 0, 2 * pi) # radians
landscape$wspeed <- abs(rnorm(ncell(landscape), 0, 2))

burnable <- matrix(rbinom(n_rows * n_cols, size = 1, prob = 0.15),
                   n_rows, byrow = T)

# provide vegetation to compute all metrics
n_veg_types <- 6
vegetation <- matrix(sample(0:(n_veg_types-1), n_rows * n_cols, replace = TRUE),
                     n_rows, byrow = T)

# turn landscape into matrix
predmat <- values(landscape)

# spread_onepix test ------------------------------------------------------

# start
burning_cell <- sample(1:ncell(landscape), size = 1)
neighs_raw <- adjacent(landscape, 1, "queen") %>% as.numeric
not_na <- !is.na(neighs_raw)
pos <- (1:8)[not_na]
neighs_id <- neighs_raw[not_na]

s <- round(runif(1, 100, 1000)) # define seed
set.seed(s) # set seed to compare the random simulation for burning
spread_result_r <- spread_one_cell_r(
  landscape_burning = predmat[burning_cell, , drop = T], # USE DROP!!
  landscape_neighbour = predmat[neighs_id[1], , drop = T],
  coef = coefs,
  position = pos[1],
  upper_limit = 1
)
spread_result_r

set.seed(s)
spread_result_cpp_burn <- spread_one_cell(
  landscape_burning = predmat[burning_cell, , drop = T], # USE DROP!!
  landscape_neighbour = predmat[neighs_id[1], , drop = T],
  coef = coefs,
  position = pos[1] - 1,
  upper_limit = 1
)

spread_result_cpp_prob <- spread_one_cell_prob(
  landscape_burning = predmat[burning_cell, , drop = T], # USE DROP!!
  landscape_neighbour = predmat[neighs_id[1], , drop = T],
  coef = coefs,
  position = pos[1] - 1,
  upper_limit = 1
)

spread_result_r["burn"]; spread_result_cpp_burn
spread_result_r["probs"]; spread_result_cpp_prob
# OK, the seed seems to work.


# simulate_fire test ------------------------------------------------------

# sample parameters
coefs <- rnorm(n_coef)
coefs[1] <- 100
# make ignition point(s)
ig_cell <- sample(1:ncell(landscape), 1)
ig_location <- rowColFromCell(landscape, ig_cell) %>% t
steps <- 2 # check this

# set.seed(s)
burn_result_r <- simulate_fire_deterministic_r(
  landscape = landscape, # SpatRaster
  burnable = burnable,
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  steps,
  plot_animation = F
)

# set.seed(s)
burn_result_cpp <- simulate_fire_deterministic(
  landscape = land_cube(landscape), # SpatRaster
  burnable = burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0,
  steps
)

res_cpp <- rast_from_mat(burn_result_cpp, landscape)
res_r <- rast_from_mat(burn_result_r, landscape)

par(mfrow = c(1, 2))
plot(res_cpp, col = c("green", "black"), main = "C++")
plot(res_r, col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))
# Perfect.



# Test predictors effects ------------------------------------------------

# simulate landscapes where one layer varies and check that results are OK.

# ignition
ig_location <- matrix(rep(round(size / 2), 2), 2, 1)
landscape_base <- landscape
values(landscape_base) <- 0
all_burnable <- matrix(rep(1, ncell(landscape)), n_rows, n_cols)

#______________________________

# wind test
wind_dir <- 35               # direction, in angles, from which the wind comes
lands_sub <- landscape_base
lands_sub$wdir <- rep(wind_dir * pi / 180, ncell(lands_sub))
lands_sub$wspeed <- 8
coefs <- c(-1, rep(0, n_coef - 1))
coefs[b_wind] <- 0.5

seed <- round(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = all_burnable,
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  landscape = land_cube(lands_sub),
  burnable = all_burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))
# the seed works!


#______________________________

# slope test
lands_sub <- landscape_base
lands_sub$elev <- rep(seq(800, 900, length.out = nrow(lands_sub)),
                      each = ncol(lands_sub))
coefs <- c(-1, rep(0, n_coef - 1))
coefs[b_slope] <- 30

seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = all_burnable,
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  landscape = land_cube(lands_sub),
  burnable = all_burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))


#______________________________

# vfi test
lands_sub <- landscape_base
lands_sub$vfi <- rep(seq(-3, 3, length.out = nrow(lands_sub)),
                     each = ncol(lands_sub))
coefs <- c(-1, rep(0, n_coef - 1))
coefs[b_vfi] <- 1

seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = all_burnable,
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  landscape = land_cube(lands_sub),
  burnable = all_burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

#______________________________

# tfi test
lands_sub <- landscape_base
lands_sub$tfi <- rep(seq(-3, 3, length.out = nrow(lands_sub)),
                     each = ncol(lands_sub))
coefs <- c(-1, rep(0, n_coef - 1))
coefs[b_tfi] <- runif(1, 1, 4)

seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = all_burnable,
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  landscape = land_cube(lands_sub),
  burnable = all_burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))


#______________________________

# burnable test
lands_sub <- landscape_base

burnable <- matrix(rbinom(n_rows * n_cols, size = 1, prob = 0.3),
                   n_rows, byrow = T)

coefs <- rnorm(n_coef)

seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = burnable,
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  landscape = land_cube(lands_sub),
  burnable = burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0
)

par(mfrow = c(1, 3))
plot(rast_from_mat(burnable, lands_sub), col = c("gray", "orange"), main = "Burnable")
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))


#______________________________

# steps test
lands_sub <- landscape_base
coefs <- c(1000, rep(0, n_coef - 1))
# steps <- sample(1:10, size = 1)
steps <- 0
seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = all_burnable,
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  steps,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  landscape = land_cube(lands_sub),
  burnable = all_burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0,
  steps
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"),
     main = paste("C++,", steps, "steps"))
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"),
     main = paste("C++,", steps, "steps"))
par(mfrow = c(1, 1))

## extra

simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = all_burnable,
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  steps = 0,
  plot_animation = T
)

# check steps used
cc_used <- simulate_fire_compare(
  landscape = land_cube(lands_sub),
  burnable = all_burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0,
  steps
)
cc_used$steps_used # Ok, in 17 steps burns everything.

# Notes on spread simulation ----------------------------------------------

# Simulations differ between c++ and R when extreme probability values are
# generated, because the precision of floats makes a difference.


# Similarity functions checks ---------------------------------------------


burnable <- matrix(rbinom(n_rows * n_cols, size = 1, prob = 0.9),
                   n_rows, byrow = T)
# provide vegetation to compute all metrics
n_veg_types <- 6
vegetation <- matrix(sample(0:(n_veg_types-1), n_rows * n_cols, replace = TRUE),
                     n_rows, byrow = T)
# reasonable coefs
coefs <- c(-1.5, rnorm(n_coef - 1))

# simulate fires
set.seed(1)
fire_1 <- simulate_fire_compare_veg(
  landscape = land_cube(landscape), # use the array
  burnable = burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0,
  vegetation = vegetation,
  n_veg_types = n_veg_types
)

set.seed(1)
fire_1_ <- simulate_fire_compare_veg(
  landscape = land_cube(landscape), # use the array
  burnable = burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0,
  vegetation = vegetation,
  n_veg_types = n_veg_types
)

set.seed(2)
fire_2 <- simulate_fire_compare_veg(
  landscape = land_cube(landscape), # use the array
  burnable = burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0,
  vegetation = vegetation,
  n_veg_types = n_veg_types
)

# a fire against itself, simulated with the same seed
identical(fire_1, fire_1_)

# different fires, r and cpp functions
similarity_cpp_1_2 <- compare_fires_try(fire_1, fire_2)
similarity_r_1_2 <- compare_fires_r(fire_1, fire_2)

(similarity_cpp_1_2 - similarity_r_1_2) < 1e-6

# a fire against itself, cpp function
similarity_1_1 <- compare_fires_try(fire_1, fire_1)
similarity_1_1

# compare overlap sp
ov_cpp_1_2 <- overlap_spatial(fire_1, fire_2)
ov_r_1_2 <- overlap_spatial_r(fire_1, fire_2)
abs(ov_cpp_1_2 - ov_r_1_2) < 1e-6

overlap_spatial(fire_1, fire_1)
overlap_spatial_r(fire_1, fire_1)


# Overlap plots -----------------------------------------------------------



# reasonable coefs
coefs <- c(-1.5, rnorm(n_coef - 1))

# simulate fires
fire_1 <- simulate_fire_compare(
  landscape = land_cube(landscape), # use the array
  burnable = burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0
)

fire_2 <- simulate_fire_compare(
  landscape = land_cube(landscape), # use the array
  burnable = burnable,
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0
)

common <- fire_1$burned_layer + fire_2$burned_layer

r1 <- rast_from_mat(fire_1$burned_layer, landscape)
r2 <- rast_from_mat(fire_2$burned_layer, landscape)
un <- rast_from_mat(pmax(fire_2$burned_layer, fire_1$burned_layer), landscape)
com <- rast_from_mat(pmin(fire_2$burned_layer, fire_1$burned_layer), landscape)

par(mfrow = c(2, 2))
plot(r1, col = c("green", "red"), main = "Fire 1")
plot(r2, col = c("green", "blue"), main = "Fire 2")
plot(com, col = c("green", "purple"), main = "Intersection")
plot(un, col = c("green", "black"), main = "Union")
par(mfrow = c(1, 1))

overlap_spatial(fire_1, fire_2)
overlap_spatial_r(fire_1, fire_2)