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

n_veg <- 3 # vegetation types
n_coef <- n_veg + 5 # it doesn't include steps
n_layers <- n_coef - n_veg + 1
layer_names <- c("vegetation", "ndvi", "north", "elev", "wdir", "wspeed")

# landscape raster
size <- 30
n_rows <- size
n_cols <- size
res <- 30
n_cells <- n_rows * n_cols

landscape <- rast(
  ncol = n_cols, nrow = n_rows, res = res, crs = "EPSG:5343",
  nlyrs = n_layers,
  xmin = 0, xmax = res * n_cols, ymin = 0, ymax = res * n_rows,
  names = layer_names
)

# fill data
landscape$vegetation <- sample(0:(n_veg-1), n_rows * n_cols, replace = TRUE)
landscape$ndvi <- rnorm(ncell(landscape))
landscape$north <- runif(ncell(landscape))
landscape$elev <- rnorm(ncell(landscape))
landscape$wdir <- runif(ncell(landscape), 0, 2 * pi) # radians
landscape$wspeed <- abs(rnorm(ncell(landscape), 0, 2))

# turn landscape into matrix
predmat <- values(landscape) # first column is vegetation

# spread_onepix test ------------------------------------------------------

# start

coefs <- rnorm(n_coef, 0, 3)

burning_cell <- sample(1:ncell(landscape), size = 1)
neighs_raw <- adjacent(landscape, 1, "queen") %>% as.numeric
not_na <- !is.na(neighs_raw)
pos <- (1:8)[not_na]
neighs_id <- neighs_raw[not_na]

s <- round(runif(1, 100, 1000)) # define seed
set.seed(s) # set seed to compare the random simulation for burning
spread_result_r <- spread_one_cell_r(
  vegetation = predmat[neighs_id[1], 1, drop = T],
  landscape_burning = predmat[burning_cell, -1, drop = T], # USE DROP!!
  landscape_neighbour = predmat[neighs_id[1], -1, drop = T],
  coef = coefs,
  position = pos[1],
  upper_limit = 1
)
spread_result_r

set.seed(s)
spread_result_cpp_burn <- spread_one_cell(
  vegetation = predmat[neighs_id[1], 1, drop = T],
  landscape_burning = predmat[burning_cell, -1, drop = T], # USE DROP!!
  landscape_neighbour = predmat[neighs_id[1], -1, drop = T],
  coef = coefs,
  position = pos[1] - 1,
  upper_limit = 1
)

spread_result_cpp_prob <- spread_one_cell_prob(
  vegetation = predmat[neighs_id[1], 1, drop = T],
  landscape_burning = predmat[burning_cell, -1, drop = T], # USE DROP!!
  landscape_neighbour = predmat[neighs_id[1], -1, drop = T],
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
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  steps,
  plot_animation = F
)

# set.seed(s)
burn_result_cpp <- simulate_fire_deterministic(
  landscape = land_cube(landscape)[, , -1], # without vegetation
  vegetation = land_cube(landscape)[, , 1],
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

all.equal(res_cpp, res_r)
# Perfect.

# Test predictors effects ------------------------------------------------

# simulate landscapes where one layer varies and check that results are OK.

# ignition
ig_location <- matrix(rep(round(size / 2), 2), 2, 1)
landscape_base <- landscape
values(landscape_base) <- 0

#______________________________

# wind test
wind_dir <- 90               # direction, in angles, from which the wind comes
lands_sub <- landscape_base
lands_sub$wdir <- rep(wind_dir * pi / 180, ncell(lands_sub))
lands_sub$wspeed <- 8
coefs <- c(-1, rep(0, n_coef - 1))
coefs[b_wind] <- 0.5

seed <- round(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub,
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  landscape = land_cube(lands_sub)[, , -1],
  vegetation = land_cube(lands_sub)[, , 1],
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
coefs[b_slope] <- -10

seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  landscape = land_cube(lands_sub)[, , -1],
  vegetation = land_cube(lands_sub)[, , 1],
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))


#______________________________

# elevation test
lands_sub <- landscape_base
lands_sub$elev <- rep(seq(-3, 3, length.out = nrow(lands_sub)),
                      each = ncol(lands_sub))
coefs <- c(-2, rep(0, n_coef - 1))
coefs[b_elev] <- 3

seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  landscape = land_cube(lands_sub)[, , -1],
  vegetation = land_cube(lands_sub)[, , 1],
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

#______________________________

# ndvi test
lands_sub <- landscape_base
lands_sub$ndvi <- rep(seq(-1, 1, length.out = nrow(lands_sub)),
                      each = ncol(lands_sub))
coefs <- c(-1, rep(0, n_coef - 1))
coefs[b_ndvi] <- runif(1, -4, -1)

seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  landscape = land_cube(lands_sub)[, , -1],
  vegetation = land_cube(lands_sub)[, , 1],
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))


#______________________________

# vegetation test
lands_sub <- landscape_base
lands_sub$vegetation <- sample(c(0, 1, 2, 99), size = ncell(lands_sub),
                               prob = c(1, 1, 1, 2), replace = TRUE)

coefs <- rep(0, n_coef)
coefs[1:n_veg] <- 30

seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  landscape = land_cube(lands_sub)[, , -1],
  vegetation = land_cube(lands_sub)[, , 1],
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0
)

par(mfrow = c(1, 3))
plot(lands_sub[["vegetation"]], col = c("black", "black", "black", "orange"), main = "Burnable")
plot(rast_from_mat(cc, lands_sub), col = c("green", "red"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "red"), main = "R")
par(mfrow = c(1, 1))


#______________________________

# steps test
lands_sub <- landscape_base
coefs <- rep(0, n_coef)
coefs[1:n_veg] <- 1000
steps <- sample(1:nrow(lands_sub), size = 1)
# steps <- 0
seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  steps,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  landscape = land_cube(lands_sub)[, , -1],
  vegetation = land_cube(lands_sub)[, , 1],
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
  ignition_cells = ig_location,
  coef = coefs,
  upper_limit = 1.0,
  steps = 0,
  plot_animation = T
)

# check steps used
cc_used <- simulate_fire_compare(
  landscape = land_cube(lands_sub)[, , -1],
  vegetation = land_cube(lands_sub)[, , 1],
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


landscape$vegetation <- sample(c(0, 1, 2, 99), n_rows * n_cols, replace = TRUE)

# reasonable coefs
coefs <- rnorm(n_coef)
coefs[1:n_veg] <- -1.5

# simulate fires
set.seed(1)
fire_1 <- simulate_fire_compare_veg(
  landscape = land_cube(landscape)[, , -1], # use the array
  vegetation = land_cube(landscape)[, , 1],
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0,
  n_veg_types = n_veg
)

set.seed(1)
fire_1_ <- simulate_fire_compare_veg(
  landscape = land_cube(landscape)[, , -1], # use the array
  vegetation = land_cube(landscape)[, , 1],
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0,
  n_veg_types = n_veg
)

set.seed(2)
fire_2 <- simulate_fire_compare_veg(
  landscape = land_cube(landscape)[, , -1], # use the array
  vegetation = land_cube(landscape)[, , 1],
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0,
  n_veg_types = n_veg
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
coefs <- rnorm(n_coef)
coefs[1:n_veg] <- -1

# simulate fires
fire_1 <- simulate_fire_compare(
  landscape = land_cube(landscape)[, , -1], # use the array
  vegetation = land_cube(landscape)[, , 1],
  ignition_cells = ig_location - 1,
  coef = coefs,
  upper_limit = 1.0
)

fire_2 <- simulate_fire_compare(
  landscape = land_cube(landscape)[, , -1], # use the array
  vegetation = land_cube(landscape)[, , 1],
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