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
n_terrain <- 3 # it doesn't include steps
n_b_terrain <- n_terrain - 1
n_layers <- n_terrain + 1
layer_names <- c("vegetation", "elev", "wdir", "wspeed")

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
landscape$elev <- rnorm(ncell(landscape), 1500, 300)
landscape$wdir <- runif(ncell(landscape), 0, 2 * pi) # radians
landscape$wspeed <- abs(rnorm(ncell(landscape), 0, 2))

# turn landscape into matrix
predmat <- values(landscape) # first column is vegetation

# spread_onepix test ------------------------------------------------------

# start

cv <- rnorm(n_veg, 0, 3)
ct <- rnorm(n_b_terrain, 0, 3)

burning_cell <- sample(1:ncell(landscape), size = 1)
neighs_raw <- adjacent(landscape, 1, "queen") %>% as.numeric
not_na <- !is.na(neighs_raw)
pos <- (1:8)[not_na]
neighs_id <- neighs_raw[not_na]

s <- round(runif(1, 100, 1000)) # define seed
set.seed(s) # set seed to compare the random simulation for burning
spread_result_r <- spread_one_cell_r(
  vegetation = predmat[neighs_id[1], 1, drop = T],
  terrain_burning = predmat[burning_cell, -1, drop = T], # USE DROP!!
  terrain_neighbour = predmat[neighs_id[1], -1, drop = T],
  coef_veg = cv,
  coef_terrain = ct,
  position = pos[1],
  upper_limit = 1
)
spread_result_r

set.seed(s)
spread_result_cpp_burn <- spread_one_cell(
  vegetation = predmat[neighs_id[1], 1, drop = T],
  terrain_burning = predmat[burning_cell, -1, drop = T], # USE DROP!!
  terrain_neighbour = predmat[neighs_id[1], -1, drop = T],
  coef_veg = cv,
  coef_terrain = ct,
  position = pos[1] - 1,
  upper_limit = 1
)

spread_result_cpp_prob <- spread_one_cell_prob(
  vegetation = predmat[neighs_id[1], 1, drop = T],
  terrain_burning = predmat[burning_cell, -1, drop = T], # USE DROP!!
  terrain_neighbour = predmat[neighs_id[1], -1, drop = T],
  coef_veg = cv,
  coef_terrain = ct,
  position = pos[1] - 1,
  upper_limit = 1
)

spread_result_r["burn"]; spread_result_cpp_burn
spread_result_r["probs"]; spread_result_cpp_prob
# OK, the seed seems to work.


# simulate_fire test ------------------------------------------------------

# sample parameters
cv <- rnorm(n_veg)
ct <- rnorm(n_b_terrain, 0, 3)

# make ignition point(s)
ig_cell <- sample(1:ncell(landscape), 1)
ig_location <- rowColFromCell(landscape, ig_cell) %>% t
steps <- 2 # check this

# set.seed(s)
burn_result_r <- simulate_fire_deterministic_r(
  landscape = landscape, # SpatRaster
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location,
  upper_limit = 1.0,
  steps,
  plot_animation = F
)

# set.seed(s)
burn_result_cpp <- simulate_fire_deterministic(
  vegetation = land_cube(landscape)[, , 1],
  terrain = land_cube(landscape)[, , -1], # without vegetation
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location - 1,
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
wind_dir <- 360               # direction, in angles, from which the wind comes
lands_sub <- landscape_base
lands_sub$wdir <- rep(wind_dir * pi / 180, ncell(lands_sub))
lands_sub$wspeed <- 8
cv <- rep(-1, n_veg)
ct <- rep(0, n_b_terrain)
ct[b_wind] <- 0.5

seed <- round(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub,
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  vegetation = land_cube(lands_sub)[, , 1],
  terrain = land_cube(lands_sub)[, , -1],
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location - 1,
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
lands_sub$elev <- rep(seq(800, 1500, length.out = nrow(lands_sub)),
                      each = ncol(lands_sub))
cv <- rep(-5, n_veg)
ct <- rep(0, n_b_terrain)
ct[b_slope] <- -10

seed <- round(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub,
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  vegetation = land_cube(lands_sub)[, , 1],
  terrain = land_cube(lands_sub)[, , -1],
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location - 1,
  upper_limit = 1.0
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))
# the seed works!


#______________________________

# vegetation test
lands_sub <- landscape_base
lands_sub$vegetation <- sample(c(0:(n_veg -1), 99), size = ncell(lands_sub),
                               prob = c(1, 1, 1, 2), replace = TRUE)

cv <- rep(30, n_veg)
ct <- rep(0, n_b_terrain)

seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub,
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location,
  upper_limit = 1.0,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  vegetation = land_cube(lands_sub)[, , 1],
  terrain = land_cube(lands_sub)[, , -1],
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location - 1,
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
cv <- rep(1000, n_veg)
ct <- rep(0, n_b_terrain)
steps <- sample(1:(nrow(lands_sub) / 2), size = 1)
steps <- 0
seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub,
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location,
  upper_limit = 1.0,
  steps,
  plot_animation = F
)

set.seed(seed)
cc <- simulate_fire(
  vegetation = land_cube(lands_sub)[, , 1],
  terrain = land_cube(lands_sub)[, , -1],
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location - 1,
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
  landscape = lands_sub,
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location,
  upper_limit = 1.0,
  steps,
  plot_animation = F
)

# check steps used
cc_used <- simulate_fire_compare(
  vegetation = land_cube(lands_sub)[, , 1],
  terrain = land_cube(lands_sub)[, , -1],
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location - 1,
  upper_limit = 1.0,
  steps
)
cc_used$steps_used # Ok, in 17 steps burns everything.


# Notes on spread simulation ----------------------------------------------

# Simulations differ between c++ and R when extreme probability values are
# generated, because the precision of floats makes a difference.


# Similarity functions checks ---------------------------------------------


lands_sub$vegetation <- sample(c(0:(n_veg -1), 99), size = ncell(lands_sub),
                               prob = c(1, 1, 1, 2), replace = TRUE)

cv <- rep(-1.5, n_veg)
ct <- rnorm(n_b_terrain)

# simulate fires
set.seed(1)
fire_1 <- simulate_fire_compare_veg(
  vegetation = land_cube(lands_sub)[, , 1],
  terrain = land_cube(lands_sub)[, , -1],
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location - 1,
  upper_limit = 1.0
)

set.seed(1)
fire_1_ <- simulate_fire_compare_veg(
  vegetation = land_cube(lands_sub)[, , 1],
  terrain = land_cube(lands_sub)[, , -1],
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location - 1,
  upper_limit = 1.0
)

set.seed(2)
fire_2 <- simulate_fire_compare_veg(
  vegetation = land_cube(lands_sub)[, , 1],
  terrain = land_cube(lands_sub)[, , -1],
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location - 1,
  upper_limit = 1.0
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
cv <- rnorm(n_veg, -1)
ct <- rnorm(n_b_terrain)

# simulate fires
fire_1 <- simulate_fire_compare(
  vegetation = land_cube(lands_sub)[, , 1],
  terrain = land_cube(lands_sub)[, , -1],
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location - 1,
  upper_limit = 1.0
)

fire_2 <- simulate_fire_compare(
  vegetation = land_cube(lands_sub)[, , 1],
  terrain = land_cube(lands_sub)[, , -1],
  coef_veg = cv,
  coef_terrain = ct,
  ignition_cells = ig_location - 1,
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