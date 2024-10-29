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

n_veg <- 5
n_nd <- 2
n_terrain <- 3 # it doesn't include steps
n_b_terrain <- n_terrain - 1
n_layers <- 1 + n_nd + n_terrain # extra 1 is the intercept
layer_names <- c("vegetation", "vfi", "tfi", "elev", "wdir", "wspeed")

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
landscape$vfi <- rep(0, ncell(landscape))
landscape$tfi <- rep(0, ncell(landscape))
landscape$elev <- rnorm(ncell(landscape), 1500, 300)
landscape$wdir <- runif(ncell(landscape), 0, 2 * pi) # radians
landscape$wspeed <- abs(rnorm(ncell(landscape), 0, 2))

# turn landscape into matrix
predmat <- values(landscape) # first column is vegetation
predmat[, 1] <- 1

landscape_base <- landscape
# spread_onepix test ------------------------------------------------------

# start
coef_intercepts <- rnorm(n_veg)
coef_nd <- rnorm(n_nd, 0 , 10)         # non-directional
coef_terrain <- rnorm(n_b_terrain, 0, 3)         # terrain

burning_cell <- sample(1:ncell(landscape), size = 1)
neighs_raw <- adjacent(landscape, 1, "queen") %>% as.numeric
not_na <- !is.na(neighs_raw)
pos <- (1:8)[not_na]
neighs_id <- neighs_raw[not_na]

s <- round(runif(1, 100, 1000)) # define seed
set.seed(s) # set seed to compare the random simulation for burning
spread_result_r <- spread_one_cell_r(
  vegetation = predmat[neighs_id[1], 1],
  data_nd = predmat[neighs_id[1], 2:(n_nd+1), drop = F],
  data_terrain_source = predmat[burning_cell, -(1:3), drop = T], # USE DROP!!
  data_terrain_target = predmat[neighs_id[1], -(1:3), drop = T],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  position = pos[1],
  upper_limit = 1
)
spread_result_r

set.seed(s)
spread_result_cpp_burn <- spread_one_cell(
  vegetation = predmat[neighs_id[1], 1],
  data_nd = predmat[neighs_id[1], 2:(n_nd+1), drop = F],
  data_terrain_source = predmat[burning_cell, -(1:3), drop = T], # USE DROP!!
  data_terrain_target = predmat[neighs_id[1], -(1:3), drop = T],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  position = pos[1] - 1,
  upper_limit = 1
)

spread_result_cpp_prob <- spread_one_cell_prob(
  vegetation = predmat[neighs_id[1], 1] - 1,
  data_nd = predmat[neighs_id[1], 2:(n_nd+1), drop = F],
  data_terrain_source = predmat[burning_cell, -(1:3), drop = T], # USE DROP!!
  data_terrain_target = predmat[neighs_id[1], -(1:3), drop = T],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  position = pos[1] - 1,
  upper_limit = 1
)

spread_result_r["burn"]; spread_result_cpp_burn
spread_result_r["probs"]; spread_result_cpp_prob
# OK, the seed seems to work.


# simulate_fire test ------------------------------------------------------

# sample parameters
s <- round(runif(1, 100, 1000)) # define seed
coef_intercepts <- rnorm(n_veg, -1)
coef_nd <- rnorm(n_nd, 0 , 5)         # non-directional
coef_terrain <- rnorm(n_b_terrain, 0, 3)         # terrain

# make ignition point(s)
ig_cell <- sample(1:ncell(landscape), 1)
ig_location <- rowColFromCell(landscape, ig_cell) %>% t
steps <- 7 # check this

set.seed(s)
burn_result_r <- simulate_fire_r(
  landscape = landscape, # SpatRaster
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location,
  upper_limit = 1.0,
  steps,
  plot_animation = F,
  det = T
)

set.seed(s)
land <- land_cube(landscape)
burn_result_cpp <- simulate_fire_deterministic(
  layer_vegetation = land[, , 1],
  layer_nd = land[, , 2:3],
  layer_terrain = land[, , 4:6],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location-1,
  upper_limit = 1.0,
  steps
)

res_cpp <- rast_from_mat(burn_result_cpp, landscape)
res_r <- rast_from_mat(burn_result_r, landscape)

par(mfrow = c(1, 2))
plot(res_cpp, col = c("green", "black"), main = "C++")
plot(res_r, col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

all.equal(values(res_cpp), values(res_r))
# Perfect.

# Test predictors effects ------------------------------------------------

# simulate landscapes where one layer varies and check that results are OK.

# ignition
ig_location <- matrix(rep(round(size / 2), 2), 2, 1)
landscape_base <- landscape
values(landscape_base) <- 0

#______________________________

# wind test
wind_dir <- 180               # direction, in angles, from which the wind comes
lands_sub <- landscape_base
lands_sub$wdir <- rep(wind_dir * pi / 180, ncell(lands_sub))
lands_sub$wspeed <- 20
land <- land_cube(lands_sub)

coef_intercepts <- rep(-1, n_veg)
coef_nd <- rep(0, n_nd)
coef_terrain <- rep(0, n_b_terrain)         # terrain
coef_terrain[b_wind] <- 5

seed <- round(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # SpatRaster
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location,
  upper_limit = 1.0,
  steps = 0,
  plot_animation = F,
  det = F
)

set.seed(seed)
cc <- simulate_fire(
  layer_vegetation = land[, , 1],
  layer_nd = land[, , 2:3],
  layer_terrain = land[, , 4:6],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location-1,
  upper_limit = 1.0,
  steps = 0
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

all.equal(rr, cc)


#______________________________

# slope test
lands_sub <- landscape_base
lands_sub$elev <- rep(seq(800, 1500, length.out = nrow(lands_sub)),
                      each = ncol(lands_sub))
land <- land_cube(lands_sub)

coef_intercepts <- rep(-2, n_veg)
coef_nd <- rep(0, n_nd)
coef_terrain <- rep(0, n_b_terrain)         # terrain
coef_terrain[b_slope] <- 5


seed <- round(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # SpatRaster
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location,
  upper_limit = 1.0,
  steps = 0,
  plot_animation = F,
  det = F
)

set.seed(seed)
cc <- simulate_fire(
  layer_vegetation = land[, , 1],
  layer_nd = land[, , 2:3],
  layer_terrain = land[, , 4:6],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location-1,
  upper_limit = 1.0,
  steps = 0
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

all.equal(rr, cc)

#______________________________

# vegetation test
lands_sub <- landscape_base
lands_sub$vegetation <- sample(c(0:(n_veg -1), 99), size = ncell(lands_sub),
                               prob = c(rep(1, n_veg), 0), replace = TRUE)
land <- land_cube(lands_sub)

coef_intercepts <- c(-30, 30, 30, -30, -30)
coef_nd <- rep(0, n_nd)
coef_terrain <- rep(0, n_b_terrain)         # terrain


seed <- round(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # SpatRaster
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location,
  upper_limit = 1.0,
  steps = 0,
  plot_animation = F,
  det = F
)

set.seed(seed)
cc <- simulate_fire(
  layer_vegetation = land[, , 1],
  layer_nd = land[, , 2:3],
  layer_terrain = land[, , 4:6],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location-1,
  upper_limit = 1.0,
  steps = 0
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

all.equal(rr, cc)


#______________________________

# vfi/tfi test
lands_sub <- landscape_base
lands_sub$vfi <- rep(seq(-1, 1, length.out = nrow(lands_sub)),
                     each = ncol(lands_sub))
land <- land_cube(lands_sub)

coef_intercepts <- rep(-1, n_veg)
coef_nd <- c(3, 0)
coef_terrain <- rep(0, n_b_terrain)         # terrain

seed <- round(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # SpatRaster
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location,
  upper_limit = 1.0,
  steps = 0,
  plot_animation = F,
  det = F
)

set.seed(seed)
cc <- simulate_fire(
  layer_vegetation = land[, , 1],
  layer_nd = land[, , 2:3],
  layer_terrain = land[, , 4:6],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location-1,
  upper_limit = 1.0,
  steps = 0
)
par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))


#______________________________

# steps test
lands_sub <- landscape_base
lands_sub$vegetation <- 0
land <- land_cube(lands_sub)

coef_intercepts <- rep(-1, n_veg)
coef_nd <- c(0, 0)
coef_terrain <- rep(0, n_b_terrain)         # terrain
steps <- sample(1:(nrow(lands_sub) / 2), size = 1)

seed <- round(runif(1, 10, 1e6))

# steps <- 10
seed <- floor(runif(1, 10, 1e6))

set.seed(seed)
rr <- simulate_fire_r(
  landscape = lands_sub, # SpatRaster
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location,
  upper_limit = 1.0,
  steps = steps,
  plot_animation = F,
  det = F
)

set.seed(seed)
cc <- simulate_fire(
  layer_vegetation = land[, , 1],
  layer_nd = land[, , 2:3],
  layer_terrain = land[, , 4:6],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location-1,
  upper_limit = 1.0,
  steps = steps
)

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"),
     main = paste("C++,", steps, "steps"))
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"),
     main = paste("R,", steps, "steps"))
par(mfrow = c(1, 1))

## extra

# check steps used and counts veg

lands_sub <- landscape_base
lands_sub$vegetation <- 4
land <- land_cube(lands_sub)

cc_used <- simulate_fire_compare_veg(
  layer_vegetation = land[, , 1],
  layer_nd = land[, , 2:3],
  layer_terrain = land[, , 4:6],
  coef_intercepts = rep(1e6, n_veg),
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location-1,
  upper_limit = 1.0,
  steps = 0,
  n_veg = 5
)
cc_used$steps_used # Ok, in 17 steps burns everything.
cc_used$counts_veg

# Notes on spread simulation ----------------------------------------------

# Simulations differ between c++ and R when extreme probability values are
# generated, because the precision of floats makes a difference.


# Similarity functions checks ---------------------------------------------

lands_sub <- landscape_base
lands_sub$vegetation <- sample(c(0:(n_veg -1), 99), size = ncell(lands_sub),
                               prob = c(rep(1, n_veg), 2), replace = TRUE)
land <- land_cube(lands_sub)

coef_intercepts <- rep(1, n_veg)
coef_nd <- rep(0, n_nd)
coef_terrain <- rnorm(n_b_terrain)

# simulate fires
set.seed(1)
fire_1 <- simulate_fire_compare_veg(
  layer_vegetation = land[, , 1],
  layer_nd = land[, , 2:3],
  layer_terrain = land[, , 4:6],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location-1,
  upper_limit = 1.0,
  steps = 0,
  n_veg = 5
)

set.seed(1)
fire_1_ <- simulate_fire_compare_veg(
  layer_vegetation = land[, , 1],
  layer_nd = land[, , 2:3],
  layer_terrain = land[, , 4:6],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location-1,
  upper_limit = 1.0,
  steps = 0,
  n_veg = 5
)

set.seed(2)
fire_2 <- simulate_fire_compare_veg(
  layer_vegetation = land[, , 1],
  layer_nd = land[, , 2:3],
  layer_terrain = land[, , 4:6],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location-1,
  upper_limit = 1.0,
  steps = 0,
  n_veg = 5
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

# and exploration of similarity metrics and their weighting.
lands_sub <- landscape_base
lands_sub$vegetation <- sample(c(0:(n_veg -1), 99), size = ncell(lands_sub),
                               prob = c(rep(1, n_veg), 0.5), replace = TRUE)
land <- land_cube(lands_sub)

# reasonable coefs
coef_intercepts <- rnorm(n_veg, 0)
coef_nd <- rep(0, n_nd)
coef_terrain <- rnorm(n_b_terrain)
steps1 <- sample(5:17, 1)

# simulate fires
fire_1 <- simulate_fire_compare_veg(
  layer_vegetation = land[, , 1],
  layer_nd = land[, , 2:3],
  layer_terrain = land[, , 4:6],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location-1,
  upper_limit = 1.0,
  steps = steps1,
  n_veg = 5
)

fire_2 <- simulate_fire_compare_veg(
  layer_vegetation = land[, , 1],
  layer_nd = land[, , 2:3],
  layer_terrain = land[, , 4:6],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = ig_location-1,
  upper_limit = 1.0,
  steps = 0,
  n_veg = 5
)

common <- fire_1$burned_layer + fire_2$burned_layer

r1 <- rast_from_mat(fire_1$burned_layer, landscape)
r2 <- rast_from_mat(fire_2$burned_layer, landscape)
un <- rast_from_mat(pmax(fire_2$burned_layer, fire_1$burned_layer), landscape)
com <- rast_from_mat(pmin(fire_2$burned_layer, fire_1$burned_layer), landscape)

ov_sp <- overlap_spatial(fire_1, fire_2)

par(mfrow = c(2, 2))
plot(r1, col = c("green", "red"), main = "Fire 1")
plot(r2, col = c("green", "blue"), main = "Fire 2")
plot(com, col = c("green", "purple"), main = "Intersection")
plot(un, col = c("green", "black"), main = "Union")
par(mfrow = c(1, 1))

ov_sp