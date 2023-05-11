# Tests for the spread functions, written in both cpp and R for double-checking.
# These tests are written to perform deep and interactive inspection of the 
# functions, necessary when the code has gone through large changes. 
# The tests.R code is a fast version.

# Packages and codes ------------------------------------------------------

library(Rcpp)
library(terra)
library(tidyverse)

# load cpp and R functions:
sourceCpp("spread_functions.cpp")
source("spread_functions.R")
sourceCpp("similarity_functions.cpp")
source("similarity_functions.R")

# constants ---------------------------------------------------------------

n_veg_types <- 6
n_terrain <- 4

# fire spread parameters (coef)
coefs <- c(rnorm(n_veg_types), 
           rnorm(n_terrain))

# terrain names
terrain_names <- c("northing", "elev", "windir")
terrain_params <- c(terrain_names, "slope")

# landscape raster
size <- 30
n_rows <- size
n_cols <- size
res <- 30
n_cells <- n_rows * n_cols

landscape <- rast(
  ncol = n_cols, nrow = n_rows, res = res, crs = "EPSG:5343",
  nlyrs = 1 + n_terrain - 1, # veg in just one layer, and slope absent
  xmin = 0, xmax = res * n_cols, ymin = 0, ymax = res * n_rows,
  names = c("vegetation", terrain_names)
)

# fill data
veg_vals <-  sample(0:(n_veg_types - 1), size = n_cells, replace = TRUE)
values(landscape)[, 1] <- veg_vals
landscape$northing <- cos(runif(ncell(landscape), 0, 2 * pi) - 315 * pi / 180)
landscape$elev <- runif(ncell(landscape), 0, 2200) / elevation_sd # scaled
landscape$windir <- runif(ncell(landscape), 0, 2 * pi) # radians

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
spread_result_r <- spread_onepix_r(
  vegetation_type = predmat[burning_cell, "vegetation"], # pass a 0-starting index
  terrain_burning = predmat[burning_cell, terrain_names, drop = T], # USE DROP!!
  terrain_neighbour = predmat[neighs_id[1], terrain_names, drop = T],
  coef_veg = coefs[1:n_veg_types],
  coef_terrain = coefs[(n_veg_types + 1) : length(coefs)],
  position = pos[1],
  upper_limit = 1
)
spread_result_r

set.seed(s)
spread_result_cpp_burn <- spread_onepix_cpp(
  vegetation_type = predmat[burning_cell, "vegetation"], # pass a 0-starting index
  terrain_burning = predmat[burning_cell, terrain_names, drop = T], # USE DROP!!
  terrain_neighbour = predmat[neighs_id[1], terrain_names, drop = T],
  coef_veg = coefs[1:n_veg_types],
  coef_terrain = coefs[(n_veg_types + 1) : length(coefs)],
  position = pos[1] - 1,
  upper_limit = 1
)

spread_result_cpp_prob <- spread_onepix_prob_cpp(
  vegetation_type = predmat[burning_cell, "vegetation"], # pass a 0-starting index
  terrain_burning = predmat[burning_cell, terrain_names, drop = T], # USE DROP!!
  terrain_neighbour = predmat[neighs_id[1], terrain_names, drop = T],
  coef_veg = coefs[1:n_veg_types],
  coef_terrain = coefs[(n_veg_types + 1) : length(coefs)],
  position = pos[1] - 1,
  upper_limit = 1
)

spread_result_r["burn"]; spread_result_cpp_burn
spread_result_r["probs"]; spread_result_cpp_prob
# OK


# simulate_fire -------------------------------------------------------

# version of simulate_fire that uses a matrix representation of the landscape
# to avoid extra computations when obtaining the neighbours.

# landscape raster
size <- 30
n_rows <- size
n_cols <- size
res <- 30
n_cells <- n_rows * n_cols

landscape <- rast(
  ncol = n_cols, nrow = n_rows, res = res, crs = "EPSG:5343",
  nlyrs = 1 + n_terrain - 1, # veg in just one layer, and slope absent
  xmin = 0, xmax = res * n_cols, ymin = 0, ymax = res * n_rows,
  names = c("vegetation", terrain_names)
)

# fill data
veg_vals <-  sample(0:(n_veg_types - 1), size = n_cells, replace = TRUE)
values(landscape)[, 1] <- veg_vals
landscape$northing <- cos(runif(ncell(landscape), 0, 2 * pi) - 315 * pi / 180)
landscape$elev <- runif(ncell(landscape), 0, 2200) / elevation_sd # scaled
landscape$windir <- runif(ncell(landscape), 0, 2 * pi) # radians

# make array landscape for cpp function
terrain_arr <- land_cube(landscape)
vegetation <- matrix(values(landscape$vegetation), n_rows, n_cols, byrow = TRUE)

# make ignition point(s)
ig_cell <- sample(1:ncell(landscape), 1)
ig_location <- rowColFromCell(landscape, ig_cell) %>% t

# sample parameters
coefs <- rnorm(n_veg_types + n_terrain, 0, 10)

# set.seed(s)
burn_result_r <- simulate_fire_deterministic_r(
  vegetation,
  terrain = landscape[[terrain_names]],
  ignition_cells = ig_location,
  coef = coefs,
  n_veg_types = n_veg_types,
  upper_limit = 1.0,
  plot_animation = TRUE
)

# set.seed(s)
burn_result_cpp <- simulate_fire_deterministic_cpp(
  vegetation,
  terrain_arr,
  ignition_cells = ig_location - 1,
  coef = coefs,
  n_veg_types = n_veg_types,
  upper_limit = 1.0
)

res_cpp <- rast_from_mat(burn_result_cpp, landscape)
res_r <- rast_from_mat(burn_result_r, landscape)

par(mfrow = c(1, 2))
plot(res_cpp, col = c("green", "black"), main = "C++")
plot(res_r, col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

# Perfect.


# Test landscape effects (matrix representation) --------------------------

# landscape raster
size <- 30
n_rows <- size
n_cols <- size
res <- 30
n_cells <- n_rows * n_cols

landscape <- rast(
  ncol = n_cols, nrow = n_rows, res = res, crs = "EPSG:5343",
  nlyrs = 1 + n_terrain - 1, # veg in just one layer, and slope absent
  xmin = 0, xmax = res * n_cols, ymin = 0, ymax = res * n_rows,
  names = c("vegetation", terrain_names)
)

# fill data
veg_vals <-  sample(0:(n_veg_types - 1), size = n_cells, replace = TRUE)
values(landscape)[, 1] <- veg_vals
landscape$northing <- cos(runif(ncell(landscape), 0, 2 * pi) - 315 * pi / 180)
landscape$elev <- runif(ncell(landscape), 0, 2200) / elevation_sd # scaled
landscape$windir <- runif(ncell(landscape), 0, 2 * pi) # radians

# make array landscape for cpp function
terrain_arr <- land_cube(landscape)
vegetation <- matrix(values(landscape$vegetation), n_rows, n_cols, byrow = TRUE)

# ignition
ig_location <- matrix(rep(round(size / 2), 2), 2, 1)

# .......................................................................

# simulate landscapes where one layer varies and check that results are OK.

# vegetation tests #
# lower portion is non-burnable.

landscape_base <- landscape
values(landscape_base) <- 0
coef_terrain <- rnorm(n_terrain)

for(v in 1:n_veg_types) {
  # v = 1
  vegetation[, ] <- 99
  vegetation[1:(1 + round((n_rows / 2))), ] <- v-1 # veg code starts at zero
  
  coef_veg <- rep(-30, n_veg_types)
  coef_veg[v] <- 0
  
  # Upper half is not flammable, lower is highly flammable
  rr <- simulate_fire_r(
    terrain = landscape_base, # use the SpatRaster
    vegetation = vegetation,
    ignition_cells = ig_location,
    coef = c(coef_veg, coef_terrain),
    n_veg_types = n_veg_types,
    upper_limit = 1.0
  )
  
  cc <- simulate_fire_cpp(
    terrain = land_cube(landscape_base), # use the SpatRaster
    vegetation = vegetation,
    ignition_cells = ig_location - 1,
    coef = c(coef_veg, coef_terrain),
    n_veg_types = n_veg_types,
    upper_limit = 1.0
  )
  
  par(mfrow = c(1, 2))
  plot(rast_from_mat(cc, landscape), col = c("green", "black"), 
       main = paste("C++,", "veg", v))
  plot(rast_from_mat(rr, landscape), col = c("green", "black"), 
       main = paste("R,", "veg", v))
  par(mfrow = c(1, 1))
}


#______________________________


# wind test
wind_dir <- 270               # direction, in angles, from which the wind comes
lands_sub <- landscape_base
lands_sub$windir <- rep(wind_dir * pi / 180, ncell(lands_sub))
coef_veg <- rep(-1, n_veg_types)
coef_terrain <- rep(0, n_terrain)
coef_terrain[terrain_params == "windir"] <- 3
vegetation[, ] <- 0

rr <- simulate_fire_r(
  terrain = lands_sub[[terrain_names]], # use the SpatRaster
  vegetation = vegetation,
  ignition_cells = ig_location,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 1.0
)

cc <- simulate_fire_cpp(
  terrain = land_cube(lands_sub), 
  vegetation = vegetation,
  ignition_cells = ig_location - 1,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 1.0
) 

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

#____________________


# slope and elevation tests

# elevation increases from below, so the fire should spread upwards because of
# slope effect (elevation effect was removed.)

lands_sub <- landscape_base
lands_sub$elev <- rep(seq(2, -2, length.out = nrow(lands_sub)),
                      each = ncol(lands_sub)) # fill values by row
# elevation values in the standardized scale!
coef_veg <- rep(-1, n_veg_types)
coef_terrain <- rep(0, n_terrain)
# coef_terrain[terrain_params == "slope"] <- 2 # to test slope
 coef_terrain[terrain_params == "elev"] <- -0.5 # to test elevation
vegetation[, ] <- 0

rr <- simulate_fire_r(
  terrain = lands_sub[[terrain_names]], # use the SpatRaster
  vegetation = vegetation,
  ignition_cells = ig_location,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 1.0
)

cc <- simulate_fire_cpp(
  terrain = land_cube(lands_sub), 
  vegetation = vegetation,
  ignition_cells = ig_location - 1,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 1.0
) 

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

# _________________________

# Northing test
lands_sub <- landscape_base
lands_sub$northing <- rep(seq(1, -1, length.out = nrow(lands_sub)),
                          each = ncol(lands_sub)) # fill values by row
coef_veg <- rep(-2, n_veg_types)
coef_terrain <- rep(0, n_terrain)
coef_terrain[terrain_params == "northing"] <- -3 # to test elevation
vegetation[, ] <- 0

rr <- simulate_fire_r(
  terrain = lands_sub[[terrain_names]], # use the SpatRaster
  vegetation = vegetation,
  ignition_cells = ig_location,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 1.0
)

cc <- simulate_fire_cpp(
  terrain = land_cube(lands_sub), 
  vegetation = vegetation,
  ignition_cells = ig_location - 1,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 1.0
) 

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

## All tests are OK


# Test simulate_fire_compare ----------------------------------------------

# this function returns a list with objects used to compare fires (simulated
# vs observed.)

lands_sub <- landscape_base
coef_veg <- rep(100, n_veg_types) # all burnable burns
coef_terrain <- rep(0, n_terrain)
vegetation <- matrix(
  sample(c(0:(n_veg_types - 1), 99), size = n_cells, replace = T),
  n_rows, n_cols
)

set.seed(12)
cc_comp <- simulate_fire_compare(
  terrain = land_cube(lands_sub), 
  vegetation = vegetation,
  ignition_cells = ig_location - 1,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 1.0
) # good

set.seed(12)
cc <- simulate_fire_cpp(
  terrain = land_cube(lands_sub), 
  vegetation = vegetation,
  ignition_cells = ig_location - 1,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 1.0
) # good

par(mfrow = c(2, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++ simple")
plot(rast_from_mat(cc_comp$burned_layer, lands_sub), col = c("green", "black"),
     main = "C++ for comparisons")
plot(rast_from_mat(vegetation, lands_sub), main = "vegetation")
par(mfrow = c(1, 1))

# Here both C++ functions should be the same.
all.equal(cc, cc_comp$burned_layer)


# compare_fires tests ------------------------------------------------------

# use the c++ function simulate_fire_compare and compute the similarity metrics
# written in c++ and R, to check they give the same result

# make landscape with varying vegetation (random)
lands_sub <- landscape_base
coef_veg <- rnorm(n_veg_types, -1) 
coef_terrain <- rnorm(n_terrain)
vegetation <- matrix(
  sample(c(0:(n_veg_types - 1), 99), size = n_cells, replace = T),
  n_rows, n_cols
)
# make ignition point(s)
ig_location <- matrix(data = c(round(nrow(landscape) / 2), round(ncol(landscape) / 2)),
                      ncol = 1)

# simulate fires to compare (c++)
fire1 <- simulate_fire_compare(
  terrain = land_cube(lands_sub), 
  vegetation = vegetation,
  ignition_cells = ig_location - 1,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 1.0
)

fire2 <- simulate_fire_compare(
  terrain = land_cube(lands_sub), 
  vegetation = vegetation,
  ignition_cells = ig_location - 1,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 1.0
)

comp_rr <- compare_fires_r(fire1, fire2)
comp_cc <- compare_fires_try(fire1, fire2)

cbind(comp_rr, comp_cc)
all(abs(comp_rr - comp_cc) < 1e-7)

# test the function with equal fires (all indices should be 1)
comp_rr <- compare_fires_r(fire1, fire1)
comp_cc <- compare_fires_try(fire1, fire1)
# should be true:
all.equal(comp_rr, comp_cc)
(c(comp_rr, comp_cc) %>% unique) == 1

# test the function with completely disjoint fires
# (overlap_sp should be 0)
ig1 <- matrix(data = c(round(nrow(landscape) / 2), 1), ncol = 1)
fire1 <- simulate_fire_compare(
  terrain = land_cube(lands_sub), 
  vegetation = vegetation,
  ignition_cells = ig1 - 1,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 0
)

ig2 <- matrix(data = c(round(nrow(landscape) / 2), ncol(landscape)), ncol = 1)
fire2 <- simulate_fire_compare(
  terrain = land_cube(lands_sub), 
  vegetation = vegetation,
  ignition_cells = ig2 - 1,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 0
)

comp_rr <- compare_fires_r(fire1, fire2)
comp_cc <- compare_fires_try(fire1, fire2)
# should be true:
all(abs(comp_rr - comp_cc) < 1e-7)
(c(comp_rr["overlap_sp"], comp_cc["overlap_sp"]) %>% unique) == 0



# emulate_loglik test ------------------------------------------------------


# landscape
lands_sub <- landscape_base
lands_sub$northing <- cos(runif(ncell(landscape), 0, 2 * pi) - 315 * pi / 180)
lands_sub$elev <- (runif(ncell(landscape), 0, 2200) - elevation_mean) / elevation_sd # scaled
lands_sub$windir <- runif(ncell(landscape), 0, 2 * pi) # radians

coef_veg <- rnorm(n_veg_types, -0.5) 
coef_terrain <- rnorm(n_terrain)
vegetation <- matrix(
  sample(c(0:(n_veg_types - 1), 99), size = n_cells, replace = T),
  n_rows, n_cols
)

# make ignition point(s)
ig_location <- matrix(data = c(round(nrow(landscape) / 2), round(ncol(landscape) / 2)),
                      ncol = 1)

# make reference fire (the observed one)
fire1 <- simulate_fire_compare(
  terrain = land_cube(lands_sub), 
  vegetation = vegetation,
  ignition_cells = ig_location - 1,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 1.0
)

loglik_fire1 <- emulate_loglik_try(
  terrain = land_cube(lands_sub), 
  vegetation = vegetation,
  ignition_cells = ig_location - 1,
  coef = c(coef_veg, coef_terrain),
  n_veg_types = n_veg_types,
  upper_limit = 1.0,

  fire_ref = fire1,
  n_replicates = 1000
)
summary(loglik_fire1)