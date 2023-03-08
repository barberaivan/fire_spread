# Tests for the spread functions, written in both cpp and R

library(Rcpp)
library(terra)
library(tidyverse)

# load cpp and R functions:
sourceCpp("spread_functions.cpp")
source("spread_functions.R")

# spread_around test ------------------------------------------------------


# create data for testing
# model coefficients (includes intercept)
coefs <- c(0.5,
           -1.5, -1.3, -0.5,
           1, 1, 1, 1, 1)
names(coefs) <- c("intercept",
                  "subalpine", "wet", "dry",
                  "fwi",
                  "aspect",
                  "wind",
                  "elev",
                  "slope")
### IMPORTANT  ---> wind and elevation parameters must be the last ones.

elev_column <- which(names(coefs) == "elev") - 1
wind_column <- which(names(coefs) == "wind") - 1 # -1 because the design matrix will have no intercept

# landscape raster
size <- 30
n_rows <- size
n_cols <- size

landscape <- rast(
  ncol = n_cols, nrow = n_rows,
  nlyrs = length(coefs) - 2, # intercept and slope absent
  xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000,
  names = names(coefs)[-c(1, length(coefs))]
)

# fill vegetation
veg_vals <-  rmultinom(ncell(landscape), size = 1, prob = rep(0.25, 4))[1:3, ] %>% t
values(landscape)[, 1:3] <- veg_vals
landscape$fwi <- rnorm(ncell(landscape))            # fwi anomalies
landscape$aspect <- cos(runif(ncell(landscape), 0, 2 * pi) - 315 * pi / 180)  # northwestyness
landscape$wind <- runif(ncell(landscape), 0, 2 * pi)    # wind direction in radians
landscape$elev <- runif(ncell(landscape), 0, 2200)  # m above sea level

# vector of distances between a cell and its neighbours
# (used to compute slope effect)
distances <- rep(res(landscape)[1], 8) # sides
distances[c(1, 3, 6, 8)] <- res(landscape)[1] * sqrt(2)

# turn landscape into matrix
predmat <- values(landscape)


# spread_onepix test ------------------------------------------------------

# start
burning_cell <- sample(1:ncell(landscape), size = 1)
neighs_raw <- adjacent(landscape, 1, "queen") %>% as.numeric
not_na <- !is.na(neighs_raw)
pos <- (1:8)[not_na]
neighs_id <- neighs_raw[not_na]

s <- 35 # define seed
set.seed(s) # set seed to compare the random simulation for burning
spread_result_r <- spread_onepix_r(
  data_burning = predmat[burning_cell, , drop = T], # USE DROP!!
  data_neighbour = predmat[neighs_id[1], , drop = T],
  coef = coefs,
  position = pos[1],
  distances = distances,
  elev_column = elev_column,
  wind_column = wind_column,
  upper_limit = 1
)

set.seed(s)
spread_result_cpp_burn <- spread_onepix_cpp(
  data_burning = predmat[burning_cell, , drop = F], # USE DROP!!
  data_neighbour = predmat[neighs_id[1], , drop = F],
  coef = coefs,
  position = pos[1] - 1, # -1 for 0-indexing
  distance = distances[pos[1]],
  upper_limit = 1,
  wind_column = which(names(landscape) == "wind") - 1,
  elev_column = which(names(landscape) == "elev") - 1  # -1 for 0-indexing
)

spread_result_cpp_prob <- spread_onepix_prob_cpp(
  data_burning = predmat[burning_cell, , drop = F], # USE DROP!!
  data_neighbour = predmat[neighs_id[1], , drop = F],
  coef = coefs,
  position = pos[1] - 1, # -1 for 0-indexing
  distance = distances[pos[1]],
  upper_limit = 1,
  wind_column = which(names(landscape) == "wind") - 1,
  elev_column = which(names(landscape) == "elev") - 1  # -1 for 0-indexing
)


spread_result_r["burn"] == spread_result_cpp_burn
spread_result_r["probs"] == spread_result_cpp_prob
# OK



# simulate_fire -------------------------------------------------------

# version of simulate_fire that uses a matrix representation of the landscape
# to avoid extra computations when obtaining the neighbours.

# create data for testing
# model coefficients (includes intercept)
coefs <- c(0.5,
           -1.5, -1.3, -0.5,
           1, 1, 1, 1, 1)
names(coefs) <- c("intercept",
                  "subalpine", "wet", "dry",
                  "fwi",
                  "aspect",
                  "wind",
                  "elev",
                  "slope")
### IMPORTANT  ---> wind, elevation and slope parameters must be the last ones.

elev_column <- which(names(coefs) == "elev") - 1
wind_column <- which(names(coefs) == "wind") - 1 # -1 because the design matrix will have no intercept

# landscape raster
size <- 10
n_rows <- size
n_cols <- size

landscape <- rast(
  ncol = n_cols, nrow = n_rows,
  nlyrs = length(coefs) - 2, # intercept and slope absent
  xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000,
  names = names(coefs)[-c(1, length(coefs))]
)

# fill vegetation
veg_vals <-  rmultinom(ncell(landscape), size = 1, prob = rep(0.25, 4))[1:3, ] %>% t
values(landscape)[, 1:3] <- veg_vals
landscape$fwi <- rnorm(ncell(landscape))            # fwi anomalies
landscape$aspect <- cos(runif(ncell(landscape), 0, 2 * pi) - 315 * pi / 180)  # northwestyness
landscape$wind <- runif(ncell(landscape), 0, 2 * pi)    # wind direction in radians
landscape$elev <- runif(ncell(landscape), 0, 2200)  # m above sea level

# make array landscape for cpp function
landscape_arr <- array(NA,
                       dim = c(nrow(landscape), ncol(landscape), nlyr(landscape)))
landscape_values <- terra::values(landscape) # get values in matrix form
for(l in 1:nlyr(landscape)) {
  landscape_arr[, , l] <- matrix(landscape_values[, l],
                                 nrow(landscape), ncol(landscape),
                                 byrow = TRUE) # byrow because terra provides
  # the values this way.
}

# vector of distances between a cell and its neighbours
# (used to compute slope effect)
distances <- rep(res(landscape)[1], 8) # sides
distances[c(1, 3, 6, 8)] <- res(landscape)[1] * sqrt(2)

# make ignition point(s)
ig_location <- matrix(data = c(round(nrow(landscape) / 2), round(ncol(landscape) / 2)),
                      ncol = 1)

ig_cell <- sample(1:ncell(landscape), 1)
ig_location <- rowColFromCell(landscape, ig_cell) %>% t

set.seed(s)
burn_result_r <- simulate_fire_deterministic_r(
  landscape = landscape,
  burnable = matrix(1, nrow(landscape), ncol(landscape)),
  ignition_cells = ig_location,
  coef = coefs,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
)

set.seed(s)
burn_result_cpp <- simulate_fire_deterministic_cpp(
  landscape = landscape_arr,
  burnable = matrix(1, nrow(landscape), ncol(landscape)),
  ignition_cells = ig_location - 1,
  coef = coefs,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
)

landscape_cpp <- landscape[[1]]
landscape_r <- landscape[[1]]

cc <- burn_result_cpp
rr <- burn_result_r

for(i in 1:ncol(cc)) {
  cc[, i] <- burn_result_cpp[i, ]
  rr[, i] <- burn_result_r[i, ]
}
values(landscape_cpp) <- as.numeric(cc)
values(landscape_r) <- as.numeric(rr)

par(mfrow = c(1, 2))
plot(landscape_cpp, col = c("green", "black"), main = "C++")
plot(landscape_r, col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

# Perfect.

# luego testear efectos de paisajes particulares y ver si las funciones
# estocásticas dan patrones similares.

# Test landscape effects (matrix representation) --------------------------

# function to turn landscape from SpatRaster to array with matrix representation
land_cube <- function(x) {
  v <- values(x)
  a <- array(NA, dim = c(nrow(x), ncol(x), nlyr(x)),
             dimnames = list(row = NULL, col = NULL, layer = names(x)))
  for(l in 1:nlyr(x)) a[, , l] <- matrix(v[, l], nrow(x), ncol(x), byrow = TRUE)
  return(a)
}

# Function to turn burned matrix into SpatRaster (for plotting)
rast_from_mat <- function(m, fill_raster) { # fill_raster is a SpatRaster from terra
  mt <- t(m)
  for(i in 1:nrow(m)) mt[, i] <- m[i, ]
  r <- fill_raster[[1]]
  values(r) <- as.numeric(mt)

  return(r)
}

par(mfrow = c(1, 2))
plot(landscape_cpp, col = c("green", "black"), main = "C++")
plot(landscape_r, col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

# create data for testing
# model coefficients (includes intercept)
coefs <- c(0.5,
           -1.5, -1.3, -0.5,
           1, 1, 1, 1, 1)
names(coefs) <- c("intercept",
                  "subalpine", "wet", "dry",
                  "fwi",
                  "aspect",
                  "wind",
                  "elev",
                  "slope")
### IMPORTANT  ---> wind and elevation parameters must be the last ones.

elev_column <- which(names(coefs) == "elev") - 1
wind_column <- which(names(coefs) == "wind") - 1 # -1 because the design matrix will have no intercept

# landscape raster
size <- 30
n_rows <- size
n_cols <- size

landscape <- rast(
  ncol = n_cols, nrow = n_rows,
  nlyrs = length(coefs) - 2, # intercept and slope absent
  xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000,
  names = names(coefs)[-c(1, length(coefs))]
)

ig_location <- matrix(rep(round(size / 2), 2), 2, 1)

# .......................................................................

# simulate landscapes where one layer varies and check that results are OK.

# vegetation tests #

landscape_base <- landscape
landscape_base$subalpine <- 0
landscape_base$wet <- 0
landscape_base$dry <- 0 # it's all shrubland
landscape_base$fwi <- 0
landscape_base$aspect <- 0
landscape_base$wind <- 0 # north wind
landscape_base$elev <- elevation_mean

# subalpine
lands_sub <- landscape_base
lands_sub$subalpine[1:(ncell(landscape)/2)] <- 1 # the northern half is subalpine
c_sub <- coefs
c_sub["subalpine"] <- -1.2 # subalpine is not flammable
c_sub["intercept"] <- 0 # shrubland is very flammable
c_sub["wind"] <- 0 # remove wind effect

land_arr <- land_cube(lands_sub)

# Upper half is not flammable, lower is highly flammable
rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location,
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) # good

cc <- simulate_fire_cpp(
  landscape = land_arr, # use the array
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location - 1,
  coef = c_sub,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
) # good

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

# ______________________

# wet
lands_sub <- landscape_base
lands_sub$wet[1:(ncell(landscape)/2)] <- 1 # the northern half is subalpine
c_sub <- coefs
c_sub["wet"] <- 4
c_sub["intercept"] <- -1
c_sub["wind"] <- 0 # remove wind effect

land_arr <- land_cube(lands_sub)

rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location,
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) # good

cc <- simulate_fire_cpp(
  landscape = land_arr, # use the array
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location - 1,
  coef = c_sub,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
) # good

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))


# ______________________

# dry
lands_sub <- landscape_base
lands_sub$dry[1:(ncell(landscape)/2)] <- 1 # the northern half is subalpine
c_sub <- coefs
c_sub["dry"] <- -12 # subalpine is not flammable
c_sub["intercept"] <- 0 # shrubland is very flammable
c_sub["wind"] <- 0 # remove wind effect

land_arr <- land_cube(lands_sub)

rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location,
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) # good

cc <- simulate_fire_cpp(
  landscape = land_arr, # use the array
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location - 1,
  coef = c_sub,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
) # good

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))


#______________________________


# wind test
wind_dir <- 135               # direction, in angles, from which the wind comes
lands_sub <- landscape_base
lands_sub$wind <- rep(wind_dir * pi / 180, ncell(lands_sub))
c_sub <- rep(0, length(coefs))
names(c_sub) <- names(coefs)
c_sub["wind"] <- 3 # increase wind effect
c_sub["intercept"] <- 0

land_arr <- land_cube(lands_sub)

rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location,
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) # good

cc <- simulate_fire_cpp(
  landscape = land_arr, # use the array
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location - 1,
  coef = c_sub,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
) # good

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

#____________________


# slope test
lands_sub <- landscape_base
lands_sub$elev <- rep(seq(2000, 1000, length.out = nrow(lands_sub)),
                      each = ncol(lands_sub)) # fill values by row
c_sub <- rep(0, length(coefs))
names(c_sub) <- names(coefs)
c_sub["slope"] <- 5 # increase slope effect
c_sub["intercept"] <- -1

# elevation increases from below, so the fire should spread upwards because of
# slope effect (elevation effect was removed.)

land_arr <- land_cube(lands_sub)

rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location,
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) # good

cc <- simulate_fire_cpp(
  landscape = land_arr, # use the array
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location - 1,
  coef = c_sub,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
) # good

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

#________________________

# elevation test
lands_sub <- landscape_base
lands_sub$elev <- rep(seq(2000, 1000, length.out = nrow(lands_sub)),
                      each = ncol(lands_sub)) # fill values by row
c_sub <- rep(0, length(coefs))
names(c_sub) <- names(coefs)
c_sub["elev"] <- -3 # high elevation effect (higher is less flammable)
c_sub["intercept"] <- 1

# elevation increases from below, so the fire should spread downwards because of
# the elevation effect, which is negative. Here the slope effect was removed.
land_arr <- land_cube(lands_sub)

rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location,
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) # good

cc <- simulate_fire_cpp(
  landscape = land_arr, # use the array
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location - 1,
  coef = c_sub,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
) # good

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))


# ______________

# aspect test
lands_sub <- landscape_base
lands_sub$aspect <- rep(seq(1, -1, length.out = nrow(lands_sub)),
                        each = ncol(lands_sub)) # fill values by row
c_sub <- rep(0, length(coefs))
names(c_sub) <- names(coefs)
c_sub["aspect"] <- 2 # increase aspect effect
c_sub["intercept"] <- -1

land_arr <- land_cube(lands_sub)

rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location,
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) # good

cc <- simulate_fire_cpp(
  landscape = land_arr, # use the array
  burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
  ignition_cells = ig_location - 1,
  coef = c_sub,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
) # good

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))

# ________________________

# Burnable test
lands_sub <- landscape_base
c_sub <- rep(0, length(coefs))
names(c_sub) <- names(coefs)
c_sub["intercept"] <- 100 # all burns, except non burnable
burnable_sim = matrix(rbinom(ncell(lands_sub), size = 1, prob = 0.5),
                      nrow(land_arr), ncol(land_arr))

land_arr <- land_cube(lands_sub)


rr <- simulate_fire_r(
  landscape = lands_sub, # use the SpatRaster
  burnable = burnable_sim,
  ignition_cells = ig_location,
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) # good

cc <- simulate_fire_cpp(
  landscape = land_arr, # use the array
  burnable = burnable_sim,
  ignition_cells = ig_location - 1,
  coef = c_sub,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
) # good

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++")
plot(rast_from_mat(rr, lands_sub), col = c("green", "black"), main = "R")
par(mfrow = c(1, 1))
# Here R and C++ should return the same
all.equal(cc, rr)

## All tests are OK


# Test simulate_fire_compare ----------------------------------------------

# this function returns a list with objects used to compare fires (simulated
# vs observed.)

lands_sub <- landscape_base
c_sub <- rep(0, length(coefs))
names(c_sub) <- names(coefs)
c_sub["intercept"] <- 100 # all burns, except non burnable

set.seed(round(runif(1, 1, 2000)))
burnable_sim = matrix(rbinom(ncell(lands_sub), size = 1, prob = 0.5),
                      nrow(land_arr), ncol(land_arr))

land_arr <- land_cube(lands_sub)

set.seed(12)
cc_comp <- simulate_fire_compare(
  landscape = land_arr, # use the SpatRaster
  burnable = burnable_sim,
  ignition_cells = ig_location - 1,
  coef = c_sub,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
) # good

set.seed(12)
cc <- simulate_fire_cpp(
  landscape = land_arr, # use the array
  burnable = burnable_sim,
  ignition_cells = ig_location - 1,
  coef = c_sub,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
) # good

par(mfrow = c(1, 2))
plot(rast_from_mat(cc, lands_sub), col = c("green", "black"), main = "C++ simple")
plot(rast_from_mat(cc_comp$burned_layer, lands_sub), col = c("green", "black"), 
     main = "C++ for comparisons")
par(mfrow = c(1, 1))

# Here both C++ functions should be the same.
all.equal(cc, cc_comp$burned_layer)


# compare_fires tests ------------------------------------------------------

# use the c++ function simulate_fire_compare and compute the similarity metrics
# written in c++ and R, to check they give the same result


# make landscape with varying vegetation (random)

# fill vegetation
veg_vals <-  rmultinom(ncell(landscape), size = 1, prob = rep(0.25, 4))[1:3, ] %>% t
values(landscape)[, 1:3] <- veg_vals
landscape$fwi <- rnorm(ncell(landscape))            # fwi anomalies
landscape$aspect <- cos(runif(ncell(landscape), 0, 2 * pi) - 315 * pi / 180)  # northwestyness
landscape$wind <- runif(ncell(landscape), 0, 2 * pi)    # wind direction in radians
landscape$elev <- runif(ncell(landscape), 0, 2200)  # m above sea level

# make array landscape for cpp function
landscape_arr <- land_cube(landscape)

# vector of distances between a cell and its neighbours
# (used to compute slope effect)
distances <- rep(res(landscape)[1], 8) # sides
distances[c(1, 3, 6, 8)] <- res(landscape)[1] * sqrt(2)

# model coefficients
coefs <- c(0.5,
           -1.5, -1.3, -0.5,
           1, 1, 1, 1, 1)
names(coefs) <- c("intercept",
                  "subalpine", "wet", "dry",
                  "fwi",
                  "aspect",
                  "wind",
                  "elev",
                  "slope")

# make ignition point(s)
ig_location <- matrix(data = c(round(nrow(landscape) / 2), round(ncol(landscape) / 2)),
                      ncol = 1)

# simulate fires to compare

fire1 <- simulate_fire_compare(
  landscape = landscape_arr, # use the SpatRaster
  burnable = matrix(1, nrow(landscape_arr), ncol(landscape_arr)),
  ignition_cells = ig_location - 1,
  coef = coefs,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
)

fire2 <- simulate_fire_compare(
  landscape = landscape_arr, # use the SpatRaster
  burnable = matrix(1, nrow(landscape_arr), ncol(landscape_arr)),
  ignition_cells = ig_location - 1,
  coef = coefs,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
)

comp_rr <- compare_fires_r(fire1, fire2)  
comp_cc <- compare_fires_try(fire1, fire2)    

cbind(comp_rr, comp_cc)

# should be true
all.equal(comp_rr, comp_cc)

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
  landscape = landscape_arr, # use the SpatRaster
  burnable = matrix(1, nrow(landscape_arr), ncol(landscape_arr)),
  ignition_cells = ig1 - 1,
  coef = coefs,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 0
)

ig2 <- matrix(data = c(round(nrow(landscape) / 2), ncol(landscape)), ncol = 1)
fire2 <- simulate_fire_compare(
  landscape = landscape_arr, # use the SpatRaster
  burnable = matrix(1, nrow(landscape_arr), ncol(landscape_arr)),
  ignition_cells = ig2 - 1,
  coef = coefs,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 0
)

comp_rr <- compare_fires_r(fire1, fire2)  
comp_cc <- compare_fires_try(fire1, fire2)  
# should be true:
all.equal(comp_rr, comp_cc)
(c(comp_rr["overlap_sp"], comp_cc["overlap_sp"]) %>% unique) == 0 



# emulate_loglik test ------------------------------------------------------


# fill vegetation
veg_vals <-  rmultinom(ncell(landscape), size = 1, prob = rep(0.25, 4))[1:3, ] %>% t
values(landscape)[, 1:3] <- veg_vals
landscape$fwi <- rnorm(ncell(landscape))            # fwi anomalies
landscape$aspect <- cos(runif(ncell(landscape), 0, 2 * pi) - 315 * pi / 180)  # northwestyness
landscape$wind <- runif(ncell(landscape), 0, 2 * pi)    # wind direction in radians
landscape$elev <- runif(ncell(landscape), 0, 2200)  # m above sea level

# make array landscape for cpp function
landscape_arr <- land_cube(landscape)

# vector of distances between a cell and its neighbours
# (used to compute slope effect)
distances <- rep(res(landscape)[1], 8) # sides
distances[c(1, 3, 6, 8)] <- res(landscape)[1] * sqrt(2)

# model coefficients
coefs <- c(0.5,
           -1.5, -1.3, -0.5,
           1, 1, 1, 1, 1)
names(coefs) <- c("intercept",
                  "subalpine", "wet", "dry",
                  "fwi",
                  "aspect",
                  "wind",
                  "elev",
                  "slope")

# make ignition point(s)
ig_location <- matrix(data = c(round(nrow(landscape) / 2), round(ncol(landscape) / 2)),
                      ncol = 1)


# make reference fire (the observed one)
fire1 <- simulate_fire_compare(
  landscape = landscape_arr, # use the SpatRaster
  burnable = matrix(1, nrow(landscape_arr), ncol(landscape_arr)),
  ignition_cells = ig_location - 1,
  coef = coefs,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
)

loglik_fire1 <- emulate_loglik_try(
  landscape = landscape_arr, # use the SpatRaster
  burnable = matrix(1, nrow(landscape_arr), ncol(landscape_arr)),
  ignition_cells = ig_location - 1,
  coef = coefs,
  wind_layer = wind_column - 1,
  elev_layer = elev_column - 1,
  distances = distances,
  upper_limit = 1.0,

  fire_ref = fire1,
  n_replicates = 10
)

loglik_fire1


# TO DO: check metrics with carefully desined landscapes.