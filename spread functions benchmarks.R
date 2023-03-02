# benchmarks

library(terra)
library(tidyverse)
library(Rcpp)
library(viridis)
library(microbenchmark)
theme_set(theme_bw())

sourceCpp("spread_functions.cpp")

# Data preparation --------------------------------------------------------

# import landscape (cholila)
land <- readRDS("data_cholila_landscape.rds")

# cholila raster file
land_raster <- rast("data_cholila_elevation.tif")
nrow(land) == ncell(land_raster)

# distances for 30 m resolution
distances <- rep(res(land_raster)[1], 8) # sides
distances[c(1, 3, 6, 8)] <- res(land_raster)[1] * sqrt(2)

# coefficients
coefs <- c(1000,
           0, 0, 0,
           0, 0, 0, 0, 0)

ig_location <- cellFromRowCol(land_raster,
                              nrow(land_raster) / 2, ncol(land_raster) / 2)

ig_location_mat <- rowColFromCell(land_raster, ig_location) %>% t


# landscape array for matrix representation
land_arr <- array(
  NA,
  dim = c(nrow(land_raster), ncol(land_raster), ncol(land)),
  dimnames = list(row = NULL, col = NULL, layer = colnames(land))
)
for(l in 1:ncol(land)) {
  land_arr[, , l] <- matrix(land[, l],
                            nrow(land_raster),
                            ncol(land_raster),
                            byrow = TRUE)
}
str(land_arr)

# Benchmark fool vs cool versions, whole landscape burnable ----------------
# cool version has vector and matrix form

# fool version
start_time_fool <- Sys.time()
fire_fool <- simulate_fire_fool_cpp(
  landscape = land[, 1:7],
  burnable = rep(1, ncell(land_raster)),
  ignition_cells = ig_location - 1,
  n_rowcol = c(nrow(land_raster), ncol(land_raster)),
  coef = coefs,
  wind_column = 6 - 1,
  elev_column = 7 - 1,
  distances = distances,
  upper_limit = 1
)
end_time_fool <- Sys.time()
time_fool <- end_time_fool - start_time_fool


# cool version vector
start_time_cool <- Sys.time()
fire_cool <- simulate_fire_cpp(
  landscape = land[, 1:7],
  burnable = rep(1, ncell(land_raster)),
  ignition_cells = ig_location - 1,
  n_rowcol = c(nrow(land_raster), ncol(land_raster)),
  coef = coefs,
  wind_column = 6 - 1,
  elev_column = 7 - 1,
  distances = distances,
  upper_limit = 1
)
end_time_cool <- Sys.time()
time_cool <- end_time_cool - start_time_cool

# cool version matrix
start_time_cool_m <- Sys.time()
fire_cool_m <- simulate_fire_mat_cpp(
  landscape = land_arr[, , 1:7],
  burnable = matrix(1, nrow(land_raster), ncol(land_raster)),
  ignition_cells = ig_location_mat - 1,
  coef = coefs,
  wind_layer = 6 - 1,
  elev_layer = 7 - 1,
  distances = distances,
  upper_limit = 1
)
end_time_cool_m <- Sys.time()
time_cool_m <- end_time_cool_m - start_time_cool_m

# check all fires were large
sum(fire_cool) / length(fire_cool)
sum(fire_fool) / length(fire_fool)
sum(as.numeric(fire_cool_m)) / length(fire_cool)

# compare
print(time_fool)
print(time_cool)
print(time_cool_m)

# quotient:
(as.numeric(time_cool / 60) / as.numeric(time_fool)) # fool is in min, cool is in sec.
# 0.164394

# cool matrix over cool vector
(as.numeric(time_cool_m) / as.numeric(time_cool))
# 0.2513233


# __________________________________

# with microbenchmark


mbm_all <- microbenchmark(
  fool = simulate_fire_fool_cpp(
    landscape = land[, 1:7],
    burnable = rep(1, ncell(land_raster)),
    ignition_cells = ig_location - 1,
    n_rowcol = c(nrow(land_raster), ncol(land_raster)),
    coef = coefs,
    wind_column = 6 - 1,
    elev_column = 7 - 1,
    distances = distances,
    upper_limit = 1
  ),
  cool_vec = simulate_fire_cpp(
    landscape = land[, 1:7],
    burnable = rep(1, ncell(land_raster)),
    ignition_cells = ig_location - 1,
    n_rowcol = c(nrow(land_raster), ncol(land_raster)),
    coef = coefs,
    wind_column = 6 - 1,
    elev_column = 7 - 1,
    distances = distances,
    upper_limit = 1
  ),
  cool_mat = simulate_fire_mat_cpp(
    landscape = land_arr[, , 1:7],
    burnable = matrix(1, nrow(land_raster), ncol(land_raster)),
    ignition_cells = ig_location_mat - 1,
    coef = coefs,
    wind_layer = 6 - 1,
    elev_layer = 7 - 1,
    distances = distances,
    upper_limit = 1
  ),
  times = 1
)


mbm_all

# Benchmark fool vs cool versions, real burnable layer --------------------

# fool version
start_time_fool <- Sys.time()
fire_fool <- simulate_fire_fool_cpp(
  landscape = land[, 1:7],
  burnable = land[, "burnable"],
  ignition_cells = ig_location - 1,
  n_rowcol = c(nrow(land_raster), ncol(land_raster)),
  coef = coefs,
  wind_column = 6 - 1,
  elev_column = 7 - 1,
  distances = distances,
  upper_limit = 1
)
end_time_fool <- Sys.time()
(time_fool <- end_time_fool - start_time_fool)


# cool version
start_time_cool <- Sys.time()
fire_cool <- simulate_fire_cpp(
  landscape = land[, 1:7],
  burnable = land[, "burnable"],
  ignition_cells = ig_location - 1,
  n_rowcol = c(nrow(land_raster), ncol(land_raster)),
  coef = coefs,
  wind_column = 6 - 1,
  elev_column = 7 - 1,
  distances = distances,
  upper_limit = 1
)
end_time_cool <- Sys.time()
(time_cool <- end_time_cool - start_time_cool)

# cool version matrix
start_time_cool_m <- Sys.time()
fire_cool_m <- simulate_fire_mat_cpp(
  landscape = land_arr[, , 1:7],
  burnable = land_arr[, , "burnable"],
  ignition_cells = ig_location_mat - 1,
  coef = coefs,
  wind_layer = 6 - 1,
  elev_layer = 7 - 1,
  distances = distances,
  upper_limit = 1
)
end_time_cool_m <- Sys.time()
time_cool_m <- end_time_cool_m - start_time_cool_m

# check all fires were large
sum(fire_cool) / sum(land[, "burnable"])
sum(fire_fool) / sum(land[, "burnable"]) # OK
sum(fire_cool_m %>% as.numeric) / sum(land[, "burnable"])
# OK

# compare
print(time_fool)
print(time_cool)
print(time_cool_m)

# quotient:
(as.numeric(time_cool / 60) / as.numeric(time_fool)) # fool is in min, cool is in sec.
# 0.05663786

# mat over vector
(as.numeric(time_cool_m) / as.numeric(time_cool)) # fool is in min, cool is in sec.



# __________________________________

# with microbenchmark


mbm_real <- microbenchmark(
  fool = simulate_fire_fool_cpp(
    landscape = land[, 1:7],
    burnable = land[, "burnable"],
    ignition_cells = ig_location - 1,
    n_rowcol = c(nrow(land_raster), ncol(land_raster)),
    coef = coefs,
    wind_column = 6 - 1,
    elev_column = 7 - 1,
    distances = distances,
    upper_limit = 1
  ),
  cool_vec = simulate_fire_cpp(
    landscape = land[, 1:7],
    burnable = land[, "burnable"],
    ignition_cells = ig_location - 1,
    n_rowcol = c(nrow(land_raster), ncol(land_raster)),
    coef = coefs,
    wind_column = 6 - 1,
    elev_column = 7 - 1,
    distances = distances,
    upper_limit = 1
  ),
  cool_mat = simulate_fire_mat_cpp(
    landscape = land_arr[, , 1:7],
    burnable = land_arr[, , "burnable"],
    ignition_cells = ig_location_mat - 1,
    coef = coefs,
    wind_layer = 6 - 1,
    elev_layer = 7 - 1,
    distances = distances,
    upper_limit = 1
  ),
  times = 1
)

mbm_all; mbm_real

# Unit: seconds

# Whole landscape is burnable
# función   tiempo
# fool      293.49499  #
# cool_vec  46.10568   #
# cool_mat  11.21112   # 11.21112 / 46.10568 = 0.2431614
#                      # 11.21112 / 293.49499 = 0.03819868

# Real burnable landscape
# función   tiempo
# fool      529.725118
# cool_vec  28.034889  #
# cool_mat  6.996406   # 6.996406 / 28.034889 = 0.2495607
                       # 6.996406 / 529.725118 = 0.01320762


