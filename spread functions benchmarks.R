# benchmarks

library(terra)
library(tidyverse)
library(Rcpp)
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



# Benchmark


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


# compare
print(time_cool_m)


# __________________________________

# with microbenchmark


mbm_all <- microbenchmark(
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
sum(fire_cool_m %>% as.numeric) / sum(land[, "burnable"])
# OK

# compare
print(time_cool_m)

# __________________________________

# with microbenchmark


mbm_real <- microbenchmark(
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


