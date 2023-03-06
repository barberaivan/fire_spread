options(scipen = 999) # turn off scientific notation

library(terra)
library(tidyverse)
library(Rcpp)
library(microbenchmark)
theme_set(theme_bw())

sourceCpp("spread_functions.cpp")

# Data preparation --------------------------------------------------------

## the "fire_spread_data" folder must be located in the same directory as the fire_spread
## repo is.

# get path for data
local_dir <- normalizePath(getwd(), winslash = "\\", mustWork = TRUE)
dir_split <- strsplit(local_dir, .Platform$file.sep)[[1]]
# replace the "fire_spread" directory by "data"
dir_split[length(dir_split)] <- "fire_spread_data"
data_path <- paste(dir_split, collapse = .Platform$file.sep)

land_path <- file.path(data_path, "focal fires data", "data_cholila_landscape.rds")
elev_path <- file.path(data_path, "focal fires data", "data_cholila_elevation.tif")

# import landscape (cholila)
land <- readRDS(land_path)

# cholila raster file

land_raster <- rast(elev_path)
assertthat::assert_that(nrow(land) == ncell(land_raster))


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

# Benchmark cholila saving results ----------------------------------------

# All terrain is burnable
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
"With all terrain burnable"
mbm_all


# Only real burnable terrain is burnable
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

"With only real burnable terrain"
mbm_real

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



# Internal benchmark ------------------------------------------------------

# To evaluate whether further improvements are achievable.

# Create landscape.
coefs <- c(1000, rep(0, 8))
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
size <- 50
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




# sourceCpp("spread_functions.cpp")
# Remove clocks from latest function, rename as _mat2 and compare with _mat,
# which is the one using "continue".

mbm2 <- microbenchmark(
  times = 1000,
  continue = simulate_fire_mat_deterministic_cpp(
    landscape = landscape_arr,
    burnable = matrix(1, nrow(landscape), ncol(landscape)),
    ignition_cells = ig_location,
    coef = coefs,
    wind_layer = wind_column - 1,
    elev_layer = elev_column - 1,
    distances = distances,
    upper_limit = 1.0
  )
)

mbm2