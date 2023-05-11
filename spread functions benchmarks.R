options(scipen = 999) # turn off scientific notation

library(terra)
library(tidyverse)
library(Rcpp)
library(microbenchmark)
theme_set(theme_bw())

sourceCpp("spread_functions.cpp")
source("spread_functions.R") # for land_cube()

# Data preparation --------------------------------------------------------

land <- readRDS(file.path("data", "focal fires data", 
                          "landscapes_ig-known_non-steppe", "2015_50.rds"))
land_raster <- rast(file.path("data", "focal fires data", "raw data from GEE",
                              "fire_data_raw_2015_50.tif"))

# cholila raster file
assertthat::assert_that(
  (ncol(land$vegetation) * nrow(land$vegetation)) == ncell(land_raster)
)

n_veg_types <- 6
n_terrain <- 4

# coefficients
coefs <- c(
  rep(1000, n_veg_types),
  rep(0, n_terrain)
)

ig_location <- cellFromRowCol(land_raster,
                              nrow(land_raster) / 2, ncol(land_raster) / 2)

ig_location_mat <- rowColFromCell(land_raster, ig_location) %>% t


# Benchmark cholila saving results ----------------------------------------

# make vegetation layer with all the terrain burnable.
vegetation_burn_all <- land$vegetation
vegetation_burn_all[vegetation_burn_all == 99] <- 0
hist(as.vector(vegetation_burn_all))
hist(as.vector(land$vegetation))
# ok

# All terrain is burnable
mbm_all <- microbenchmark(
  burns_all = {
    fire_all <- simulate_fire_cpp(
      terrain = land$terrain, 
      vegetation = vegetation_burn_all,
      ignition_cells = ig_location_mat - 1,
      coef = coefs,
      n_veg_types = n_veg_types,
      upper_limit = 1
    )
  },
  times = 5
)
"With all terrain burnable"
mbm_all

# Only real burnable terrain is burnable
mbm_real <- microbenchmark(
  burns_plants = {
    fire_real <- simulate_fire_cpp(
      terrain = land$terrain, 
      vegetation = land$vegetation,
      ignition_cells = ig_location_mat - 1,
      coef = coefs,
      n_veg_types = n_veg_types,
      upper_limit = 1
    )
  },
  times = 5
)

"With only real burnable terrain"
mbm_real

# Unit: seconds

# "With all terrain burnable"
# > mbm_all
# expr          min       lq     mean   median       uq      max    neval
# burns_all 1.806337 1.809189 1.818991 1.811668 1.819008 1.848751     5

# "With only real burnable terrain"
# > mbm_real
# expr            min       lq     mean   median       uq      max     neval
# burns_plants 1.162128 1.169062 1.168091 1.169391 1.169741 1.170131     5

# check:
sum(as.vector(fire_all)) / ncell(land_raster)
sum(as.vector(fire_real)) / ncell(land_raster)

# Comments: before using floats and representing vegetation in a single layer,
# it took 11 s to burn the whole landscape, and ~7 s to burn all vegetated 
# pixels. The times were reduced by a factor of
1.811668 / 11 # 0.1646971
1.169391 / 7  # 0.1670559

# For the worst-timing particles in this landscape, simulating 10 fires in 1000 
# particles would take approximately
1.169391 * 10 * 1000 / 3600 # 3.248308 h

# If parallelized, assuming a reduction factor of 0.2, it would take
1.169391 * 10 * 1000 * 0.2 / 60 # 38.9797 min