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
  burns_all_raw = {
    fire_all <- simulate_fire_cpp(
      terrain = land$terrain, 
      vegetation = vegetation_burn_all,
      ignition_cells = ig_location_mat - 1,
      coef = coefs,
      n_veg_types = n_veg_types,
      upper_limit = 1
    )
  },
  
  burns_all_pak = {
    fire_all <- FireSpread::simulate_fire(
      terrain = land$terrain, 
      vegetation = vegetation_burn_all,
      ignition_cells = ig_location_mat - 1,
      coef = coefs,
      n_veg_types = n_veg_types,
      upper_limit = 1
    )
  },
  
  times = 10
)
mbm_all
# With all terrain burnable
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval cld
# burns_all_raw 1.724945 1.731835 1.757202 1.746229 1.759198 1.848571    10   b
# burns_all_pak 1.488216 1.494799 1.502134 1.499918 1.509534 1.525945    10  a 
1.494799 / 1.731835 # 0.86 


# Only real burnable terrain is burnable
mbm_real <- microbenchmark(
  burns_plants_raw = {
    fire_real <- simulate_fire_cpp(
      terrain = land$terrain, 
      vegetation = land$vegetation,
      ignition_cells = ig_location_mat - 1,
      coef = coefs,
      n_veg_types = n_veg_types,
      upper_limit = 1
    )
  },
  
  burns_plants_pak = {
    fire_real <- FireSpread::simulate_fire(
      terrain = land$terrain, 
      vegetation = land$vegetation,
      ignition_cells = ig_location_mat - 1,
      coef = coefs,
      n_veg_types = n_veg_types,
      upper_limit = 1
    )
  },
  
  times = 50
)

mbm_real
# With only real burnable terrain
# Unit: milliseconds
# expr       min        lq     mean    median       uq      max neval cld
# burns_plants_raw 1101.3957 1115.6209 1128.355 1123.6498 1136.985 1206.459    50   b
# burns_plants_pak  972.1404  982.6578 1005.346  996.8097 1009.942 1198.779    50  a 
982.6578 / 1101.3957

# las optimizaciones con floats bajaron el tiempo en un factor de ~ 0.87

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
0.9826578 * 10 * 1000 / 3600 # 2.72 h

# If parallelized, assuming a reduction factor of 0.2, it would take
0.9826578 * 10 * 1000 * 0.2 / 60 # 32.76 min
