# benchmarks

library(terra)
library(tidyverse)
library(Rcpp)
library(viridis)
theme_set(theme_bw())

sourceCpp("spread_functions.cpp")

# Data preparation --------------------------------------------------------

# import landscape (cholila)
data_dir <- "/home/ivan/Insync/Fire spread modelling/data/focal fires data/"
land <- readRDS(paste(data_dir, "data_cholila_landscape.R", sep = ""))

# cholila raster file
land_raster <- rast(paste(data_dir, "data_cholila_elevation.tif", sep = ""))
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


# Benchmark fool vs cool versions, whole landscape burnable ----------------

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


# cool version
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

# check both fires were large
sum(fire_cool) / length(fire_cool)
sum(fire_fool) / length(fire_fool)

# compare
print(time_fool)
print(time_cool)

# quotient:
(as.numeric(time_cool / 60) / as.numeric(time_fool)) # fool is in min, cool is in sec.
# 0.164394



# Benchmark fool vs cool versions, whole landscape burnable ----------------

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

# check both fires were large
sum(fire_cool) / sum(land[, "burnable"])
sum(fire_fool) / sum(land[, "burnable"]) # OK

# compare
print(time_fool)
print(time_cool)

# quotient:
(as.numeric(time_cool / 60) / as.numeric(time_fool)) # fool is in min, cool is in sec.
# 0.05663786

# La mejora es bestial, y más aún en paisajes reales, con obstáculos

