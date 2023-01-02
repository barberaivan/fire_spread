# This script prepares the landscapes where fires are gonna be simulated. 
# First, I use only the landscape for Cholila, which is the largest. I will test
# whether its size is sufficient to detect large discrepancies.

library(terra)
data_dir <- "/home/ivan/Insync/Fire spread modelling/data/focal fires data/"

# Cholila fire ------------------------------------------------------------

land0 <- rast(paste(data_dir, "data_cholila.tif", sep = ""))
elev <- rast(paste(data_dir, "data_cholila_elevation.tif", sep = ""))
angles_90m <- rast(paste(data_dir, "data_cholila_elevation_284_10_90m_ang.asc", sep = ""))

# interpolate windir to 30 m
windir <- project(angles_90m, elev, method = "cubicspline")

# make vegetation data 
veg_codes <- values(land0$veg)
unique(veg_codes) 
# replace nan with 1 (non burnable)
veg_codes[is.nan(veg_codes)] <- 1

