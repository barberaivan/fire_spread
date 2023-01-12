# This script prepares the landscapes where fires are gonna be simulated. 
# First, I use only the landscape for Cholila, which is the largest. I will test
# whether its size is sufficient to detect large discrepancies.


#´ The layers in data_ are:
#´   subalpine forest, {0, 1} (shrubland goes in the intercept)
#´   wet forest,       {0, 1}
#´   dry forest,       {0, 1}
#´   fwi,              (-Inf, +Inf) (standardized anomaly at pixel level)
#´   aspect            [-1, 1] # Northwestyness
#´   wind direction    [-1, 1] (windward vs leeward)
#´   elevation,        (standardized in the code)


library(terra)
library(tidyverse)
data_dir <- "/home/ivan/Insync/Fire spread modelling/data/focal fires data/"

# Cholila fire ------------------------------------------------------------

land0 <- rast(paste(data_dir, "data_cholila.tif", sep = ""))
elev <- rast(paste(data_dir, "data_cholila_elevation.tif", sep = ""))
angles_90m <- rast(paste(data_dir, "data_cholila_elevation_284_10_90m_ang.asc", sep = ""))

# interpolate windir to 30 m
windir <- project(angles_90m, elev, method = "cubicspline")

# make vegetation data 

#                                 class code
# 1                        Non burnable    1
# 2                    Subalpine forest    2
# 3                          Wet forest    3
# 4                          Dry forest    4
# 5                           Shrubland    5
# 6                           Grassland    6 (turn into shrubland)
# 7 Anthropogenic prairie and shrubland    7 (turned into shrubland)
# 8                          Plantation    8 (turned into dry forest)

veg_codes <- values(land0$veg)
unique(veg_codes) 
# replace nan with 1 (non burnable)
veg_codes[is.nan(veg_codes)] <- 1

veg_mat <- matrix(0, length(veg_codes), 3) # 3 cols for wet, subalpine, and dry forests
                                            # shrubland is the reference
colnames(veg_mat) <- c("subalpine", "wet", "dry")

veg_mat[veg_codes == 3, 1] <- 1 # subalpine 
veg_mat[veg_codes == 3, 2] <- 1 # wet
veg_mat[veg_codes == 4, 3] <- 1 # dry

burnable <- as.numeric(veg_codes != 1)

# make landscape
land0_mat <- values(land0)
colnames(land0_mat)

# IMPORTANT: ADD BURNABLE AND BURNED LAYERS IN THE END

landscape <- cbind(
  veg_mat, 
  fwi = rep(0, nrow(land0_mat)),
  aspect = cos(land0_mat[, "aspect"] * pi / 180), # northing
  windir = values(windir)[, 1],
  elev = land0_mat[, "elev"],
  burnable = burnable,
  burned = land0_mat[, "burned"]
)

head(landscape)
dim(landscape)
object.size(landscape) / 1e6

sum(landscape[, "burnable"]) / nrow(landscape) # 68 % quemable

# save cholila landscape

saveRDS(landscape, paste(data_dir, "data_cholila_landscape.R", sep = ""))





# 2021-865 landscape ------------------------------------------------------

data_dir <- "/home/ivan/Insync/Fire spread modelling/data/focal fires data/"

land0 <- rast(paste(data_dir, "data_2021_865.tif", sep = ""))
# elev <- rast(paste(data_dir, "data_2021_865_elevation.tif", sep = ""))
# angles_90m <- rast(paste(data_dir, "data_2021_865_elevation_284_10_90m_ang.asc", sep = ""))

# interpolate windir to 30 m
# windir <- project(angles_90m, elev, method = "cubicspline")

# make vegetation data 

#                                 class code
# 1                        Non burnable    1
# 2                    Subalpine forest    2
# 3                          Wet forest    3
# 4                          Dry forest    4
# 5                           Shrubland    5
# 6                           Grassland    6 (turn into shrubland)
# 7 Anthropogenic prairie and shrubland    7 (turned into shrubland)
# 8                          Plantation    8 (turned into dry forest)

veg_codes <- values(land0$veg)
unique(veg_codes) 
# replace nan with 1 (non burnable)
veg_codes[is.nan(veg_codes)] <- 1

veg_mat <- matrix(0, length(veg_codes), 3) # 3 cols for wet, subalpine, and dry forests
# shrubland is the reference
colnames(veg_mat) <- c("subalpine", "wet", "dry")

veg_mat[veg_codes == 3, 1] <- 1 # subalpine 
veg_mat[veg_codes == 3, 2] <- 1 # wet
veg_mat[veg_codes == 4, 3] <- 1 # dry

burnable <- as.numeric(veg_codes != 1)

# make landscape
land0_mat <- values(land0)
colnames(land0_mat)

# IMPORTANT: ADD BURNABLE AND BURNED LAYERS IN THE END

landscape <- cbind(
  veg_mat, 
  fwi = rep(0, nrow(land0_mat)),
  aspect = rep(0, nrow(land0_mat)), #cos(land0_mat[, "aspect"] * pi / 180), # northing
  windir = rep(0, nrow(land0_mat)), #values(windir)[, 1],
  elev = land0_mat[, "elev"],
  burnable = burnable,
  burned = land0_mat[, "burned"]
)

head(landscape)
dim(landscape)
object.size(landscape) / 1e6

sum(landscape[, "burnable"]) / nrow(landscape) # 84 % quemable

# save 2021_865 landscape

saveRDS(landscape, paste(data_dir, "data_2021_865_landscape.R", sep = ""))


