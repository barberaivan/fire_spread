# Using cholila landscape I'll evaluate two aspects:
# 1) computation time to burn all the landscape as a function of landscape size,
# 2) overlap of the largest fire with cholila, as a function of landscape size.

# Landscapes of varying size are created the following way:
# 1) compute the fire's bounding box,
# 2) make a buffer = sqrt(area(bounding_box(fire))) * p,
# 3) vary landscape size with p = seq(1, 0.2, by = 0.1).

# That was made on GEE, and the polygons of varying size are in the file
# landscapes_cholila.shp.
# They will be used to crop the full landscape (p = 1), which data is in the
# file data_cholila_landscape.R. For its raster we'll use the 
# data_cholila_elevation.tif, so we can make spatial masking with the polygons.

library(terra)
library(tidyverse)
library(Rcpp)
library(viridis)
theme_set(theme_bw())

sourceCpp("spread_functions.cpp")
sourceCpp("discrepancy_functions.cpp")


# Data preparation --------------------------------------------------------


# import landscape
data_dir <- "/home/ivan/Insync/Fire spread modelling/data/focal fires data/"
land <- readRDS(paste(data_dir, "data_cholila_landscape.R", sep = ""))

# vectors to crop landscapes of different sizes
lands_vec0 <- vect(paste(data_dir, "landscapes_cholila.shp", sep = ""))
# plot(lands_vec)

# reproject to posgar 2007 faja1 (flat)
lands_vec <- project(lands_vec0, "EPSG:5343")

# get raster with the largest landscape (prop = 1) 
land_raster <- rast(paste(data_dir, "data_cholila_elevation.tif", sep = ""))
# plot(land_raster)

# Create burnable layers for subsets of the landscape
np <- length(lands_vec)
burnable_mat <- matrix(land[, "burnable"], nrow(land), np)
colnames(burnable_mat) <- lands_vec$proportion

# define mask value to set those pixels as non-burnable
out_value <- -999999

# create burnable layers with variying-size landscapes
nrow(land) == ncell(land_raster)

for(p in 1:np) {
  print(p)
  land_raster_masked <- mask(
    land_raster, 
    lands_vec[p, ], 
    updatevalue = out_value
  )
  
  vals <- values(land_raster_masked)
  masked <- which(vals == out_value)
  burnable_mat[masked, p] <- 0
}

# # check:
# land_raster_1 <- land_raster
# land_raster_2 <- land_raster
# values(land_raster_1) <- burnable_mat[, 1]
# values(land_raster_2) <- burnable_mat[, np]
# plot(land_raster_1)
# plot(land_raster_2)
# # ok

# Data to simulate fires

# distances for 30 m resolution
distances <- rep(res(land_raster)[1], 8) # sides
distances[c(1, 3, 6, 8)] <- res(land_raster)[1] * sqrt(2)

# coefficients 
coefs <- c(1000, 
           0, 0, 0, 
           0, 0, 0, 0, 0)

ig_location <- cellFromRowCol(land_raster, 
                              nrow(land_raster) / 2, ncol(land_raster) / 2)

# see fire and burnable (it's important to use the burnable layer to 
# avoid huge fires). Note that cholila was mainly contained by non-burnable
# landscape.
# values(land_raster) <- land[, "burnable"]
# plot(land_raster, col = c("black", "green"))
# values(land_raster)[land[, "burned"] == 1] <- 2
# plot(land_raster, col = c("black", "green", "red")) # burnable

# check ig_location falls in burnable area
land[ig_location, "burnable"] == 1

# cholila size
# (land[, "burned"] %>% sum()) / nrow(land) * 100 
# 3.047275 % del paisaje



# Simulation time as a function of landscape size -------------------------

times <- numeric(np)

# matrix with largest fires
burned_mat <- matrix(NA, nrow(land), np)

for(p in 2:np) {
  
  print(lands_vec$proportion[p])
  start_time <- Sys.time()
  
  set.seed(124)
  burned_mat[, p] <- simulate_fire_cpp(
      landscape = land[, 1:7],      # to cpp we pass the values, not the raster
                                    # and remove burnable and burned layers
      burnable = burnable_mat[, p],
      ignition_cells = ig_location - 1,
      n_rowcol = c(nrow(land_raster), ncol(land_raster)),
      coef = coefs,
      wind_column = 6 - 1,
      elev_column = 7 - 1,
      distances = distances,
      upper_limit = 1
    )
  
  end_time <- Sys.time()

  times[p] <- as.difftime(end_time - start_time, format = "%X", units = "sec")
  
}

# saveRDS(burned_mat, "/home/ivan/Insync/Fire spread modelling/data/simulations/burned_mat.R")
# saveRDS(times, "/home/ivan/Insync/Fire spread modelling/data/simulations/times.R")

colnames(burned_mat) <- colnames(burnable_mat)
timesize <- data.frame(proportion = lands_vec$proportion,
                       time = times,
                       overlap = NA,
                       burned_size = colSums(burned_mat))

# Compute overlaps
for(p in 1:np) {
  timesize$overlap[p] <- overlap_sp(land[, "burned"], burned_mat[, p])
}

# compute time relative to max
timesize$time_q <- timesize$time / timesize$time[1]

par(mfrow = c(2, 2))
plot(overlap ~ proportion, data = timesize,
     ylim = c(0, max(timesize$overlap)))
plot(time ~ proportion, data = timesize, ylab = "time (min)",
     ylim = c(0, max(timesize$time)))
plot(time_q ~ proportion, data = timesize,
     ylim = c(0, 1))
plot(time_q ~ burned_size, data = timesize,
     ylim = c(0, 1))
par(mfrow = c(1, 1))

# 0.2 is a very high overlap to note differences. perhaps p = 0.6 is good, but
# the computation time is too similar.

# Plot contrasting landscapes
land_raster_1 <- land_raster
land_raster_06 <- land_raster
land_raster_05 <- land_raster
values(land_raster_1) <- burnable_mat[, "1"]
values(land_raster_06) <- burnable_mat[, "0.6"]
values(land_raster_05) <- burnable_mat[, "0.5"]

# # paint cholila
# values(land_raster_1)[land[, "burned"] == 1] <- 3
# values(land_raster_06)[land[, "burned"] == 1] <- 3
# values(land_raster_05)[land[, "burned"] == 1] <- 3

# paint sims
values(land_raster_1)[burned_mat[, "1"] == 1] <- 2
values(land_raster_06)[burned_mat[, "0.6"] == 1] <- 2
values(land_raster_05)[burned_mat[, "0.5"] == 1] <- 2


par(mfrow = c(1, 3))
plot(land_raster_1, main = "1")
plot(land_raster_06, main = "0.6")
plot(land_raster_05, main = "0.5")
par(mfrow = c(1, 1))


# with size = 0.5 it gets an overlap almost = 0.1

# The jump in computation time as a function of proportion is not explained by
# burned_size. But when proportion <= 0.5, a large area of the landscape that is 
# burnable cannot be reached by fires starting in the middle. 
# I will measure the computation time in landscapes where all the landscape is 
# burnable, to test whether the burnable structure affects the time-burned
# relationship.


# Simulation time as a function of landscape size (all burnable) ---------

burnable_mat_all <- matrix(1, nrow(land), np)
colnames(burnable_mat_all) <- lands_vec$proportion

# create burnable layers with variying-size landscapes
for(p in 1:np) {
  print(p)
  land_raster_masked <- mask(
    land_raster, 
    lands_vec[p, ], 
    updatevalue = out_value
  )
  
  vals <- values(land_raster_masked)
  masked <- which(vals == out_value)
  burnable_mat_all[masked, p] <- 0
}

colSums(burnable_mat_all)
colSums(burnable_mat)
ncell(land_raster)
# matrix with largest fires
burned_mat_all <- matrix(NA, nrow(land), np)

times_all <- numeric(np)

for(p in 1:np) {
  
  print(lands_vec$proportion[p])
  start_time <- Sys.time()
  
  set.seed(124)
  burned_mat_all[, p] <- simulate_fire_cpp(
    landscape = land[, 1:7],      # to cpp we pass the values, not the raster
    # and remove burnable and burned layers
    burnable = burnable_mat_all[, p],
    ignition_cells = ig_location - 1,
    n_rowcol = c(nrow(land_raster), ncol(land_raster)),
    coef = coefs,
    wind_column = 6 - 1,
    elev_column = 7 - 1,
    distances = distances,
    upper_limit = 1
  )
  
  end_time <- Sys.time()
  
  print(end_time - start_time) # 3.26 min
  
  
  times_all[p] <- as.difftime(end_time - start_time, format = "%X", units = "sec")
  
}

saveRDS(burned_mat_all, "/home/ivan/Insync/Fire spread modelling/data/simulations/burned_mat_all.R")
saveRDS(times_all, "/home/ivan/Insync/Fire spread modelling/data/simulations/times_all.R")

colnames(burned_mat_all) <- colnames(burnable_mat_all)
timesize_all <- data.frame(proportion = lands_vec$proportion,
                       time = times_all,
                       overlap = NA,
                       burned_size = colSums(burned_mat_all))

# Compute overlaps
for(p in 1:np) {
  timesize_all$overlap[p] <- overlap_sp(land[, "burned"], burned_mat_all[, p])
}

# compute time relative to max
timesize_all$time_q <- timesize_all$time / timesize_all$time[1]

par(mfrow = c(2, 2))
plot(overlap ~ proportion, data = timesize_all,
     ylim = c(0, max(timesize_all$overlap)))
plot(time ~ proportion, data = timesize_all, ylab = "time (min)",
     ylim = c(0, max(timesize_all$time)))
plot(time_q ~ proportion, data = timesize_all,
     ylim = c(0, 1))
plot(time_q ~ burned_size, data = timesize_all,
     ylim = c(0, 1))
par(mfrow = c(1, 1))

# computation time is linear as a function of burned size, and it's much faster 
# (3 min vs 9 min) when the whole landscape is burnable.
# An explanation might be that more burn cycles (while iterations) are needed
# when the terrain is more complex. So, the largest overhead in the function
# is where the while does things independent of burning_size.
# It is probably the subsetting:

# burning[burning == 2] = 0; // turn off previously burned cells
# burned[burning == 1] = 1;  // set as burned the currently burning.

# here a loop runs over the whole burning and burned vectors, which is probably
# unnecessary.




# Time to burn all the 2021_865 fire --------------------------------------

land2 <- readRDS(paste(data_dir, "data_2021_865_landscape.R", sep = ""))
land_raster2 <- rast(paste(data_dir, "data_2021_865.tif", sep = ""))
ig_location2 <- cellFromRowCol(land_raster2, 
                               nrow(land_raster2) / 2, ncol(land_raster2) / 2)
land2[ig_location2, "burnable"] # OK

start_time <- Sys.time()
set.seed(124)
fire2 <- simulate_fire_cpp(
  landscape = land2[, 1:7],      # to cpp we pass the values, not the raster
  # and remove burnable and burned layers
  burnable = land2[, "burnable"],
  ignition_cells = ig_location2 - 1,
  n_rowcol = c(nrow(land_raster2), ncol(land_raster2)),
  coef = coefs,
  wind_column = 6 - 1,
  elev_column = 7 - 1,
  distances = distances,
  upper_limit = 1
)
end_time <- Sys.time()

t2 <- end_time - start_time
# 1.27 min

# size quotient between 2021_865 and cholila:
sum(land2[, "burnable"]) / colSums(burnable_mat)

# 2021_865 overlap with full landscape
(ov2 <- overlap_sp(land2[, "burned"], fire2)) # 0.02650215 # quite low :)

# proportion burned is 1?
sum(fire2) / sum(land2[, "burnable"]) # yes

# plot
lr2 <- land_raster2["burned"]
values(lr2) <- 0
values(lr2)[land2[, "burnable"] == 1] <- 1
values(lr2)[land2[, "burned"] == 1] <- 2
plot(lr2)

