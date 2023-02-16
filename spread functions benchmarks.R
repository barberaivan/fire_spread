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





#### OLD CODE BELOW, NOT TESTED ###########################################
# (BENCHMARKS) ------------------------------------------------------------
# When the data for all fires is downloaded and ready, check the RAM usage of 
# the simulations. If we have RAM unused, we can test the function with 
# precomputed neighbours.


# Compare spread with and without precomputed neighbours ------------------

neighs_mat <- adjacent(r_predictors, 1:ncell(r_predictors),
                       directions = "queen") %>% as.matrix()
typeof(neighs_mat)
# replace NA with cpp NA_INTEGER
na_value <- nrow(neighs_mat) + 2
for(i in 1:ncol(neighs_mat)) {
  # i = 1
  neighs_mat[, i] <- as.integer(neighs_mat[, i])
  neighs_mat[is.na(neighs_mat[, i]), i] <- na_value
}
typeof(neighs_mat)

# object.size(neighs_mat) / 1e6 # 327 Mb for ncol = nrow = 1600

set.seed(seed)
system.time(
  burn_sim_cpp <- simulate_fire_cpp(landscape = values(r_predictors),
                                    burnable = rep(1, ncell(r_predictors)),
                                    ignition_cells = burning_cells - 1,
                                    n_rowcol = c(nrow(r_predictors), ncol(r_predictors)),
                                    coef = coefs,
                                    moves = moves,
                                    wind_column = wind_column - 1,
                                    elev_column = elev_column - 1,
                                    distances = distances,
                                    angles = as.numeric(angles),
                                    upper_limit = 1.0)
)


# sourceCpp("spread_functions.cpp")
set.seed(seed)
system.time(
  burn_sim_cpp_notadj <- simulate_fire_cpp_notadj(
    landscape = values(r_predictors),
    neighbours_matrix = neighs_mat - 1,
    burnable = rep(1, ncell(r_predictors)),
    ignition_cells = burning_cells - 1,
    n_rowcol = c(nrow(r_predictors), ncol(r_predictors)),
    coef = coefs,
    moves = moves,
    wind_column = wind_column - 1,
    elev_column = elev_column - 1,
    distances = distances,
    angles = as.numeric(angles),
    upper_limit = 1.0)
)

# size   //  computing neighs // neighs precomputed   //   quotient
#   600  //     5.475         //       4.549          //  1.203561
#   1200 //     32.013        //       31.906         //  1.003354
#
all.equal(burn_sim_cpp, burn_sim_cpp_notadj)
sum(burn_sim_cpp) / nrow(neighs_mat)

# bug: el notadj no quema nada, por eso es tan fast.




# Serious bechmarking -----------------------------------------------------

set_landscape <- function(size) {
  
  n_rows <- size
  n_cols <- size
  
  r_veg <- rast(ncol = n_cols, nrow = n_rows,
                xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000)
  r_elev <<- r_veg
  r_aspect <<- r_veg
  r_wind <<- r_veg
  r_pp <<- r_veg
  r_temp <<- r_veg
  
  values(r_veg) <- rbinom(ncell(r_veg), prob = 0.5, size = 1)
  # # to test elevation:
  # elevs <- matrix(NA, nrow(r_veg), ncol(r_veg))
  # elevs[1, ] <- 2000
  # for(i in 2:nrow(elevs)) elevs[i, ] <- elevs[i-1, ] - 2000 / nrow(r_veg)
  # as.numeric(t(elevs))#
  values(r_elev) <<- rnorm(ncell(r_veg), 1000, 100) %>% abs
  values(r_aspect) <<- cos(runif(ncell(r_veg), 0, 2 * pi) - 315 * pi / 180)
  values(r_wind) <<- runif(ncell(r_veg), 0, 2 * pi)
  values(r_pp) <<- runif(ncell(r_veg), 0, 1)
  values(r_temp) <<- runif(ncell(r_veg), 0, 1)
  
  r_predictors <<- c(r_veg, r_aspect, r_pp, r_temp, r_elev, r_wind)
  names(r_predictors) <<- names(coefs)[-1]
  
  neighs_mat <<- adjacent(r_predictors, 1:ncell(r_predictors),
                          directions = "queen") %>% as.matrix()
  typeof(neighs_mat)
  # replace NA with cpp NA_INTEGER
  na_value <<- nrow(neighs_mat) + 2
  for(i in 1:ncol(neighs_mat)) {
    # i = 1
    neighs_mat[, i] <<- as.integer(neighs_mat[, i])
    neighs_mat[is.na(neighs_mat[, i]), i] <<- na_value
  }
  
  # create fire raster
  burn <<- r_veg
  values(burn) <<- 0 #rpois(ncell(burn), lambda = 1)
  
  # ignition point:
  ig_location <<- cellFromRowCol(burn, nrow(burn) / 2, ncol(burn) / 2)
  # set as burning
  values(burn)[ig_location:(ig_location+5)] <- 1
  # plot(r_predictors)
  
  # initialize burning_cells cell id
  burning_cells <<- which(values(burn) == 1) # cell number id
}


run_cpp_notadj <- function() {
  
  burn_sim_cpp_notadj <- simulate_fire_cpp_notadj(
    landscape = values(r_predictors),
    neighbours_matrix = neighs_mat - 1,
    burnable = rep(1, ncell(r_predictors)),
    ignition_cells = burning_cells - 1,
    n_rowcol = c(nrow(r_predictors), ncol(r_predictors)),
    coef = coefs,
    wind_column = wind_column - 1,
    elev_column = elev_column - 1,
    distances = distances,
    upper_limit = 1.0)
  
  return(burn_sim_cpp_notadj)
}

run_cpp_adj <- function() {
  burn_sim_cpp <- simulate_fire_cpp(landscape = values(r_predictors),
                                    burnable = rep(1, ncell(r_predictors)),
                                    ignition_cells = burning_cells - 1,
                                    n_rowcol = c(nrow(r_predictors), ncol(r_predictors)),
                                    coef = coefs,
                                    wind_column = wind_column - 1,
                                    elev_column = elev_column - 1,
                                    distances = distances,
                                    upper_limit = 1.0)
  
  return(burn_sim_cpp)
}
# run_cpp_adj
coefs[1] <- -6 # initially set as 0.5, but increased to burn everything
set_landscape(600)
mbm <- microbenchmark(
  not_adj = run_cpp_notadj(),
  adj = run_cpp_adj(),
  times = 500
)
# mbm
autoplot(mbm)

# con size 600, p ~ 1:
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# not_adj 3.455129 3.522688 3.724965 3.714426 3.869174 4.145065    10
# adj     4.300220 4.453728 4.695394 4.575410 4.826309 5.339491    10

# Es más rápido, pero no es una locura la mejora.

# con size 600, p ~ 0:
# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval
# not_adj 29.35553 29.80398 39.09290 31.20053 49.48754 278.7073   500
# adj     18.57764 18.83971 26.06469 19.04619 25.33430 281.7620   500

# calcular los vecinos es faster cuando el fuego no quema casi nada.
