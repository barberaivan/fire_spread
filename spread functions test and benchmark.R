# Tests for the spread functions, written in both cpp and R

library(Rcpp)
library(terra)
library(tidyverse)
library(microbenchmark)

# load cpp and R functions:
sourceCpp("spread_functions.cpp")
source("spread_functions.R")


# Test adjacency related functions ----------------------------------------

# compare the functions written in R by me, the ones from terra and the cpp ones
# written with 0 and 1 indexing.

# data for testing
rowcol <- c(5, 6) # number of rows and columns of the landscape
rtest <- rast(ncol = rowcol[2], nrow = rowcol[1],
              xmin = -1000, xmax = 1000,
              ymin = -1000, ymax = 1000)

# cell to rowcol

# check visually whether results make sense or not, running the following lines 
# a few times

# start
print("--------------------------------------------------------------------")
matrix(1:prod(rowcol), byrow = TRUE, ncol = rowcol[2]) # matrix with cell ids
matrix(1:prod(rowcol), byrow = TRUE, ncol = rowcol[2]) - 1 # with 0-indexing

cell <- sample(1:prod(rowcol), size = 1)
print(paste("cell: ", cell))

cell_to_rowcol(cell, rowcol)                 # R
cell_to_rowcol_cpp(cell, rowcol)             # cpp
terra::rowColFromCell(rtest, cell) %>% t     # terra
cell_to_rowcol_cpp0(cell - 1, rowcol)        # cpp 0-indexing
# end


# rowcol to cell

print("--------------------------------------------------------------------")
matrix(1:prod(rowcol), byrow = TRUE, ncol = rowcol[2])
matrix(1:prod(rowcol), byrow = TRUE, ncol = rowcol[2]) - 1

row_col_id <- matrix(
  c(sample(1:rowcol[1], size = 1),
    sample(1:rowcol[2], size = 1)),
  ncol = 1
)
print(paste("row:", row_col_id[1,], ", col:", row_col_id[2,]))

rowcol_to_cell(row_col_id, rowcol)
rowcol_to_cell_cpp(row_col_id, rowcol)
terra::cellFromRowCol(rtest, row_col_id[1], row_col_id[2]) 
rowcol_to_cell_cpp0(row_col_id - 1, rowcol)



# adjacent

print("--------------------------------------------------------------------")
cell <- sample(1:prod(rowcol), size = 1)
print(paste("cell: ", cell))

matrix(1:prod(rowcol), byrow = TRUE, ncol = rowcol[2])
matrix(1:prod(rowcol), byrow = TRUE, ncol = rowcol[2]) - 1

adjacent_r(cell, rowcol)
adjacent_cpp(cell, rowcol)
terra::adjacent(rtest, cell, directions = "queen")
adjacent_cpp0(cell - 1, rowcol)



# spread_around test ------------------------------------------------------


# create data for testing
# model coefficients (includes intercept)
coefs <- c(0.5, 
           -1.5, -1.3, -0.5, 
           1, 1, 1, 1, 1)
names(coefs) <- c("intercept", 
                  "subalpine", "wet", "dry",
                  "fwi", 
                  "aspect", 
                  "wind", 
                  "elev",
                  "slope")
### IMPORTANT  ---> wind and elevation parameters must be the last ones.

elev_column <- which(names(coefs) == "elev") - 1
wind_column <- which(names(coefs) == "wind") - 1 # -1 because the design matrix will have no intercept

# landscape raster
size <- 30
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

# vector of distances between a cell and its neighbours
# (used to compute slope effect)
distances <- rep(res(landscape)[1], 8) # sides
distances[c(1, 3, 6, 8)] <- res(landscape)[1] * sqrt(2)

# turn landscape into matrix
predmat <- values(landscape)


# TEST probability output

# start
burning_cell <- sample(1:ncell(landscape), size = 1)
neighs_raw <- adjacent(landscape, 1, "queen") %>% as.numeric
not_na <- !is.na(neighs_raw)
pos <- (1:8)[not_na]
neighs_id <- neighs_raw[not_na]

s <- round(runif(1, 0, 10000)) # get seed

set.seed(s) # set seed to compare the random simulation for burning
spread_result_r <- spread_around_r(
  data_burning = predmat[burning_cell, , drop = F], # USE DROP!!
  data_neighbours = predmat[neighs_id, , drop = F],
  coef = coefs, 
  positions = pos,
  distances = distances,
  upper_limit = 1
)

set.seed(s)
spread_result_cpp_burn <- spread_around_cpp(
  data_burning = predmat[burning_cell, , drop = F], # USE DROP!!
  data_neighbours = predmat[neighs_id, , drop = F],
  coef = coefs, 
  positions = pos - 1, # -1 for 0-indexing
  distances = distances,
  upper_limit = 1,
  wind_column = which(names(landscape) == "wind") - 1,
  elev_column = which(names(landscape) == "elev") - 1  # -1 for 0-indexing
)

spread_result_cpp_prob <- spread_around_prob_cpp(
  data_burning = predmat[burning_cell, , drop = F], # USE DROP!!
  data_neighbours = predmat[neighs_id, , drop = F],
  coef = coefs, 
  positions = pos - 1, # -1 for 0-indexing
  distances = distances,
  upper_limit = 1,
  wind_column = which(names(landscape) == "wind") - 1,
  elev_column = which(names(landscape) == "elev") - 1  # -1 for 0-indexing
)

# compare
all.equal(spread_result_r[, "probs"], spread_result_cpp_prob)
all.equal(spread_result_r[, "burn"], spread_result_cpp_burn)
# end


# simulate_fire tests -----------------------------------------------------


# create data for testing
# model coefficients (includes intercept)
coefs <- c(0.5, 
           -1.5, -1.3, -0.5, 
           1, 1, 1, 1, 1)
names(coefs) <- c("intercept", 
                  "subalpine", "wet", "dry",
                  "fwi", 
                  "aspect", 
                  "wind", 
                  "elev",
                  "slope")
### IMPORTANT  ---> wind and elevation parameters must be the last ones.

elev_column <- which(names(coefs) == "elev") - 1
wind_column <- which(names(coefs) == "wind") - 1 # -1 because the design matrix will have no intercept

# landscape raster
size <- 30
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

# vector of distances between a cell and its neighbours
# (used to compute slope effect)
distances <- rep(res(landscape)[1], 8) # sides
distances[c(1, 3, 6, 8)] <- res(landscape)[1] * sqrt(2)

# make ignition point(s)
ig_location <- cellFromRowCol(landscape, nrow(landscape) / 2, ncol(landscape) / 2)

s <- round(runif(1, 1, 20000))

set.seed(s)
burn_result_r <- simulate_fire_r(
  landscape = landscape,
  burnable = rep(1, ncell(landscape)),
  ignition_cells = ig_location,
  n_rowcol = c(nrow(landscape), ncol(landscape)),
  coef = coefs,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
)

set.seed(s)
# This function plots the spread, only useful for testing purposes
burn_result_r_plot <- simulate_fire_plot(
  landscape = landscape,
  burnable = rep(1, ncell(landscape)),
  ignition_cells = ig_location,
  n_rowcol = c(nrow(landscape), ncol(landscape)),
  coef = coefs,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
)

set.seed(s)
burn_result_cpp <- simulate_fire_cpp(
  landscape = values(landscape),      # to cpp we pass the values, not the raster
  burnable = rep(1, ncell(landscape)),
  ignition_cells = ig_location - 1,
  n_rowcol = c(nrow(landscape), ncol(landscape)),
  coef = unname(coefs),
  wind_column = wind_column - 1,
  elev_column = elev_column - 1,
  distances = distances,
  upper_limit = 1.0
)

all.equal(burn_result_r, burn_result_r_plot)
all.equal(burn_result_r, burn_result_cpp)


# Testing for effects in simulate_fire -------------------------------------

# simulate landscapes where one layer varies and check that results are OK.

# vegetation tests # 

landscape_base <- landscape
landscape_base$subalpine <- 0
landscape_base$wet <- 0
landscape_base$dry <- 0 # it's all shrubland
landscape_base$fwi <- 0
landscape_base$aspect <- 0
landscape_base$wind <- 0 # north wind
landscape_base$elev <- elevation_mean

# subalpine
lands_sub <- landscape_base
lands_sub$subalpine[1:(ncell(landscape)/2)] <- 1 # the northern half is subalpine
c_sub <- coefs
c_sub["subalpine"] <- -30 # subalpine is not flammable
c_sub["intercept"] <- 10 # shrubland is very flammable
c_sub["wind"] <- 0 # remove wind effect

# Upper half is not flammable, lower is highly flammable
simulate_fire_plot(
  landscape = lands_sub,
  burnable = rep(1, ncell(lands_sub)),
  ignition_cells = ig_location,
  n_rowcol = c(nrow(lands_sub), ncol(lands_sub)),
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) # good

# wet
lands_sub <- landscape_base
lands_sub$wet[1:(ncell(landscape)/2)] <- 1 # the northern half is subalpine
c_sub <- coefs
c_sub["wet"] <- -12 # subalpine is not flammable
c_sub["intercept"] <- 10 # shrubland is very flammable
c_sub["wind"] <- 0 # remove wind effect

# Upper half is not flammable, lower is highly flammable
simulate_fire_plot(
  landscape = lands_sub,
  burnable = rep(1, ncell(lands_sub)),
  ignition_cells = ig_location,
  n_rowcol = c(nrow(lands_sub), ncol(lands_sub)),
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) # good


# wind test
wind_dir <- 270               # direction, in angles, from which the wind comes
lands_sub <- landscape_base
lands_sub$wind <- rep(wind_dir * pi / 180, ncell(lands_sub))
c_sub <- rep(0, length(coefs))
names(c_sub) <- names(coefs)
c_sub["wind"] <- 3 # increase wind effect
c_sub["intercept"] <- 0 

simulate_fire_plot(
  landscape = lands_sub,
  burnable = rep(1, ncell(lands_sub)),
  ignition_cells = ig_location,
  n_rowcol = c(nrow(lands_sub), ncol(lands_sub)),
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) 


# slope test
lands_sub <- landscape_base
lands_sub$elev <- rep(seq(2000, 1000, length.out = nrow(lands_sub)),
                      each = ncol(lands_sub)) # fill values by row
c_sub <- rep(0, length(coefs))
names(c_sub) <- names(coefs)
c_sub["slope"] <- 3 # increase slope effect
c_sub["intercept"] <- -1 

# elevation increases from below, so the fire should spread upwards because of 
# slope effect (elevation effect was removed.)
simulate_fire_plot(
  landscape = lands_sub,
  burnable = rep(1, ncell(lands_sub)),
  ignition_cells = ig_location,
  n_rowcol = c(nrow(lands_sub), ncol(lands_sub)),
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) 


# elevation test
lands_sub <- landscape_base
lands_sub$elev <- rep(seq(2000, 1000, length.out = nrow(lands_sub)),
                      each = ncol(lands_sub)) # fill values by row
c_sub <- rep(0, length(coefs))
names(c_sub) <- names(coefs)
c_sub["elev"] <- -3 # high elevation effect (higher is less flammable)
c_sub["intercept"] <- 1

# elevation increases from below, so the fire should spread downwards because of 
# the elevation effect, which is negative. Here the slope effect was removed.
simulate_fire_plot(
  landscape = lands_sub,
  burnable = rep(1, ncell(lands_sub)),
  ignition_cells = ig_location,
  n_rowcol = c(nrow(lands_sub), ncol(lands_sub)),
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) 

# aspect test
lands_sub <- landscape_base
lands_sub$aspect <- rep(seq(1, -1, length.out = nrow(lands_sub)),
                        each = ncol(lands_sub)) # fill values by row
c_sub <- rep(0, length(coefs))
names(c_sub) <- names(coefs)
c_sub["aspect"] <- 2 # increase aspect effect
c_sub["intercept"] <- -1 

# Fire should spread upwards
simulate_fire_plot(
  landscape = lands_sub,
  burnable = rep(1, ncell(lands_sub)),
  ignition_cells = ig_location,
  n_rowcol = c(nrow(lands_sub), ncol(lands_sub)),
  coef = c_sub,
  wind_column = wind_column,
  elev_column = elev_column,
  distances = distances,
  upper_limit = 1.0
) 


# Tests are OK. 





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