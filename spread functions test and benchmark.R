# Tests for the spread functions written in cpp (through Rcpp), in
# spread_functions.cpp

library(Rcpp)
library(terra)
library(tidyverse)
library(microbenchmark)

sourceCpp("spread_functions.cpp")

# Create data -------------------------------------------------------------

# model coefficients (includes intercept)
coefs <- c(0.5, -1.5, 0.5, 0.5, 0.5, 0.5, 0.5)
names(coefs) <- c("Intercept", "veg", "aspect", "pp", "temp", "elev", "wind")
### IMPORTANT  ---> wind and elevation columns must be the last.

elev_column <- which(names(coefs) == "elev") - 1
wind_column <- which(names(coefs) == "wind") - 1 # -1 because the design matrix will have no intercept

# test wind effect
# coefs <- c(0.5, 0, 0, 0, 3, 0, 0)
# test slope effect
# coefs <- c(0.5, 0, 10, 0, 0, 0, 0)

# predictors raster

size <- 1200

n_rows <- size
n_cols <- size

r_veg <- rast(ncol = n_cols, nrow = n_rows,
              xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000)
r_elev <- r_veg
r_aspect <- r_veg
r_wind <- r_veg
r_pp <- r_veg
r_temp <- r_veg

values(r_veg) <- rbinom(ncell(r_veg), prob = 0.5, size = 1)
# # to test elevation:
# elevs <- matrix(NA, nrow(r_veg), ncol(r_veg))
# elevs[1, ] <- 2000
# for(i in 2:nrow(elevs)) elevs[i, ] <- elevs[i-1, ] - 2000 / nrow(r_veg)
# as.numeric(t(elevs))#
values(r_elev) <- rnorm(ncell(r_veg), 1000, 100) %>% abs
values(r_aspect) <- cos(runif(ncell(r_veg), 0, 2 * pi) - 315 * pi / 180)
values(r_wind) <- runif(ncell(r_veg), 0, 2 * pi)
values(r_pp) <- runif(ncell(r_veg), 0, 1)
values(r_temp) <- runif(ncell(r_veg), 0, 1)

r_predictors <- c(r_veg, r_aspect, r_pp, r_temp, r_elev, r_wind)
names(r_predictors) <- names(coefs)[-1]

# values(r_predictors)[10, ]

# vector of distances and angles between a cell and its neighbours
# (used to compute slope and wind effects)

distances <- rep(res(r_predictors)[1], 8) # sides
distances[c(1, 3, 6, 8)] <- res(r_predictors)[1] * sqrt(2)

angles_matrix <- matrix(
  c(315, 0, 45,
    270, NA, 90,
    225, 180, 135),
  3, 3, byrow = TRUE) * pi / 180
angles <- as.numeric(angles_matrix %>% t) %>% na.omit # in radians

# turn predictors into matrix
predmat <- values(r_predictors)
head(predmat)

# create fire raster
burn <- r_veg
# fill raster
# 0: unburned
# 1: burning
# 2: burned
# 3: not flammable
values(burn) <- 0 #rpois(ncell(burn), lambda = 1)

# ignition point:
#ig_location <- sample(1:ncell(burn), size = 1, replace = FALSE) #floor(ncell(burn) / 2 - 2) # cell number
ig_location <- cellFromRowCol(burn, nrow(burn) / 2, ncol(burn) / 2)
# set as burning
values(burn)[ig_location:(ig_location+5)] <- 1
# plot(r_predictors)

# initialize burning_cells cell id
burning_cells <- which(values(burn) == 1) # cell number id


# Adjacency-related functions ---------------------------------------------

# define queen neighbours
moves <- matrix(c(-1,-1,-1,  0,0,  1,1,1,
                  -1, 0, 1, -1,1, -1,0,1),
                nrow = 2, byrow = TRUE)
# in the eight-neighbour setting (1 step queen moves), to arrive at neighbours
# we can make eight moves of rows and columns. moves indicates the row movements
# (moves[1, ]) and the column movements (moves[2, ]) to get each neighbour.
# neighbours are ordered row-wise:
# (1, 2, 3,
#  4, NA, 5,
#  6, 7, 8)

# row-col-cell functions

# functions to translate cell id to row-col (analogue of terra::rowColFromCell)
# Warning: there are no range restrictions, so use carefully.

cell_to_rowcol <- function(cell, n_rowcol) {
  # n_rowcol[1] = rows; n_rowcol[2] = cols.
  row_position <- ceiling(cell / n_rowcol[2])
  col_position <- n_rowcol[2] - (row_position * n_rowcol[2] - cell)
  return(rbind(row_position, col_position))
}

# rowcol has to be a matrix (rowcol[1, ] = rows, rowcol[2, ] = columns)
rowcol_to_cell <- function(rowcol, n_rowcol) {
  return((rowcol[1, ] - 1) * n_rowcol[2] + rowcol[2, ])
}

# Adjacent  ----------------------------------------------------------------

adjacent_r <- function(cells, n_rowcol) {

  # (data to test)
  # cells = rep(burning_cells, 3)
  # cells = burning_cells
  # n_rowcol = rc

  # get row and col from cell id
  row_col <- cell_to_rowcol(cells, n_rowcol)

  # neighbours row_col
  neigh_rc <- array(NA, dim = c(2, 8, length(cells)))
  for(i in 1:length(cells)) {
    neigh_rc[, , i] <- row_col[, i] + moves # neighbours row-column pairs
  }

  # Get position of values out of range
  valid_id <- matrix(0, length(cells), 8)
  for(c in 1:length(cells)) {
    for(i in 1:8) {
      if(neigh_rc[1, i, c] > 0 & neigh_rc[1, i, c] <= n_rowcol[1] &
         neigh_rc[2, i, c] > 0 & neigh_rc[2, i, c] <= n_rowcol[2]) {
        valid_id[c, i] <- 1
      }
    }
  }

  # get cell id
  neigh_cell <- matrix(NA, length(cells), 8)
  for(c in 1:length(cells)) {
    for(i in 1:8) {
      if(valid_id[c, i] == 1) {
        tmp <- matrix(neigh_rc[, i, c], 2, 1)
        neigh_cell[c, i] <- rowcol_to_cell(tmp, n_rowcol)
      } else {
        neigh_cell[c, i] <- NA
      }
    }
  }

  return(neigh_cell)
}


# Test adjacency related functions ----------------------------------------

# sourceCpp("Spread function/spread_functions_2022-11-29.cpp")

rowcol <- c(5, 6)
rtest <- rast(ncol = rowcol[2], nrow = rowcol[1],
              xmin = -1000, xmax = 1000,
              ymin = -1000, ymax = 1000)


# cell to rowcol

matrix(1:prod(rowcol), byrow = TRUE, ncol = rowcol[2])
matrix(1:prod(rowcol), byrow = TRUE, ncol = rowcol[2]) - 1

cell <- sample(1:prod(rowcol), size = 1)
print(paste("cell: ", cell))

cell_to_rowcol(cell, rowcol)
cell_to_rowcol_cpp(cell, rowcol)
terra::rowColFromCell(rtest, cell) %>% t
cell_to_rowcol_cpp0(cell - 1, rowcol)

# rowcol to cell

rowcol <- c(5, 6)
rast <- rast()


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
terra::cellFromRowCol(rtest, row_col_id[1], row_col_id[2]) # terra
rowcol_to_cell_cpp0(row_col_id - 1, rowcol)

# adjacent

cell <- sample(1:prod(rowcol), size = 1)
print(paste("cell: ", cell))

matrix(1:prod(rowcol), byrow = TRUE, ncol = rowcol[2])
matrix(1:prod(rowcol), byrow = TRUE, ncol = rowcol[2]) - 1

adjacent_r(cell, rowcol)
adjacent_cpp(cell, rowcol, moves)
terra::adjacent(rtest, cell, directions = "queen")
adjacent_cpp0(cell - 1, rowcol, moves)


# Spread around function test -----------------------------------------------


spread_around_r <- function(data_burning,
                            data_neighbours,
                            coef = coefs,
                            positions = 1:8,
                            distances = distances,
                            angles = angles,
                            upper_limit = 1) {

  # recompute wind and elev columns (will be slope and wind effects)
  data_neighbours[, "elev"] <- sin(atan((data_neighbours[, "elev"] - data_burning[1, "elev"]) / distances[positions]))
  # CAREFUL: data_burning was a df, so we need to select its unic
  # and first row: [1, "elev"]

  data_neighbours[, "wind"] <- cos(angles[positions] - data_burning[1, "wind"])
  # the same here, watch the [1, "wind"]

  # careful with upper limit
  probs <- plogis(coefs[1] + as.matrix(data_neighbours) %*% coefs[2:length(coefs)]) * upper_limit
  probs <- as.numeric(probs)

  # set.seed(seed)
  burn <- rbinom(length(probs), size = 1, prob = probs)

  # return both the probability and the burn to check with the cpp function
  result <- cbind(probs, burn)
  colnames(result) <- c("probs", "burn")
  # return(burn)
  return(result)
}

# Test (use burn raster, created above)

# get cell id from neighbours
rc = c(nrow(burn), ncol(burn))
burning_cell <- sample(1:prod(rc), size = 1)

neighbours_matrix <- adjacent_r(cells = burning_cell, n_rowcol = rc)
notna_neighs <- which(!is.na(neighbours_matrix[1, ]))

neighbours_matrix <- neighbours_matrix[1, notna_neighs, drop = F]
posits <- (1:8)[notna_neighs]

db <- r_predictors[burning_cell] # data burning
dn <- r_predictors[neighbours_matrix[1, ]] # data neighbours

# For cpp we have to provide arguments in its simplest form, so turn
# neighbours data into matrix:
m <- matrix(NA, nrow = nrow(dn), ncol = ncol(dn))
for(i in 1:ncol(m)) m[, i] <- dn[, i] %>% as.numeric
# simulate a few burns
m
# sourceCpp("Spread function/spread_functions_2022-11-29.cpp")

spread_around_prob_cpp(data_burning = as.numeric(db),
                  data_neighbours = m,
                  coef = unname(coefs),
                  positions = posits - 1, # position of valid neighbours
                  wind_column = wind_column - 1,
                  elev_column = elev_column - 1,
                  distances = distances,
                  angles = as.numeric(angles), # in radians
                  upper_limit = 1)
m
# la función de cpp modifica la matrix de data_neighbours, entonces cuando
# la vuelvo a invocar usa otros datos y da cualquiera.
# arreglar eso (duplicar no funcionó)

spread_around_r(data_burning = db,
                data_neighbours = dn, # with column names!
                coef = unname(coefs),
                positions = posits, # position of valid neighbours
                distances = distances,
                angles = as.numeric(angles), # in radians
                upper_limit = 1)[, "probs"]
# probabilities are OK

# test stochastic simulation
seed <- round(runif(1, 1, 1e6))
print(paste("seed = ", seed))

set.seed(seed)
sim_r <- replicate(10,
spread_around_r(data_burning = db,
                data_neighbours = dn, # with column names!
                coef = unname(coefs),
                positions = posits, # position of valid neighbours
                distances = distances,
                angles = as.numeric(angles), # in radians
                upper_limit = 1)[, "burn"]
)

set.seed(seed)
sim_cpp <- replicate(10,
spread_around_cpp(data_burning = as.numeric(db),
                  data_neighbours = m,
                  coef = unname(coefs),
                  positions = posits - 1, # position of valid neighbours
                  wind_column = wind_column - 1,
                  elev_column = elev_column - 1,
                  distances = distances,
                  angles = as.numeric(angles), # in radians
                  upper_limit = 1)
)

all.equal(sim_r, sim_cpp)
# perfect.


# simulate_fire function test -----------------------------------------------

# create fire raster
r <- r_veg
# fill raster
# 0: unburned
# 1: burning
# 2: burned
# 3: not flammable
values(r) <- 0 #rpois(ncell(r), lambda = 1)

# ignition point:
#ig_location <- sample(1:ncell(r), size = 1, replace = FALSE) #floor(ncell(r) / 2 - 2) # cell number
ig_location <- cellFromRowCol(r, nrow(r) / 2, ncol(r) / 2)
# set as burning
values(r)[ig_location] <- 1
# plot(r_predictors)

# burned vector:
burned <- numeric(nrow(predmat))
rc <- c(nrow(r), ncol(r))

# function to spread

simulate_fire_r <- function(landscape = r_predictors,
                            burnable = rep(1, nrow(r_predictors)),
                            ignition_cells = burning_cells,
                            n_rowcol = c(nrow(r_predictors), ncol(r_predictors)),
                            coef = coefs,
                            moves = moves,
                            wind_column = wind_column,
                            elev_column = elev_column,
                            distances = distances,
                            angles = as.numeric(angles),
                            upper_limit = 1.0) {

  ### definitions for testing
  landscape = r_predictors
  burnable = rep(1, nrow(r_predictors))
  ignition_cells = burning_cells
  n_rowcol = c(nrow(r_predictors), ncol(r_predictors))
  coef = coefs
  moves = moves
  wind_column = wind_column
  elev_column = elev_column
  distances = distances
  angles = as.numeric(angles)
  upper_limit = 1.0
  ###

  n_row <- n_rowcol[1]
  n_col <- n_rowcol[2]
  n_cell <- n_row * n_col

  # // Create burn layer, which will be exported.
  burn = rep(0, n_cell)
  # // Fill the burn vector with 3 in the non-burnable pixels
  burn[burnable == 0] <- 3
  # // Set to 1 the ignition point
  burn[ignition_cells] <- 1

  # // Get burning cells (initialized at the ignition point)
  burning_cells <- ignition_cells

  # burn is the vector to hold the burned and burnable state;
  # predictors_matrix si the raster of predictors into matrix form;
  # n_rowcol is the number of rows and columns of the raster;

  # spread
  j = 1
  while(length(burning_cells) > 0) {
    # print(paste("cycle ", j, sep = ""))

    # get cell id from neighbours
    neighbours_matrix <- terra::adjacent(landscape,
                                         cells = burning_cells,
                                         directions = "queen")
    #colnames(neighbours_matrix) <- 1:8

    # spread from burning pixels
    for(b in 1:length(burning_cells)) {
      # b = 1
      # print(paste("burning cell ", b, sep = ""))

      #b = 1
      # get neighbours available to burn (cell ids), while getting rid of
      # non-existent neighbours

      filter <- !is.na(neighbours_matrix[b, ]) &
                burn[neighbours_matrix[b, ]] == 0
      cols_use <- which(filter)
      neighbours <- neighbours_matrix[b, cols_use]

      # Subset required data from landscape
      data_burning <- values(landscape)[burning_cells[b], , drop = FALSE]
      data_neighbours <- values(landscape)[neighbours, , drop = FALSE]

      # simulate spread
      if(length(neighbours) > 0) {

        burn[neighbours] <- spread_around_r(data_burning = data_burning,
                                            data_neighbours = data_neighbours, # with column names!
                                            coef = coef,
                                            positions = cols_use, # position of valid neighbours
                                            distances = distances,
                                            angles = angles, # in radians
                                            upper_limit = upper_limit)[, "burn"]
      }

    } # end loop over burning pixels

    # update cycle step
    j <- j + 1

    # update: burning to burned
    burn[burning_cells] <- 2

    # update burning_cells
    burning_cells <- which(burn == 1)
  }

  return(burn)
}



# Test
# sourceCpp("Spread function/spread_functions_2022-11-29.cpp")
seed = runif(1, 1, 20000) %>% round

# par(mfrow = c(1, 2))

# set.seed(seed)
# system.time(
# burn_sim_r <- simulate_fire_r(landscape = r_predictors,
#                               burnable = rep(1, nrow(r_predictors)),
#                             ignition_cells = burning_cells,
#                             n_rowcol = c(nrow(r_predictors), ncol(r_predictors)),
#                             coef = coefs,
#                             moves = moves,
#                             wind_column = wind_column,
#                             elev_column = elev_column,
#                             distances = distances,
#                             angles = as.numeric(angles),
#                             upper_limit = 1.0)
# )
# burn_sim_rast_r <- r_veg
# values(burn_sim_rast_r) <- burn_sim_r
# plot(burn_sim_rast_r, main = "R")

# sourceCpp("spread_functions.cpp")
# set.seed(seed)
# system.time(
# burn_sim_cpp <- simulate_fire_cpp(landscape = values(r_predictors),
#                                   burnable = rep(1, ncell(r_predictors)),
#                                   ignition_cells = burning_cells - 1,
#                                   n_rowcol = c(nrow(r_predictors), ncol(r_predictors)),
#                                   coef = coefs,
#                                   moves = moves,
#                                   wind_column = wind_column - 1,
#                                   elev_column = elev_column - 1,
#                                   distances = distances,
#                                   angles = as.numeric(angles),
#                                   upper_limit = 1.0)
# )

# burn_sim_rast_cpp <- r_veg
# values(burn_sim_rast_cpp) <- burn_sim_cpp
# plot(burn_sim_rast_cpp, main = "C++")
# burn_sim_cpp
# par(mfrow = c(1, 1))

# con size = 25, cpp tarda solo un 1 % de lo de R.
# con size = 200, cpp tarda un 0.12 % de R. (833 : 1 relation)
# con size = 100, cpp tarda un 0.42 % de R. (238 : 1 relation)
#                                            258 : 1
# con size = 1600 (~ 50 km de lado), cpp tarda un ___ % de R.
#   más alla de lo que tarde R, lleva al menos 9 h y sigue corriendo
#   (un solo fuego.) Tiene p re alta, o sea que quema casi todo el paisaje, pero
#   igual no hay need de que tarde tanto. Y la RAM está por 4 Gb de uso,
#   pero R avisa que usa ~ 1.4 por ahora.

# Tiempos cpp (s), quemando todo, paisaje cuadrado de size == ncol == nrow
# 800,  11.477
# 1600, 70.230  (es menos que al cuadrado)


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
    moves = moves,
    wind_column = wind_column - 1,
    elev_column = elev_column - 1,
    distances = distances,
    angles = as.numeric(angles),
    upper_limit = 1.0)

  return(burn_sim_cpp_notadj)
}

run_cpp_adj <- function() {
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