# Write in base R language 2 spread functions, with recursive and 
# repetitive updating, without using terra functions.

library(bench)
library(terra)
library(tidyverse)


# Create data -------------------------------------------------------------


# model coefficients (includes intercept)
coefs <- c(0.5, -1.5, 0.5, 0.5, 0.5, 0.5, 0.5)
names(coefs) <- c("Intercept", "veg", "elev", "aspect", "wind", "pp", "temp")

# test wind effect
# coefs <- c(0.5, 0, 0, 0, 3, 0, 0)
# test slope effect
# coefs <- c(0.5, 0, 10, 0, 0, 0, 0)

# predictors raster

size <- 100

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

r_predictors <- c(r_veg, r_elev, r_aspect, r_wind, r_pp, r_temp)
names(r_predictors) <- names(coefs)[-1]

# values(r_predictors)[10, ]

# vector of distances and angles between a cell and its neighbours
# (used to compute slope and wind effects)

d <- rep(res(r_predictors)[1], 8) # sides
d[c(1, 3, 6, 8)] <- res(r_predictors)[1] * sqrt(2)

angles_matrix <- matrix(
  c(315, 0, 45,
    270, NA, 90,
    225, 180, 135), 
  3, 3, byrow = TRUE) * pi / 180
angles <- as.numeric(angles_matrix %>% t) %>% na.omit # in radians

# turn predictors into matrix
predmat <- values(r_predictors)
head(predmat)

# Spread probability function -----------------------------------------------

spread_prob <- function(data_burning, data_neighbours, upper_limit = 0.5) {
  
  # recompute wind and elev columns (will be slope and wind effects)
  data_neighbours[, "elev"] <- sin(atan((data_neighbours[, "elev"] - data_burning["elev"]) / d[cols_use]))
  data_neighbours[, "wind"] <- cos(angles[cols_use] - data_burning["wind"])
  
  # careful with upper limit
  probs <- plogis(coefs[1] + data_neighbours %*% coefs[2:length(coefs)]) * upper_limit
  return(probs)
}


# Adjacency-related functions ---------------------------------------------

# define queen neighbours
moves <- matrix(c(-1,-1,-1,0,0,1,1,1,
                  -1,0,1,-1,1,-1,0,1),
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


# test
# (rrr <- rbind(
#   sample(1:nrow(r), size = 4, replace = F),
#   sample(1:ncol(r), size = 4, replace = F)
# ))
# rowcol_to_cell(rrr,
#                c(nrow(r), ncol(r)))
# r
# 
# (sss <- sample(1:(ncol(r) * nrow(r)), size = 4, replace = F) )
# cell_to_rowcol(sss, c(nrow(r), ncol(r)))
# r


# adjacent_vector 

# This function will return the cell ids adjacent to a given set of cells,
# as a n_cell X 8 matrix (queen's neighbours). 
# The input data will be the nrow and ncol of the raster that is represented
# as a matrix


adjacent_vector <- function(cells, n_rowcol) {
  
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


# Recursive updating (test loop) --------------------------------------------

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

# initialize burning_cells cell id
burning_cells <- which(values(r) == 1) # cell number id
j = 1
plot(r, main = paste("step = ", j), col = c("green", "red"))

# burned vector:
burned <- numeric(nrow(predmat))
rc <- c(nrow(r), ncol(r))

# spread
while(length(burning_cells) > 0) {
  
  print(paste("cycle ", j, sep = ""))
  
  # update: burning to burned
  burned[burning_cells] <- 2
  
  # get cell id from neighbours
  neighbours_matrix <- adjacent_vector(cells = burning_cells, n_rowcol = rc) 
  #colnames(neighbours_matrix) <- 1:8
  
  # spread from burning pixels
  for(b in 1:length(burning_cells)) {
    print(paste("burning cell ", b, sep = ""))
    
    #b = 1
    # get neighbours available to burn (cell ids), while getting rid of 
    # non-existent neighbours
    
    filter <- !is.na(predmat[neighbours_matrix[b, ]]) &
              burned[neighbours_matrix[b, ]] == 0
    cols_use <- which(filter)
    neighbours <- neighbours_matrix[b, cols_use]
    
    # Subset required data from predictors_matrix
    data_burning <- predmat[burning_cells[b], ]
    data_neighbours <- predmat[neighbours, , drop = FALSE]
    
    # simulate spread
    if(length(neighbours) > 0) {
      
      burned[neighbours] <- rbinom(
        n = length(neighbours), size = 1,
        prob = spread_prob(data_burning, data_neighbours)
      )
    }
    
  }
  
  values(r) <- burned
  j <- j + 1
  plot(r, main = paste("step = ", j), col = c("green", "red", "black"))
  # update burning_cells
  burning_cells <- which(burned == 1)
}
#plot(r, main = paste("step = ", j), col = c("green", "red", "black"))

# Recursive updating (function) ----------------------------------------------

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
spread_redundant <- function(burning_cells_init, burned, 
                             predictors_matrix, n_rowcol) {
  
  # burning_cells_init is the cell_id for the initial burning cells
  # (ignition location);
  # burned is the vector to hold the burned and burnable state;
  # predictors_matrix si the raster of predictors into matrix form;
  # n_rowcol is the number of rows and columns of the raster;
  
  # Initialize burning cells:
  burning_cells <- burning_cells_init
  
  # spread
  while(length(burning_cells) > 0) {
    
    # print(paste("cycle ", j, sep = ""))
    
    # update: burning to burned
    burned[burning_cells] <- 2
    
    # get cell id from neighbours
    neighbours_matrix <- adjacent_vector(cells = burning_cells, n_rowcol = n_rowcol) 
    #colnames(neighbours_matrix) <- 1:8
    
    # spread from burning pixels
    for(b in 1:length(burning_cells)) {
      # print(paste("burning cell ", b, sep = ""))
      
      #b = 1
      # get neighbours available to burn (cell ids), while getting rid of 
      # non-existent neighbours
      
      filter <- !is.na(predictors_matrix[neighbours_matrix[b, ]]) &
                burned[neighbours_matrix[b, ]] == 0
      cols_use <- which(filter)
      neighbours <- neighbours_matrix[b, cols_use]
      
      # Subset required data from predictors_matrix
      data_burning <- predictors_matrix[burning_cells[b], ]
      data_neighbours <- predictors_matrix[neighbours, , drop = FALSE]
      
      # simulate spread
      if(length(neighbours) > 0) {
        
        burned[neighbours] <- rbinom(
          n = length(neighbours), size = 1,
          prob = spread_prob(data_burning, data_neighbours)
        )
      }
      
    }
    
    # update cycle step
    j <- j + 1
    # update burning_cells
    burning_cells <- which(burned == 1)
  }
  
  # when burning_cells are out, the while loop ends, so they have to become burned:
  burned[burning_cells] <- 2
  
  return(burned)
}

spread_redundant(burning_cells_init = ig_location,
                 burned = numeric(nrow(predmat)),
                 predictors_matrix = predmat,
                 n_rowcol = n_rowcol)


# Redundant updating (test loop) ------------------------------------------



# Redundant updating function ---------------------------------------------



# benchmark: redundant vs recursive ---------------------------------------
# benchmark: adjacent terra and mine --------------------------------------

bench::mark(
  adjacent_vector(cells = 1:10, n_rowcol = rc),
  adjacent(r, cells = 1:10, "queen") 
)

system.time(adjacent_vector(cells = 1:10, n_rowcol = rc))
system.time(adjacent(r, cells = 1:10, "queen"))

# el de terra is faster.


m <- matrix(NA, 10000, 10000)
m2 <- matrix(NA, 100, 100)
system.time(m[1:100, 1:100] <- 1:10000)
system.time(m2[1:100, 1:100] <- 1:10000)
# accessing a large matrix takes longer than accessing a smaller one. 
# It would make sense to use the redundant update in R, 
# but we don't know if that holds in cpp.

# Probar esas cosas pequeñas antes de pasar toda la función. 