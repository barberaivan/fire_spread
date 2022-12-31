# Functions to spread fire, written in R to check results with C++.

# The fire spread model is a cellular automata that spreads towards the 8-
# neighbouring pixels. Rasters are treated as vectors, with a cell id by pixel.
# We need functions to translate between the vectorized cell id to row-column
# and the other way around.
# In this script all functions are created for 1-starting-indexes.

# Functions to find neighbours of each target cell:
#  cell_to_rowcol
#  rowcol_to_cell
#  adjacent

# Functions to actually spread fire:
#  spread_around_cpp
#    (spreads fire towards the neighbours of one burning cell. 
#    IMPORTANT: this function holds the spread model itself, so it has to be 
#    edited if the model is changed.)
#  simulate_fire_cpp
#    (basically, a while loop running spread_around_cpp as long as there
#    are burning cells)



# Packages ----------------------------------------------------------------

library(terra)
library(tidyverse) # for pipe operator %>% 

# Define a few constants ------------------------------------------------

# Elevation metrics to standardize predictor
elevation_mean <- 1163.3
elevation_sd <- 399.5

# Angles between cells to compute wind effect. As the wind direction is 
# the direction from which the wind comes, these angles must represent where the
# fire would come from if from the neighbours we look at the central pixel.
angles <- c(
    135, 180, 225,
    90,       270,
    45,   0,  315
    ) * pi / 180   # in radians

# queen neighbours movements
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

# Adjacency-related functions ---------------------------------------------

#' @title cell_to_rowcol
#' @description Translates cell id to row and column.
#' @return IntegerMatrix(2, n_cells): matrix with row and column ids, with each
#'   column corresponding to a cell.
#' @param NumericVector cells: cell ids.
#' @param NumericVector n_rowcol: number of row and columns of the landscape.

cell_to_rowcol <- function(cell, n_rowcol) {
  row_position <- ceiling(cell / n_rowcol[2])
  col_position <- n_rowcol[2] - (row_position * n_rowcol[2] - cell)
  return(rbind(row_position, col_position))
}

#' @title rowcol_to_cell
#' @description Translates row and column id to cell id.
#' @return IntegerVector(n_cells): matrix with row and column ids, with each
#'   column corresponding to a cell.
#' @param IntegerMatrix rowcol: matrix with row and column ids, with each
#'   column corresponding to a cell.
#' @param NumericVector n_rowcol: number or row and columns of the landscape.

rowcol_to_cell <- function(rowcol, n_rowcol) {
  return((rowcol[1, ] - 1) * n_rowcol[2] + rowcol[2, ])
}


#' @title adjacent_r
#' @description Gets the cell id of the neighbours of focal (burning) cells,
#'   considering only "queen" direction following terra naming. It
#'   takes into account that neighbours can't fall outside the landscape.
#' @return IntegerMatrix(number of focal cells, 8): focal cells in rows,
#'   neighbours in columns. (8-pixels neighbourhood).
#' @param IntegerVector cells: focal (burning) cell ids.
#' @param NumericVector n_rowcol: number or row and columns of the landscape.

adjacent_r <- function(cells, n_rowcol) {
  
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



# Spread functions --------------------------------------------------------

#' @title spread_around
#' @description Spreads fire from one burning cell to its neighbours. Most of
#'   the work implies computing the burn probability given the covariates. This
#'   is the function that should be edited when alternative spread functions
#'   are considered. After the spread probability is computed, the only
#'   remaining step is to simulate a Bernoulli process in each neighbour.
#' @return IntegerVector(number of target neighbours): binary vector indicating
#'   burned (1) or not (0) for every analyzed neighbour. Note that invalid
#'   or already burned neighbours should not be included, but that subsetting
#'   is achieved by another function.

#' @param NumericVector data_burning: environmental data from burning cell
#'   (the most important here are wind direction and elevation).
#' @param NumericMatrix data_neighbours: environmental data from target
#'   neighbours with a column by landscape layer.
#' @param NumericVector coef: parameters in logistic regression to compute the
#'   spread probability as a function of covariates.
#' @param IntegerVector positions: relative position of each neighbour in
#'   relation to the burning cell. The eight neighbours are labelled from 0 to
#'   7 beggining from the upper-left one (by row):
#'   0 1 2
#'   3   4
#'   5 6 7.
#'   This is necessary to compute the elevation and wind effects, as they
#'   depend on the angle and distance between burning and target pixels.
#'   (Sometimes not all neighbours are burnable, so this vector will not
#'   always be complete.)
#' @param int wind_column: column in the data (landscape) with wind value.
#' @param int elev_column: column in the data (landscape) with elevation value.
#'   Wind and elevation columns must be the last 2.
#' @param NumericVector distances: distances (m) between burning and target cells,
#'   in the same order as positions. Used to compute the elevation effect.
#'   This vector depends on the neighbourhood design and on the pixel scale.
#'   If unchanged, it's always the same.
#' @param NumericVector angles: angles (radians) between burning and target cells,
#'   in the same order as positions. Used to compute the wind effect. This
#'   vector is always the same, unless the neighborhood is changed.
#' @param double upper_limit: upper limit for spread probability (setting to
#'   1 makes absurdly large fires).

#' The layers in data_ are:
#'   subalpine forest, {0, 1} (shrubland goes in the intercept)
#'   wet forest,       {0, 1}
#'   dry forest,       {0, 1}
#'   fwi,              (-Inf, +Inf) (standardized anomaly at pixel level)
#'   aspect            [-1, 1] # Northwestyness
#'   wind direction    [-1, 1] (windward vs leeward)
#'   elevation,        (standardized in the code)

#' The coefficients vector has a parameter by layer plus an intercept and the
#' slope effect:
#'   [Intercept] shrubland
#'   subalpine forest,
#'   wet forest,
#'   dry forest,
#'   fwi,
#'   aspect,
#'   wind,
#'   elevation,  (note slope comes after elevation)
#'   [slope],    (downhill or uphill, (-1, 1): 0 = flat, 1 = above, -1 = below)

spread_around_r <- function(data_burning,
                            data_neighbours,
                            coef,
                            positions = 1:8,
                            distances = distances,
                            upper_limit = 1) {
  
  # compute wind, elevation and slope terms 
  slope_term <- sin(atan((data_neighbours[, "elev"] - data_burning[1, "elev"]) / distances[positions]))
  
  # standardize elevation
  data_neighbours[, "elev"] <- (data_neighbours[, "elev"] - elevation_mean) / elevation_sd
  
  #modify wind term
  data_neighbours[, "wind"] <- cos(angles[positions] - data_burning[1, "wind"])
  
  # bind slope term
  data_neighbours <- cbind(data_neighbours, slope_term)
  
  # compute linear predictor
  linpred <- as.numeric(coef[1] + as.matrix(data_neighbours) %*% coef[2:length(coef)])
  
  # compute probability 
  probs <- plogis(linpred) * upper_limit
  
  
  burn <- rbinom(length(probs), size = 1, prob = probs)
  
  # return both the probability and the burn to check with the cpp function
  result <- cbind(probs, burn)
  colnames(result) <- c("probs", "burn")

  return(result)
}

# function to simulate a fire spread given the landscape,
# model coefficients and ignition points.

#' @title simulate_fire_r
#' @description function to simulate a fire spread given the landscape,
#'   model coefficients and ignition points.
#' @return IntegerVector(n_row * n_col): updated burn layer, coded as
#'   0 = not burned but burnable,
#'   1 = burning (only occurs before the function runs),
#'   2 = burned,
#'   3 = not burnable.

#' @param NumericMatrix landscape: environmental data from the whole landscape.
#' @param IntegerVector ignition_cell: id for the cell(s) where the fire begun.
#' @param IntegerVector burnable: vector indicating if each pixel is burnable (1)
#'   or not (0).
#' @param NumericVector n_rowcol: number or row and columns of the landscape.
#' @param NumericVector coef: parameters in logistic regression to compute the
#'   spread probability as a function of covariates.
#' @param IntegerVector positions: position of each neighbour in
#'   relation to the burning cell. The eight neighbours are labelled from 0 to
#'   7 beggining from the upper-left corner (by row):
#'   1 2 3
#'   4   5
#'   6 7 8.
#'   This is necessary to compute the elevation and wind effects, as they
#'   depend on the angle and distance between burning and target pixels.
#'   (Sometimes not all neighbours are burnable, so this vector will not
#'   always be complete.)
#' @param int wind_column: column in the data (landscape) with wind value.
#' @param int elev_column: column in the data (landscape) with elevation value.
#'   Wind and elevation columns must be the last 2.
#' @param NumericVector distances: distances (m) between burning and target cells,
#'   in the same order as positions. Used to compute the elevation effect.
#'   This vector depends on the neighbourhood design and on the pixel scale.
#'   If unchanged, it's always the same.
#' @param NumericVector angles: angles (radians) between burning and target cells,
#'   in the same order as positions. Used to compute the wind effect. This
#'   vector is always the same, unless the neighbourhood is changed.
#' @param double upper_limit: upper limit for spread probability (setting to
#'   1 makes absurdly large fires).

simulate_fire_r <- function(landscape,
                            burnable,
                            ignition_cells,
                            n_rowcol,
                            coef,
                            wind_column = wind_column,
                            elev_column = elev_column,
                            distances = distances,
                            upper_limit = 1.0) {
  
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
  
  # make burned binary layer to return
  burned <- rep(0, length(burn))
  burned[burn == 2] <- 1
  return(burned)
}


# The same function but plotting the intermmediate result.

simulate_fire_plot <- function(landscape,
                            burnable,
                            ignition_cells,
                            n_rowcol,
                            coef,
                            wind_column = wind_column,
                            elev_column = elev_column,
                            distances = distances,
                            upper_limit = 1.0) {
  
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
  
  # Fire raster for plotting
  burn_raster <- landscape[[1]]
  values(burn_raster) <- 0
  values(burn_raster)[burning_cells] <- 1
  # plot_colors <- data.frame(value = 0:3, color = c("green", "red", "black", "grey"))
  # plot(burn_raster, col = plot_colors) ## it does not work
  plot(burn_raster, col = c("green", "red"), main = "step 0")
  
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
                                            upper_limit = upper_limit)[, "burn"]
      }
      
    } # end loop over burning pixels
    
    # update: burning to burned
    burn[burning_cells] <- 2
    values(burn_raster)[burning_cells] <- 2
    
    # update burning_cells
    burning_cells <- which(burn == 1)
    values(burn_raster)[burning_cells] <- 1
    
    # plot(burn_raster, col = plot_colors) # not worked
    plot(burn_raster, col = c("green", "red", "black"), main = paste("step", j))

    # update cycle step
    j <- j + 1
  }
  
  # make burned binary layer to return
  burned <- rep(0, length(burn))
  burned[burn == 2] <- 1
  return(burned)
}
