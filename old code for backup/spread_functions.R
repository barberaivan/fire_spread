# Functions to spread fire, written in R to check results with C++.

# The fire spread model is a cellular automata that spreads towards the 8-
# neighbouring pixels.

# For a deeper explanation of functions, see <spread_functions.cpp>.

library(terra)

# Define a few constants ------------------------------------------------

# Elevation data to standardize distance between pixels
elevation_mean <- 1163.3
elevation_sd <- 399.5

# Distance between pixels, in m / elevation_sd.
# 30 m is the landscape resolution. 
distances_raw <- rep(30, 8)
distances_raw[c(1, 3, 6, 7)] <- distances_raw[c(1, 3, 6, 7)] * sqrt(2)
distances <- distances_raw / elevation_sd

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

# Terrain coefficients names
northing_id <- 1
elev_id     <- 2 
windir_id   <- 3
slope_id    <- 4

# Vegetation coefficients names
shrubland_id <- 1
subalpine_id <- 2 
wet_id       <- 3
dry_a_id     <- 4 # araucaria
dry_b_id     <- 5 # cypress
steppe_id    <- 6


# Raster-matrix-array conversion functions --------------------------------

# Terrain 3D array from landscape SpatRaster. The first layer is vegetation, 
# so it's ignored
land_cube <- function(x) {
  v <- values(x)
  a <- array(NA, dim = c(nrow(x), ncol(x), nlyr(x)-1),
             dimnames = list(row = NULL, col = NULL, layer = names(x)[-1]))
  for(l in 2:nlyr(x)) a[, , l-1] <- matrix(v[, l], nrow(x), ncol(x), 
                                           byrow = TRUE)
  return(a)
}

# Function to turn matrix into SpatRaster (for plotting)
rast_from_mat <- function(m, fill_raster) { # fill_raster is a SpatRaster from terra
  mt <- t(m)
  for(i in 1:nrow(m)) mt[, i] <- m[i, ]
  r <- fill_raster[[1]]
  values(r) <- as.numeric(mt)
  return(r)
}

# Spread functions --------------------------------------------------------

#' @title spread_onepix_prob_cpp
#' @description Calculates the probability of a cell spreading fire to another.
#' @return float [0, 1] indicating the probability.
#'
#' @param int vegetation_type: vegetation type of the target cell.
#' @param arma::frowvec terrain_burning: terrain data from burning cell. 
#' @param arma::frowvec terrain_neighbour: terrain data from target neighbour.
#' @param arma::frowvec coef_veg: intercepts in logistic regression associated
#'   to each vegetation type. To aggregate vegetation types (e.g., dry_forest =
#'   dry_forest_a or dry_forest_b), the same values are assigned to the 
#'   aggregated classes.
#' @param arma::frowvec coef_terrain: slopes in logistic regression related
#'   multiplying the terrain predictors. 
#' @param int position: relative position of the neighbour in relation to the
#' burning cell. The eight neighbours are labelled from 0 to 7 beginning from
#' the upper-left one (by row):
#'   0 1 2
#'   3   4
#'   5 6 7.
#'   This is necessary to compute the slope and wind effects, as they
#'   depend on the angle and distance between burning and target pixels.
#' @param float upper_limit: upper limit for spread probability (setting to
#'   1 makes absurdly large fires; 0.5 is preferred).

spread_onepix_r <- function(
    vegetation_type, # starts at 0
    terrain_burning,
    terrain_neighbour,
    coef_veg,
    coef_terrain,
    position,
    upper_limit = 1.0
  ) {

  # wind term
  wind_term = cos(angles[position] - terrain_burning[windir_id])
  
  # slope term (from elevation and distance)
  slope_term = sin(atan(
    (terrain_neighbour[elev_id] - terrain_burning[elev_id]) / distances[position]
  ))
  
  # compute linear predictor
  linpred <- 
    coef_veg[vegetation_type + 1] + # vegetation intercept, add 1 because 
                                    # vegetation classes start at 0
    coef_terrain[northing_id] * terrain_neighbour[northing_id] +
    coef_terrain[elev_id]     * terrain_neighbour[elev_id] +
    coef_terrain[windir_id]   * wind_term +
    coef_terrain[slope_id]    * slope_term;
  
  # burn probability
  probs <- plogis(linpred) * upper_limit

  burn <- rbinom(1, size = 1, prob = probs)

  # return both the probability and the burn to check with the cpp function
  result <- c(probs, burn)
  names(result) <- c("probs", "burn")

  return(result)
}

# .........................................................................

#' @title simulate_fire_r
#' @description function to simulate a fire spread given the landscape,
#'   model coefficients and ignition points. 
#' @return burned_res: struct containing the following objects:
#'   IntegerMatrix burned_bin, a binary matrix indicating the burned pixels
#'   IntegerMatrix burned_ids, a matrix with a column by burned pixel, 
#'     indicating its row (row1) and column (row2) in the landscape,
#'   int end, the number of burned pixels.

#' @param IntegerMatrix vegetation: integer matrix representing the vegetation 
#'   type. 99 is non-burnable, and valid values are {0, ..., n_veg_types - 1}.
#' @param SpatRaster[terra] terrain: raster with terrain data, where each layer
#'   is a predictor: {northing, elevation, wind direction}. Slope is absent
#'   because it's directional, so it's computed during the simulation.
#'   It is meant to be turned into a 3D-array [rows, cols, layers], but having 
#'   the terra object is convenient for plotting the progress.
#' @param int n_veg_types: integer indicating the number of vegetation types
#'   considered by the model (not just those present in the focal landscape).
#'   used to read properly the vegetation- and non-vegetation parameters in 
#'   coef.  
#'   
#' @param IntegerMatrix ignition_cells(2, burning_cells): row and column id for
#'   the cell(s) where the fire begun. First row has the row_id, second row has
#'   the col_id.
#' @param arma::frowvec coef: parameters in logistic regression to compute the
#'   spread probability as a function of covariates.
#' @param float upper_limit: upper limit for spread probability (setting to
#'   1 makes absurdly large fires).
                                                                                                                                    #'   testing.)
#' @param bool plot_animation: whether to plot the fire progress while running
#'   or not (set to FALSE by default).

simulate_fire_r <- function(
    vegetation,
    terrain,
    ignition_cells,
    coef,
    n_veg_types = 6,
    upper_limit = 1.0,
    plot_animation = FALSE
  ) {
  
  # separate vegetation and terrain coefficients
  coef_veg = coef[1:n_veg_types]
  coef_terrain = coef[(n_veg_types + 1) : length(coef)];
  
  # define landscape dimensions
  n_row <- nrow(vegetation)
  n_col <- ncol(vegetation)
  n_cell <- n_row * n_col
  n_layers <- nlyr(terrain)

  # turn terrain into numeric array
  terrain_arr <- array(NA, dim = c(n_row, n_col, n_layers))
  terrain_values <- terra::values(terrain) # get values in matrix form
  for(l in 1:n_layers) {
    terrain_arr[, , l] <- matrix(terrain_values[, l], n_row, n_col,
                                   byrow = TRUE) # byrow because terra provides
    # the values this way.
  }

  # Create burn layer, which will be exported.
  burned_bin = matrix(0, n_row, n_col)

  # Make burning_ids matrix
  burning_ids <- matrix(NA, 2, n_cell)

  # Initialize burning_ids and burned_bin
  for(i in 1:ncol(ignition_cells)) {
    burning_ids[1, i] <- ignition_cells[1, i]
    burning_ids[2, i] <- ignition_cells[2, i]

    burned_bin[ignition_cells[1, i], ignition_cells[2, i]] <- 1
  }

  # positions from where to read the ids of the currently burning cells
  start <- 1
  end <- ncol(ignition_cells)

  # Fire raster for plotting
  burn_raster <- terrain[[1]]
  values(burn_raster) <- 0

  # get burning cells ids
  burning_cells <- cellFromRowCol(burn_raster,
                                  row = burning_ids[1, start:end],
                                  col = burning_ids[2, start:end])

  values(burn_raster)[burning_cells] <- 1
  if (plot_animation) {
    plot(burn_raster, col = c("green", "red"), main = "step 0")
  }

  # spread
  j <- 1
  burning_size <- length(burning_cells)

  while(burning_size > 0) {
    # Loop over all the burning cells to burn their neighbours. Use end_forward
    # to update the last position in burning_ids within this loop, without
    # compromising the loop's integrity.
    end_forward <- end

    # Loop over burning cells in the cycle

    # b is going to keep the position in burning_ids that have to be evaluated
    # in this burn cycle

    # spread from burning pixels
    for(b in start:end) {
      # Get burning_cells' data
      terrain_burning <- terrain_arr[burning_ids[1, b], burning_ids[2, b], ];

      # get neighbours (adjacent computation here)
      neighbours <- burning_ids[, b] + moves

      # Loop over neighbours of the focal burning cell
      for(n in 1:8) {
        
        # Is the cell in range?
        out_of_range <- (
          (neighbours[1, n] < 1) | (neighbours[1, n] > n_row) | # check rows
          (neighbours[2, n] < 1) | (neighbours[2, n] > n_col)   # check cols
        )
        if(out_of_range) next # (jumps to next iteration if TRUE)
        
        # Get vegetation class to know whether it's burnable
        veg_target <- vegetation[neighbours[1, n], neighbours[2, n]]
        
        # Is the cell burnable?
        burnable_cell <- 
          (burned_bin[neighbours[1, n], neighbours[2, n]] == 0) &
          (veg_target < 99) # 99 is non-burnable
        if(!burnable_cell) next

        # obtain data from the neighbour
        terrain_neighbour = terrain_arr[neighbours[1, n], neighbours[2, n], ];

        # simulate fire
        burn <- spread_onepix_r(
          veg_target,
          terrain_burning,
          terrain_neighbour,
          coef_veg,
          coef_terrain,
          n, # neighbour position (in 1:8)
          upper_limit
        )["burn"] # because it returns also the probability
        
        if(burn == 0) next

        # If burned,
        # store id of recently burned cell and
        # set 1 to burned_bin
        # (but advance end_forward first)
        end_forward <- end_forward + 1
        burning_ids[1, end_forward] = neighbours[1, n]
        burning_ids[2, end_forward] = neighbours[2, n]
        burned_bin[neighbours[1, n], neighbours[2, n]] <- 1
      } # end loop over neighbours of burning cell b

    } # end loop over burning cells from this cycle (j)

    # update start and end
    start <- end + 1
    end <- end_forward
    burning_size <- end - start + 1

    # update: burning to burned
    values(burn_raster)[burning_cells] <- 2

    if(burning_size > 0) {
      # update burning_cells (this correspond to the next cycle)
      burning_cells <- cellFromRowCol(burn_raster,
                                      row = burning_ids[1, start:end],
                                      col = burning_ids[2, start:end])
      values(burn_raster)[burning_cells] <- 1
    }

    if (plot_animation) {
      plot(burn_raster, col = c("green", "red", "black"), main = paste("step", j))
    }

    # update cycle step
    j <- j + 1
  }

  return(burned_bin)
}

# .........................................................................
# The same function but deterministic, to test if the discrepancy between R and
# cpp is caused by seed problems

simulate_fire_deterministic_r <- function(
    vegetation,
    terrain,
    ignition_cells,
    coef,
    n_veg_types = 6,
    upper_limit = 1.0,
    plot_animation = FALSE
) {
  
  # separate vegetation and terrain coefficients
  coef_veg = coef[1:n_veg_types]
  coef_terrain = coef[(n_veg_types + 1) : length(coef)];
  
  # define landscape dimensions
  n_row <- nrow(vegetation)
  n_col <- ncol(vegetation)
  n_cell <- n_row * n_col
  n_layers <- nlyr(terrain)
  
  # turn terrain into numeric array
  terrain_arr <- array(NA, dim = c(n_row, n_col, n_layers))
  terrain_values <- terra::values(terrain) # get values in matrix form
  for(l in 1:n_layers) {
    terrain_arr[, , l] <- matrix(terrain_values[, l], n_row, n_col,
                                 byrow = TRUE) # byrow because terra provides
    # the values this way.
  }
  
  # Create burn layer, which will be exported.
  burned_bin = matrix(0, n_row, n_col)
  
  # Make burning_ids matrix
  burning_ids <- matrix(NA, 2, n_cell)
  
  # Initialize burning_ids and burned_bin
  for(i in 1:ncol(ignition_cells)) {
    burning_ids[1, i] <- ignition_cells[1, i]
    burning_ids[2, i] <- ignition_cells[2, i]
    
    burned_bin[ignition_cells[1, i], ignition_cells[2, i]] <- 1
  }
  
  # positions from where to read the ids of the currently burning cells
  start <- 1
  end <- ncol(ignition_cells)
  
  # Fire raster for plotting
  burn_raster <- terrain[[1]]
  values(burn_raster) <- 0
  
  # get burning cells ids
  burning_cells <- cellFromRowCol(burn_raster,
                                  row = burning_ids[1, start:end],
                                  col = burning_ids[2, start:end])
  
  values(burn_raster)[burning_cells] <- 1
  if (plot_animation) {
    plot(burn_raster, col = c("green", "red"), main = "step 0")
  }
  
  # spread
  j <- 1
  burning_size <- length(burning_cells)
  
  while(burning_size > 0) {
    # Loop over all the burning cells to burn their neighbours. Use end_forward
    # to update the last position in burning_ids within this loop, without
    # compromising the loop's integrity.
    end_forward <- end
    
    # Loop over burning cells in the cycle
    
    # b is going to keep the position in burning_ids that have to be evaluated
    # in this burn cycle
    
    # spread from burning pixels
    for(b in start:end) {
      # Get burning_cells' data
      terrain_burning <- terrain_arr[burning_ids[1, b], burning_ids[2, b], ];
      
      # get neighbours (adjacent computation here)
      neighbours <- burning_ids[, b] + moves
      
      # Loop over neighbours of the focal burning cell
      for(n in 1:8) {
        
        # Is the cell in range?
        out_of_range <- (
          (neighbours[1, n] < 1) | (neighbours[1, n] > n_row) | # check rows
          (neighbours[2, n] < 1) | (neighbours[2, n] > n_col)   # check cols
        )
        if(out_of_range) next # (jumps to next iteration if TRUE)
        
        # Get vegetation class to know whether it's burnable
        veg_target <- vegetation[neighbours[1, n], neighbours[2, n]]
        
        # Is the cell burnable?
        burnable_cell <- 
          (burned_bin[neighbours[1, n], neighbours[2, n]] == 0) &
          (veg_target < 99) # 99 is non-burnable
        if(!burnable_cell) next
        
        # obtain data from the neighbour
        terrain_neighbour = terrain_arr[neighbours[1, n], neighbours[2, n], ];
        
        # simulate fire
        burn <- spread_onepix_r(
          veg_target,
          terrain_burning,
          terrain_neighbour,
          coef_veg,
          coef_terrain,
          n, # neighbour position (in 1:8)
          upper_limit
        )#["burn"] # because it returns also the probability
        
        # make deterministic!
        if(burn["probs"] < 0.5000000000) next
        
        # If burned,
        # store id of recently burned cell and
        # set 1 to burned_bin
        # (but advance end_forward first)
        end_forward <- end_forward + 1
        burning_ids[1, end_forward] = neighbours[1, n]
        burning_ids[2, end_forward] = neighbours[2, n]
        burned_bin[neighbours[1, n], neighbours[2, n]] <- 1
      } # end loop over neighbours of burning cell b
      
    } # end loop over burning cells from this cycle (j)
    
    # update start and end
    start <- end + 1
    end <- end_forward
    burning_size <- end - start + 1
    
    # update: burning to burned
    values(burn_raster)[burning_cells] <- 2
    
    if(burning_size > 0) {
      # update burning_cells (this correspond to the next cycle)
      burning_cells <- cellFromRowCol(burn_raster,
                                      row = burning_ids[1, start:end],
                                      col = burning_ids[2, start:end])
      values(burn_raster)[burning_cells] <- 1
    }
    
    if (plot_animation) {
      plot(burn_raster, col = c("green", "red", "black"), main = paste("step", j))
    }
    
    # update cycle step
    j <- j + 1
  }
  
  return(burned_bin)
}