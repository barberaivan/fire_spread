# Functions to spread fire, written in R to check results with C++.

# The fire spread model is a cellular automata that spreads towards the 8-
# neighbouring pixels.

library(terra)

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


# Spread functions --------------------------------------------------------

#' @title spread_onepix_r
#' @description Spreads fire from one burning cell to one of its neighbours.
#'   Most of the work implies computing the burn probability given the
#'   covariates. This is the function that should be edited when alternative
#'   spread functions are considered.
#' @return c("probs", "burn") with the probability and the burn

#' @param NumericVector data_burning: environmental data from burning cell
#'   (the most important here are wind direction and elevation).
#' @param NumericMatrix data_neighbours: environmental data from target
#'   neighbours with a column by landscape layer.
#' @param NumericVector coef: parameters in logistic regression to compute the
#'   spread probability as a function of covariates.
#' @param IntegerVector position: relative position of the neighbour in
#'   relation to the burning cell. The eight neighbours are labelled from 0 to
#'   7 beggining from the upper-left one (by row):
#'   0 1 2
#'   3   4
#'   5 6 7.
#'   This is necessary to compute the elevation and wind effects, as they
#'   depend on the angle and distance between burning and target pixels.
#' @param int wind_column: column in the data (landscape) with wind value.
#' @param int elev_column: column in the data (landscape) with elevation value.
#'   Wind and elevation columns must be the last 2.
#' @param NumericVector distances: distances (m) between burning and target cells,
#'   in the same order as positions. Used to compute the elevation effect.
#'   This vector depends on the neighbourhood design and on the pixel scale.
#'   If unchanged, it's always the same.
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


# The same function but to be used for a single pixel
spread_onepix_r <- function(data_burning,
                            data_neighbour,
                            coef,
                            position,
                            distances,
                            elev_column,
                            wind_column,
                            upper_limit = 1) {

  # compute wind, elevation and slope terms
  slope_term <- sin(atan(
    (data_neighbour[elev_column] - data_burning[elev_column]) / distances[position]
    ))

  # standardize elevation
  data_neighbour[elev_column] <- (data_neighbour[elev_column] - elevation_mean) / elevation_sd

  #modify wind term
  data_neighbour[wind_column] <- cos(angles[position] - data_burning[wind_column])

  # bind slope term
  data_neighbour <- c(data_neighbour, slope_term)

  # compute linear predictor
  linpred <- as.numeric(coef[1] + t(as.matrix(data_neighbour)) %*% coef[2:length(coef)])

  # compute probability
  probs <- plogis(linpred) * upper_limit

  burn <- rbinom(1, size = 1, prob = probs)

  # return both the probability and the burn to check with the cpp function
  result <- c(probs, burn)
  names(result) <- c("probs", "burn")

  return(result)
}

# function to simulate a fire spread given the landscape,
# model coefficients and ignition points.

#' @title simulate_fire_mat_r
#' @description function to simulate a fire spread given the landscape,
#'   model coefficients and ignition points. This uses a matrix representation
#'   of the landscape, using row-col ids for cells, instead of cell number.
#'   In addition, this function plots the fire progress.
#' @return IntegerMatrix(n_row, n_col): burned layer, a matrix with zeroes
#'   (unburned) and ones (burned).

#' @param SpatRaster[terra] landscape: raster with environmental data from the
#'   whole landscape. It is meant to be turned into an array
#'   [rows, cols, layers], but having the terra object is convenient for
#'   plotting.
#' @param IntegerMatrix ignition_cell: row-col id for the cell(s) where the fire
#'   begun, with a cell by column.
#' @param IntegerMatrix burnable: matrix indicating if each pixel is burnable (1)
#'   or not (0).
#' @param NumericVector coef: parameters in logistic regression to compute the
#'   spread probability as a function of covariates.
#' @param int wind_layer: layer in the data (landscape) with wind values.
#' @param int elev_layer: layer in the data (landscape) with elevation values.
#'   Wind and elevation layers must be the last 2.
#' @param NumericVector distances: distances (m) between burning and target cells,
#'   in the same order as positions. Used to compute the elevation effect.
#'   This vector depends on the neighbourhood design and on the pixel scale.
#'   If unchanged, it's always the same.
#' @param double upper_limit: upper limit for spread probability (setting to
#'   1 makes absurdly large fires).

simulate_fire_mat_r <- function(
    landscape,
    burnable,
    ignition_cells,
    coef,
    wind_column,
    elev_column,
    distances,
    upper_limit = 1.0,
    plot_animation = FALSE) {

  n_row <- nrow(burnable)
  n_col <- ncol(burnable)
  n_cell <- n_row * n_col
  n_layers <- nlyr(landscape)

  # turn landscape into numeric array
  landscape_arr <- array(NA, dim = c(n_row, n_col, n_layers))
  landscape_values <- terra::values(landscape) # get values in matrix form
  for(l in 1:n_layers) {
    landscape_arr[, , l] <- matrix(landscape_values[, l], n_row, n_col,
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
  burn_raster <- landscape[[1]]
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
      data_burning <- landscape_arr[burning_ids[1, b], burning_ids[2, b], ];

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

        # Is the cell burnable?
        burnable_cell <- (burned_bin[neighbours[1, n], neighbours[2, n]] == 0) &
                         (burnable[neighbours[1, n], neighbours[2, n]] == 1)
        if(!burnable_cell) next

        # obtain data from the neighbour
        data_neighbour = landscape_arr[neighbours[1, n], neighbours[2, n], ];

        # simulate fire
        burn <- spread_onepix_r(
          data_burning,
          data_neighbour,
          coef,
          n,
          distances,
          elev_column,
          wind_column,
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

# ...........................................................................
# The same function but deterministic, to test if the discrepancy between R and
# cpp is caused by seed problems

simulate_fire_mat_deterministic_r <- function(
    landscape,
    burnable,
    ignition_cells,
    coef,
    wind_column,
    elev_column,
    distances,
    upper_limit = 1.0,
    plot_animation = FALSE) {

  n_row <- nrow(burnable)
  n_col <- ncol(burnable)
  n_cell <- n_row * n_col
  n_layers <- nlyr(landscape)

  # turn landscape into numeric array
  landscape_arr <- array(NA, dim = c(n_row, n_col, n_layers))
  landscape_values <- terra::values(landscape) # get values in matrix form
  for(l in 1:n_layers) {
    landscape_arr[, , l] <- matrix(landscape_values[, l], n_row, n_col,
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

  burn_raster <- landscape[[1]]
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
      data_burning <- landscape_arr[burning_ids[1, b], burning_ids[2, b], ];

      # get neighbours (adjacent computation here)
      neighbours <- burning_ids[, b] + moves

      # Loop over neighbours of the focal burning cell

      for(n in 1:8) {

        # Is the cell in range?
        out_of_range <- (
          (neighbours[1, n] < 1) | (neighbours[1, n] > n_row) | # check rows
          (neighbours[2, n] < 1) | (neighbours[2, n] > n_col)   # check cols
        );
        if(out_of_range) next # (jumps to next iteration if TRUE)

        # Is the cell burnable?
        burnable_cell <- (burned_bin[neighbours[1, n], neighbours[2, n]] == 0) &
                         (burnable[neighbours[1, n], neighbours[2, n]] == 1)
        if(!burnable_cell) next

        # obtain data from the neighbour
        data_neighbour = landscape_arr[neighbours[1, n], neighbours[2, n], ];

        # simulate fire
        burn <- spread_onepix_r(
          data_burning,
          data_neighbour,
          coef,
          n,
          distances,
          elev_column,
          wind_column,
          upper_limit
        )#["burn"] # because it returns also the probability

        ## make deterministic!!
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