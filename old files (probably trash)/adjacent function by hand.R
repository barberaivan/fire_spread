library(tidyverse)
library(terra)

# Adjacency function

# This function will be used to get the adjacent pixels for a burning one. It
# will have some common functionalities with the terra::adjacent function, but
# I will write it in Rcpp.

# Rasters can be thought as square matrices. Following terra, we label their cells 
# filling by row. Another way to represent them is in vector form, with the matrix
# being vectorized row-wise (not as usually done). 

# The vector form might be simpler because we reduce one dimension, but it
# requires to evaluate the neighborhood of every pixel in a (probably) more 
# convoluted way compared to the case of a really spatial matrix. 
# Anyway, in the array or matrix form, we will still need to compute the 
# neighborhood cells to follow the id of burning cells.

# We will try both approaches.

# Example data and inputs ----------------------------------------------

# create example data
n_row <- 5
n_col <- 7
n_cell <- n_col * n_row
r <- matrix(1:n_cell, n_row, n_col, byrow = TRUE)
(cell <- sample(1:n_cell, 1))

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

# row-col-cell functions ------------------------------------------------

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


# adjacent_vector ---------------------------------------------------------

# This function will return the cell ids adjacent to a given cell. 
# The input data will be the nrow and ncol of the raster that is represented
# as a matrix


adjacent_vector <- function(cell, n_rowcol) {
  
  #cell = cell; n_rowcol = ras_rc
  
  # get row and col from cell id
  row_col <- cell_to_rowcol(cell, n_rowcol) %>% as.numeric
  
  # neighbours row_col 
  neigh_rc <- row_col + moves # neighbours row-column pairs

  # Get position of values out of range
  valid_id <- numeric(8)
  for(i in 1:8) {
   if(neigh_rc[1, i] > 0 & neigh_rc[1, i] <= n_rowcol[1] &
      neigh_rc[2, i] > 0 & neigh_rc[2, i] <= n_rowcol[2]) {
     valid_id[i] <- 1
   }
  }
  # neigh_rc[1, neigh_rc[1, ] < 1 | neigh_rc[1, ] > n_row] <- NA
  # neigh_rc[2, neigh_rc[2, ] < 1 | neigh_rc[2, ] > n_col] <- NA

  # get cell id
  neigh_cell <- numeric(8)
  for(i in 1:8) {
    if(valid_id[i] == 1) {
      neigh_cell[i] <- rowcol_to_cell(neigh_rc[, i, drop = F], n_rowcol)
    } else {
      neigh_cell[i] <- NA
    }
  }

  return(neigh_cell)
}

# create raster to visualize
ras <- rast(ncol = n_col, nrow = n_row,
            xmin = 0, xmax = n_col * 30, ymin = 0, ymax = n_row * 30)
ras_rc <- c(nrow(ras), ncol(ras))

values(ras) <- 0
cell <- sample(1:(ncol(ras)*nrow(ras)), size = 1)
values(ras)[cell] <- 1
values(ras)[adjacent_vector(cell, ras_rc)] <- 2
plot(ras, col = c("black", "red", "gray"), main = paste("cell =", cell))
adjacent_vector(cell, ras_rc)

# Done


# adjacent_matrix ---------------------------------------------------------

# create raster to visualize
ras <- rast(ncol = n_col, nrow = n_row,
            xmin = 0, xmax = n_col * 30, ymin = 0, ymax = n_row * 30)
ras_rc <- c(nrow(ras), ncol(ras))
cell <- sample(1:(ncol(ras)*nrow(ras)), size = 1)
rc <- cell_to_rowcol(cell, ras_rc)

# this function takes as input the [row, col] indexes of a burning pixel
adjacent_matrix <- function(row_col, n_rowcol) {
  
  #n_rowcol = ras_rc
  
  # neighbours row_col 
  neigh_rc <- row_col + moves # neighbours row-column pairs
 
  # Make NA the invalid positions
  for(i in 1:8) {
    if(! (neigh_rc[1, i] > 0 & neigh_rc[1, i] <= n_rowcol[1] &
          neigh_rc[2, i] > 0 & neigh_rc[2, i] <= n_rowcol[2])) {
       neigh_rc[, i] <- NA
    }
  }
  
  return(neigh_rc)
}

# create raster to visualize
ras <- rast(ncol = n_col, nrow = n_row,
            xmin = 0, xmax = n_col * 30, ymin = 0, ymax = n_row * 30)
ras_rc <- c(nrow(ras), ncol(ras))

cell <- sample(1:(ncol(ras)*nrow(ras)), size = 1)
row_col <- rowColFromCell(ras, cell) %>% as.numeric()

tmp <- adjacent_matrix(row_col, ras_rc)
tmp2 <- sapply(1:ncol(tmp), function(i) cellFromRowCol(ras, tmp[1, i], tmp[2, i]))

values(ras) <- 0
values(ras)[cell] <- 1
values(ras)[tmp2] <- 2
plot(ras, col = c("black", "red", "gray"), main = paste("cell =", cell))

# Vector, matrix and array sizes -------------------------------------------
# (it's all the same)

draws <- rnorm(1900 * 1900 * 8)
v1 <- draws
v2 <- matrix(draws, 1900, 1900 * 8)
v3 <- array(draws, dim = c(1900, 1900, 8))

