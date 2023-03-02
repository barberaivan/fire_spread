library(Rcpp)
library(terra)
library(tidyverse)

# empezaré traduciendo a cpp las funciones de adyacencia.
# igual creo datos para ir probando las funciones.

# Create data -------------------------------------------------------------

# model coefficients (includes intercept)
coefs <- c(0.5, -1.5, 0.5, 0.5, 0.5, 0.5, 0.5)
names(coefs) <- c("Intercept", "veg", "elev", "aspect", "wind", "pp", "temp")
elev_column <- which(names(coefs) == "elev") - 1
wind_column <- which(names(coefs) == "wind") - 1 # -1 because the design matrix will have no intercept

# test wind effect
# coefs <- c(0.5, 0, 0, 0, 3, 0, 0)
# test slope effect
# coefs <- c(0.5, 0, 10, 0, 0, 0, 0)

# predictors raster

size <- 25

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
values(r)[ig_location:(ig_location+5)] <- 1
# plot(r_predictors)

# initialize burning_cells cell id
burning_cells <- which(values(r) == 1) # cell number id


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

sourceCpp("Spread function/rcpp_spread_functions.cpp")
n_rowcol <- c(5, 5)
cell_test <- sample(1:prod(n_rowcol), size = 3, replace = F); cell_test
matrix(1:prod(n_rowcol), n_rowcol[1], n_rowcol[2], byrow = TRUE)
ttt <- cell_to_rowcol_cpp(cell_test, n_rowcol = n_rowcol)
ttt
rowcol_to_cell_cpp(ttt, n_rowcol = n_rowcol)




# rowcol has to be a matrix (rowcol[1, ] = rows, rowcol[2, ] = columns)
rowcol_to_cell <- function(rowcol, n_rowcol) {
  return((rowcol[1, ] - 1) * n_rowcol[2] + rowcol[2, ])
}


# Adjacent (not still in c++) ---------------------------------------------

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



# Spread probability function -----------------------------------------------

spread_prob <- function(data_burning, data_neighbours, upper_limit = 1, cols_use) {
  
  # recompute wind and elev columns (will be slope and wind effects)
  data_neighbours[, "elev"] <- sin(atan((data_neighbours[, "elev"] - data_burning[1, "elev"]) / distances[cols_use]))
  # CAREFUL: data_burning was a df, so we need to select its unic 
  # and first row: [1, "elev"]
  
  data_neighbours[, "wind"] <- cos(angles[cols_use] - data_burning[1, "wind"])
  # the same here, watch the [1, "wind"]
  
  # careful with upper limit
  probs <- plogis(coefs[1] + as.matrix(data_neighbours) %*% coefs[2:length(coefs)]) * upper_limit
  return(probs %>% as.numeric)
}


# Test spread_prob_cpp ----------------------------------------------------

# get cell id from neighbours
rc = c(nrow(r), ncol(r))
neighbours_matrix <- adjacent_vector(cells = burning_cells, n_rowcol = rc) 
# for burning_cell[1]:
db <- r_predictors[burning_cells[1]] # data burning
dn <- r_predictors[neighbours_matrix[1, ]] # data neighbours
wind_column

#colnames(neighbours_matrix) <- 1:8


burning_cells

m <- matrix(NA, nrow = nrow(dn), ncol = ncol(dn))
for(i in 1:ncol(m)) m[, i] <- dn[, i] %>% as.numeric
# matrix(as.numeric(as.matrix(dn)), nrow = nrow(dn))

sourceCpp("Spread function/rcpp_spread_functions2.cpp")
spread_prob_cpp(data_burning = as.numeric(db), 
                data_neighbours = m,
                coef = unname(coefs),
                columns = (1:8) - 1,
                wind_column = as.integer(3),
                elev_column = as.integer(1), 
                distance = distances, 
                angles = as.numeric(angles),
                upper_limit = 1
                )
# careful: give the arguments in its simplest form, so c++ does not complain.

sourceCpp("Spread function/rcpp_spread_functions.cpp")
s <- rnbinom(1, 10, 0.50)
set.seed(s)
fff2 <- spread_prob_cpp(data_burning = as.numeric(db), 
                data_neighbours = as.matrix(dn),
                coef = unname(coefs),
                columns = (1:8) - 1,
                wind_column = as.integer(3),
                elev_column = as.integer(1), 
                distance = distances, 
                angles = as.numeric(angles),
                upper_limit = 1
)

pp <- spread_prob(data_burning = db, data_neighbours = dn, upper_limit = 1, cols_use = 1:8)
# now they match in probs.

set.seed(s)
fff <- rbinom(8, size = 1, prob = pp)
all.equal(fff, fff2)

# perfect, the stocastic part is the same (testing spread prob functions.)
# The approaches are comparable.



# Adjacent function -------------------------------------------------------

# This function will return the cell ids adjacent to a given set of burning cells,
# as a n_cell X 8 matrix (queen's neighbours). 
# The input data will be the cell ids of burning cells and nrow and ncol of the 
# landscape raster stack that is represented as a matrix

# R function

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

# The Rcpp version (made for R indexes, starting at 1):
sourceCpp("Spread function/rcpp_spread_functions_nesting.cpp")
adjacent_cpp(cells = 5, n_rowcol = n_rowcol, moves = moves)
adjacent_vector(cells = 5, n_rowcol = n_rowcol)

n_rowcol <- sample(5:100, 2, replace = TRUE)
ccc <- sample(1:prod(n_rowcol), 
              size = sample(1:round(prod(n_rowcol) * 0.1), size = 1), 
              replace = FALSE)
all.equal(
ad1 <- adjacent_cpp(cells = ccc, n_rowcol = n_rowcol, moves = moves),
ad2 <- adjacent_vector(cells = ccc, n_rowcol = n_rowcol)
)

# plot rasters
r_ex1 <- r_ex2 <- rast(ncol = n_rowcol[2], nrow = n_rowcol[1],
              xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000)
values(r_ex1) <- values(r_ex2) <- 0
values(r_ex1)[ccc] <- values(r_ex2)[ccc] <- 1

n1 <- unique(as.numeric(ad1)) %>% na.omit
n1 <- n1[!n1 %in% ccc]

n2 <- unique(as.numeric(ad2)) %>% na.omit
n2 <- n2[!n2 %in% ccc]

values(r_ex1)[n1] <- 2
values(r_ex2)[n2] <- 2

par(mfrow = c(1, 2))
plot(r_ex1, col = c("green", "red", "black"))
plot(r_ex2, col = c("green", "red", "black"))
par(mfrow = c(1, 1))

# Next step is simulate spread. As everything will be done in c++, indexes will 
# have to start at 0. A very simple way is to substract 1 from the neigh_cells
# result, but we can't be sure. 

# test cell to rowcol
sourceCpp("Spread function/rcpp_spread_functions_nesting_i0.cpp")

n_rowcol <- sample(4:100, size = 2, replace = T)
cell_test <- sample(1:prod(n_rowcol), size = 3, replace = F); cell_test
#matrix(1:prod(n_rowcol), n_rowcol[1], n_rowcol[2], byrow = TRUE)
ttt <- cell_to_rowcol_cpp(cell_test, n_rowcol = n_rowcol)
ttt
ttt0 <- cell_to_rowcol_cpp0(cell_test - 1, n_rowcol = n_rowcol)
ttt0
ttt1 <- cell_to_rowcol(cell_test, n_rowcol = n_rowcol) %>% as.matrix
all.equal(ttt, ttt0+1)
all.equal(ttt, ttt1)
# done

ttt
cell_test
rowcol_to_cell(ttt, n_rowcol)
rowcol_to_cell_cpp(ttt, n_rowcol)
rowcol_to_cell_cpp0(ttt0, n_rowcol) + 1

all.equal(rowcol_to_cell(ttt, n_rowcol),
          rowcol_to_cell_cpp0(ttt0, n_rowcol) + 1)
all.equal(rowcol_to_cell(ttt, n_rowcol),
          rowcol_to_cell_cpp(ttt, n_rowcol))


# Test adjacent0

# The Rcpp version (made for R indexes, starting at 1):
sourceCpp("Spread function/rcpp_spread_functions_nesting_i0.cpp")

n_rowcol <- sample(5:100, 2, replace = TRUE)
ccc <- sample(1:prod(n_rowcol), 
              size = sample(1:round(prod(n_rowcol) * 0.1), size = 1), 
              replace = FALSE)
all.equal(
  ad1 <- adjacent_cpp(cells = ccc, n_rowcol = n_rowcol, moves = moves),
  ad2 <- adjacent_vector(cells = ccc, n_rowcol = n_rowcol)
)
all.equal(
  ad1 <- adjacent_cpp(cells = ccc, n_rowcol = n_rowcol, moves = moves),
  ad3 <- adjacent_cpp0(cells = ccc - 1, n_rowcol = n_rowcol, moves = moves) + 1
)

# plot rasters
r_ex1 <- r_ex2 <- r_ex3 <- rast(ncol = n_rowcol[2], nrow = n_rowcol[1],
                       xmin = -1000, xmax = 1000, ymin = -1000, ymax = 1000)
values(r_ex1) <- values(r_ex2) <- values(r_ex3) <- 0
values(r_ex1)[ccc] <- values(r_ex2)[ccc] <- values(r_ex3)[ccc] <- 1

n1 <- unique(as.numeric(ad1)) %>% na.omit
n1 <- n1[!n1 %in% ccc]

n2 <- unique(as.numeric(ad2)) %>% na.omit
n2 <- n2[!n2 %in% ccc]

n3 <- unique(as.numeric(ad3)) %>% na.omit
n3 <- n3[!n3 %in% ccc]

values(r_ex1)[n1] <- 2
values(r_ex2)[n2] <- 2
values(r_ex3)[n3] <- 2

par(mfrow = c(1, 3))
plot(r_ex1, col = c("green", "red", "black"))
plot(r_ex2, col = c("green", "red", "black"))
plot(r_ex3, col = c("green", "red", "black"))
par(mfrow = c(1, 1))

# good. 
