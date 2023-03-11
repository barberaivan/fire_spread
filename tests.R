library(testthat)

library(Rcpp)
library(terra)


sourceCpp("spread_functions.cpp")
source("spread_functions.R")
sourceCpp("similarity_functions.cpp")
source("similarity_functions.R")

# create data for testing
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


ig_location <- matrix(rep(round(size / 2), 2), 2, 1)


# simulate landscapes where one layer varies and check that results are OK.

land_cube <- function(x) {
  v <- values(x)
  a <- array(NA, dim = c(nrow(x), ncol(x), nlyr(x)),
             dimnames = list(row = NULL, col = NULL, layer = names(x)))
  for(l in 1:nlyr(x)) a[, , l] <- matrix(v[, l], nrow(x), ncol(x), byrow = TRUE)
  return(a)
}

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
c_sub["subalpine"] <- -1.2 # subalpine is not flammable
c_sub["intercept"] <- 0 # shrubland is very flammable
c_sub["wind"] <- 0 # remove wind effect

land_arr <- land_cube(lands_sub)

distances <- rep(res(landscape)[1], 8) # sides
distances[c(1, 3, 6, 8)] <- res(landscape)[1] * sqrt(2)


test_that("Fire spread functions", {

  set.seed(30)
  fire_r <- simulate_fire_r(
    landscape = lands_sub,
    burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
    ignition_cells = ig_location,
    coef = c_sub,
    wind_column = wind_column,
    elev_column = elev_column,
    distances = distances,
    upper_limit = 1.0
  )

  set.seed(30)
  fire_cpp <- simulate_fire_cpp(
    landscape = land_arr,
    burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
    ignition_cells = ig_location - 1,
    coef = c_sub,
    wind_layer = wind_column - 1,
    elev_layer = elev_column - 1,
    distances = distances,
    upper_limit = 1.0
  )

  set.seed(30)
  fire_compare_cpp <- simulate_fire_compare(
    landscape = land_arr,
    burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
    ignition_cells = ig_location - 1,
    coef = c_sub,
    wind_layer = wind_column - 1,
    elev_layer = elev_column - 1,
    distances = distances,
    upper_limit = 1.0
  )

  expect_equal(fire_r, fire_cpp)
  expect_equal(fire_r, fire_compare_cpp$burned_layer)
})

test_that("Deterministic fire spread functions", {
  fire_r <- simulate_fire_deterministic_r(
    landscape = lands_sub,
    burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
    ignition_cells = ig_location,
    coef = c_sub,
    wind_column = wind_column,
    elev_column = elev_column,
    distances = distances,
    upper_limit = 1.0
  )

  fire_cpp <- simulate_fire_deterministic_cpp(
    landscape = land_arr,
    burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
    ignition_cells = ig_location - 1,
    coef = c_sub,
    wind_layer = wind_column - 1,
    elev_layer = elev_column - 1,
    distances = distances,
    upper_limit = 1.0
  )

  expect_equal(fire_r, fire_cpp)
})

test_that("Similarity functions", {
  set.seed(1)
  fire_1 <- simulate_fire_compare(
    landscape = land_arr,
    burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
    ignition_cells = ig_location - 1,
    coef = c_sub,
    wind_layer = wind_column - 1,
    elev_layer = elev_column - 1,
    distances = distances,
    upper_limit = 1.0
  )

  set.seed(1)
  fire_1_ <- simulate_fire_compare(
    landscape = land_arr,
    burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
    ignition_cells = ig_location - 1,
    coef = c_sub,
    wind_layer = wind_column - 1,
    elev_layer = elev_column - 1,
    distances = distances,
    upper_limit = 1.0
  )

  set.seed(2)
  fire_2 <- simulate_fire_compare(
    landscape = land_arr,
    burnable = matrix(1, nrow(land_arr), ncol(land_arr)),
    ignition_cells = ig_location - 1,
    coef = c_sub,
    wind_layer = wind_column - 1,
    elev_layer = elev_column - 1,
    distances = distances,
    upper_limit = 1.0
  )

  expect_equal(fire_1, fire_1_)

  similarity_cpp_1_2 <- compare_fires_try(fire_1, fire_2)
  similarity_r_1_2 <- compare_fires_r(fire_1, fire_2)

  expect_equal(similarity_cpp_1_2, similarity_r_1_2)

  similarity_1_1 <- compare_fires_try(fire_1, fire_1)

  expect_equal(unname(similarity_1_1), rep(1, length(similarity_1_1)))
})
