library(testthat)

library(Rcpp)
library(terra)

sourceCpp("spread_functions.cpp")
source("spread_functions.R")
sourceCpp("similarity_functions.cpp")
source("similarity_functions.R")

# create data for testing
n_veg_types <- 6
n_terrain <- 4

set.seed(2345)
coefs <- c(
  rnorm(n_veg_types, 0),
  rnorm(n_terrain)
)

# terrain names
terrain_names <- c("northing", "elev", "windir")
terrain_params <- c(terrain_names, "slope")

# landscape raster
size <- 30
n_rows <- size
n_cols <- size
res <- 30
n_cells <- n_rows * n_cols

landscape <- rast(
  ncol = n_cols, nrow = n_rows, res = res, crs = "EPSG:5343",
  nlyrs = 1 + n_terrain - 1, # veg in just one layer, and slope absent
  xmin = 0, xmax = res * n_cols, ymin = 0, ymax = res * n_rows,
  names = c("vegetation", terrain_names)
)

# fill data
veg_vals <-  sample(c(0:(n_veg_types - 1), 99), size = n_cells, replace = TRUE)
vegetation <- matrix(veg_vals, n_rows, byrow = T)
values(landscape)[, 1] <- veg_vals
landscape$northing <- cos(runif(ncell(landscape), 0, 2 * pi) - 315 * pi / 180)
landscape$elev <- runif(ncell(landscape), 0, 2200) / elevation_sd # scaled
landscape$windir <- runif(ncell(landscape), 0, 2 * pi) # radians

ig_location <- matrix(rep(round(size / 2), 2), 2, 1)

test_that("Fire spread functions", {

  set.seed(30)
  fire_r <- simulate_fire_r(
    terrain = landscape[[terrain_names]], # use the SpatRaster
    vegetation = vegetation,
    ignition_cells = ig_location,
    coef = coefs,
    n_veg_types = n_veg_types,
    upper_limit = 1.0
  )

  set.seed(30)
  fire_cpp <- simulate_fire_cpp(
    terrain = land_cube(landscape), # use the array
    vegetation = vegetation,
    ignition_cells = ig_location - 1,
    coef = coefs,
    n_veg_types = n_veg_types,
    upper_limit = 1.0
  )

  set.seed(30)
  fire_compare_cpp <- simulate_fire_compare(
    terrain = land_cube(landscape), # use the array
    vegetation = vegetation,
    ignition_cells = ig_location - 1,
    coef = coefs,
    n_veg_types = n_veg_types,
    upper_limit = 1.0
  )

  expect_equal(fire_r, fire_cpp)
  expect_equal(fire_r, fire_compare_cpp$burned_layer)
})

test_that("Deterministic fire spread functions", {
  fire_r <- simulate_fire_deterministic_r(
    terrain = landscape[[terrain_names]], # use the SpatRaster
    vegetation = vegetation,
    ignition_cells = ig_location,
    coef = coefs,
    n_veg_types = n_veg_types,
    upper_limit = 1.0
  )

  fire_cpp <- simulate_fire_deterministic_cpp(
    terrain = land_cube(landscape), # use the array
    vegetation = vegetation,
    ignition_cells = ig_location - 1,
    coef = coefs,
    n_veg_types = n_veg_types,
    upper_limit = 1.0
  )

  expect_equal(fire_r, fire_cpp)
})

test_that("Similarity functions", {
  set.seed(1)
  fire_1 <- simulate_fire_compare(
    terrain = land_cube(landscape), # use the array
    vegetation = vegetation,
    ignition_cells = ig_location - 1,
    coef = coefs,
    n_veg_types = n_veg_types,
    upper_limit = 1.0
  )

  set.seed(1)
  fire_1_ <- simulate_fire_compare(
    terrain = land_cube(landscape), # use the array
    vegetation = vegetation,
    ignition_cells = ig_location - 1,
    coef = coefs,
    n_veg_types = n_veg_types,
    upper_limit = 1.0
  )

  set.seed(2)
  fire_2 <- simulate_fire_compare(
    terrain = land_cube(landscape), # use the array
    vegetation = vegetation,
    ignition_cells = ig_location - 1,
    coef = coefs,
    n_veg_types = n_veg_types,
    upper_limit = 1.0
  )

  expect_equal(fire_1, fire_1_)

  similarity_cpp_1_2 <- compare_fires_try(fire_1, fire_2)
  similarity_r_1_2 <- compare_fires_r(fire_1, fire_2)

  expect_equal(similarity_cpp_1_2, similarity_r_1_2, tolerance = 1e-6)

  similarity_1_1 <- compare_fires_try(fire_1, fire_1)

  expect_equal(unname(similarity_1_1), rep(1, length(similarity_1_1)))
})