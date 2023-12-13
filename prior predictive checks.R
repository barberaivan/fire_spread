# Code to define priors, at least placeholder ones to perform the
# MVN-ABC pseudolikelihood estimation.

# The priors should try to be similar across parameters.

# vfi and tfi are standardized, although not at the fire scale, but at the
# regional scale.
# slope ranges from -1 to 1, and wind between 0 and 15.

# Packages and data -------------------------------------------------------

library(terra)
library(tidyverse)
library(logitnorm)
library(LaplacesDemon)
library(FireSpread)

source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for rast_from_mat and a few constants

# load Cholila fire data. The raster is used as a template for plotting.
land <- readRDS(file.path("data", "focal fires data",
                          "landscapes_ig-known", "2015_50.rds"))
land_raster <- rast(file.path("data", "focal fires data", "raw data from GEE",
                              "fire_data_raw_2015_50.tif"))



# Get windspeed scale in all landscapes -----------------------------------

files <- list.files(file.path("data", "focal fires data", "landscapes_ig-known"))
wsumm <- function(x) {
 c(
   "min" = min(x, na.rm = TRUE),
   "max" = max(x, na.rm = TRUE),
   "sd" = sd(x, na.rm = TRUE),
   "mean" = mean(x, na.rm = TRUE),
   "median" = median(x, na.rm = TRUE)
 )
}

wind_values <- do.call("rbind", lapply(1:length(files), function(i) {
  # i = 1
  print(i)
  ll <- readRDS(file.path("data", "focal fires data", "landscapes_ig-known", files[i]))
  return(wsumm(ll$landscape[, , "wspeed"]))
}))
rm(ll); gc()

summary(wind_values) # summary across fires
# most vary between 0 and 10


# Constants ---------------------------------------------------------------

upper_limit <- 0.5

# Logistic curves for standardized variables ------------------------------

sd_z <- 5
r_z <- 0.15
sd_log <- 1.5
mean_log <- -0.5

par(mfrow = c(1, 2))
#normal
curve(upper_limit * plogis(abs(rnorm(1, 0, sd_z)) * x),
      main = paste("Normal(sd = ", sd_z, ")", sep = ""),
      from = -5, to = 5, col = rgb(0, 0, 0, 0.1),
      ylab = "burn prob",
      xlab = "standardized variable",
      ylim = c(0, upper_limit))
for(i in 1:1000) {
  curve(upper_limit * plogis(abs(rnorm(1, 0, sd_z)) * x), add = TRUE,
        col = rgb(0, 0, 0, 0.1))
}

# exponential
curve(upper_limit * plogis(rexp(1, r_z) * x),
      main = paste("exponential (r = ", r_z, ")", sep = ""),
      from = -5, to = 5, col = rgb(0, 0, 0, 0.1),
      ylab = "burn prob",
      xlab = "standardized variable",
      ylim = c(0, upper_limit))
for(i in 1:1000) {
  curve(upper_limit * plogis(rexp(1, r_z) * x), add = TRUE,
        col = rgb(0, 0, 0, 0.1))
}

# # lognormal # HARD
# curve(upper_limit * plogis(rlnorm(1, mean_log, sd_log) * x),
#       main = paste("lognormal (mean = ", mean_log, " sd = ", sd_log, ")", sep = ""),
#       from = -5, to = 5, col = rgb(0, 0, 0, 0.1),
#       ylab = "burn prob",
#       xlab = "standardized variable",
#       ylim = c(0, upper_limit))
# for(i in 1:1000) {
#   curve(upper_limit * plogis(rlnorm(1, mean_log, sd_log) * x), add = TRUE,
#         col = rgb(0, 0, 0, 0.1))
# }

par(mfrow = c(1, 1))

# densities
curve(dnorm(x, 0, sd_z) * 2, from = 0, to = 20, ylim = c(0, 0.3))
curve(dexp(x, r_z), col = 2, add = TRUE)
curve(dlnorm(x, 3, 2), col = 4, add = TRUE)   # hard to mimic with log-normal

# fit log-normal to exponential:
rsamps <- rexp(1e6, 0.1)
log_samps <- log(rsamps)
mm <- lm(log_samps ~ 1)
curve(dexp(x, 0.25), col = 2, from = 0, to = 20)
curve(dlnorm(x, log(coef(mm)), sigma(mm)), col = 4, add = TRUE) # feo feo el fit

# Logistic curves for [-1, 1] variables ---------------------------------

sd_01 <- 10
r_01 <- 0.04

par(mfrow = c(1, 2))

#normal
curve(upper_limit * plogis(abs(rnorm(1, 0, sd_01)) * x),
      main = paste("Normal(sd = ", sd_01, ")", sep = ""),
      from = -1, to = 1, col = rgb(0, 0, 0, 0.1),
      ylab = "burn prob",
      xlab = "[-1, 1] variable",
      ylim = c(0, upper_limit))
for(i in 1:1000) {
  curve(upper_limit * plogis(abs(rnorm(1, 0, sd_01)) * x), add = TRUE,
        col = rgb(0, 0, 0, 0.1))
}

# exponential
curve(upper_limit * plogis(rexp(1, r_01) * x),
      main = paste("exponential (r = ", r_01, ")", sep = ""),
      from = -1, to = 1, col = rgb(0, 0, 0, 0.1),
      ylab = "burn prob",
      xlab = "[-1, 1] variable",
      ylim = c(0, upper_limit))
for(i in 1:1000) {
  curve(upper_limit * plogis(rexp(1, r_01) * x), add = TRUE,
        col = rgb(0, 0, 0, 0.1))
}

par(mfrow = c(1, 1))


# densities
curve(dnorm(x, 0, sd_01) * 2, from = 0, to = 20,
      ylim = c(0, 0.18))
curve(dexp(x, r_01), col = 2, add = TRUE)

# Logistic curves for wind [0, 15] ---------------------------------

# 0.5 is frequently the wind sd
sd_wind <- sd_z * 0.5
r_wind <- r_z / 0.5
wupr <- 10

par(mfrow = c(1, 2))

#normal
curve(upper_limit * plogis(abs(rnorm(1, 0, sd_wind)) * x),
      main = paste("Normal(sd = ", sd_wind, ")", sep = ""),
      from = -wupr, to = wupr, col = rgb(0, 0, 0, 0.1),
      ylab = "burn prob",
      xlab = "wind",
      ylim = c(0, upper_limit))
for(i in 1:1000) {
  curve(upper_limit * plogis(abs(rnorm(1, 0, sd_wind)) * x), add = TRUE,
        col = rgb(0, 0, 0, 0.1))
}

# exponential
curve(upper_limit * plogis(rexp(1, r_wind) * x),
      main = paste("exponential (r = ", r_wind, ")", sep = ""),
      from = -wupr, to = wupr, col = rgb(0, 0, 0, 0.1),
      ylab = "burn prob",
      xlab = "wind",
      ylim = c(0, upper_limit))
for(i in 1:1000) {
  curve(upper_limit * plogis(rexp(1, r_wind) * x), add = TRUE,
        col = rgb(0, 0, 0, 0.1))
}

par(mfrow = c(1, 1))



# Intercepts distribution -------------------------------------------------

sd_int <- 5
mu_int <- 0#logit(0.2 / upper_limit)
curve(dlogitnorm(x, mu_int, sd_int), n = 1000)

# Function to simulate from the prior -------------------------------------

prior_dist <- function(mu_int = 0, sd_int = 10,
                       r_slope = 0.04,
                       r_fi = 0.15,
                       r_wind = 0.3,
                       type = "rng", # or "quantile"
                       n = 1,
                       p = NULL) {

  if(type == "rng") {
    b <- cbind(
      "intercept" = rnorm(n, mu_int, sd_int),
      "b_vfi" = rexp(n, r_fi),
      "b_tfi" = rexp(n, r_fi),
      "b_slope" = rexp(n, r_slope),
      "b_wind" = rexp(n, r_wind)
    )

    if(n == 1) {
      parnames <- colnames(b)
      b <- as.numeric(b)
      names(b) <- parnames
    }
    return(b)
  }

  if(type == "quantile") {

    if(is.null(p)) stop("Provide probability to compute quantiles.")
    if(length(p) %% 5 != 0) stop("p must be a multiple of the number of parameters (5).")

    if(length(dim(p)) < 2) p <- matrix(p, ncol = 5)
    if(length(dim(p)) == 2) {
      if(ncol(p) != 5) p <- matrix(as.numeric(p), ncol = 5)
    }

    q <- cbind(
      "intercept" = qnorm(p[, 1], mu_int, sd_int),
      "b_vfi" = qexp(p[, 2], r_fi),
      "b_tfi" = qexp(p[, 3], r_fi),
      "b_slope" = qexp(p[, 4], r_slope),
      "b_wind" = qexp(p[, 5], r_wind)
    )

    return(q)
  }
}

prior_dist(n = 1)
prior_dist(n = 1, type = "quantile")
prior_dist(n = 1, type = "quantile", p = runif(4))
prior_dist(n = 1, type = "quantile", p = runif(10))
prior_dist(n = 1, type = "quantile", p = runif(20))
prior_dist(n = 1, type = "quantile", p = matrix(runif(20), 5, 4))

# Graphical prior check ---------------------------------------------------

# matrix to fill
hollow <- land$burnable - 1

(pp <- prior_dist())

set.seed(1)
fire_prior <- simulate_fire(
  landscape = land$landscape,
  burnable = land$burnable,
  ignition_cells = land$ig_rowcol,
  coef = pp,
  upper_limit = 0.5
)

hollow[fire_prior == 1] <- 1
r <- rast_from_mat(hollow, land_raster[[1]])
plot(r, col = c("black", "forestgreen", "orange"))
