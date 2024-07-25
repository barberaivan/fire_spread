# This code inherits from
# <single fire posterior - surrogate vs simulator - binary gam - many intercepts.R>.

# Here I run the simulation of small fires, to see how GAMs work.
# Evaluate 40000 particles, fit GAM, test gam on 2000 particles.

# The model with many intercepts has
# 5 veg types
# slope, wind
# steps

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(ggdensity) # geom_hdr

library(Rcpp)
library(terra)
library(abind)         # manage arrays

library(randtoolbox)   # sobol sequences
library(mgcv)          # fit GAM
library(MGMM)          # fit MVN distribution
library(Matrix)        # nearPD, searches the closest PD matrix.
library(stutils)       # tryNULL, for search of nearPD
library(logitnorm)
library(rvinecopulib)  # density estimation

library(sn)
library(trialr) # rlkjcorr
library(truncnorm)

library(adaptMCMC)
library(posterior)     # manage posterior samples
library(tidybayes)     # not sure if this was used
library(bayesplot)     # visualize posteriors

library(foreach)       # parallelization
library(doMC)          # doRNG had problems with cpp functions

library(FireSpread)    # spread and similarity functions

library(microbenchmark)

source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for rast_from_mat and a few constants

# source("estimation_functions.R") # prior_dist and other stuff

# Constants --------------------------------------------------------------

data_dir <- file.path("data", "focal fires data", "landscapes")
filenames <- list.files(data_dir)
n_fires <- length(filenames)

# dir to save output
target_dir <- file.path("files", "overlaps")

# load file with constants to standardize
ndvi_params <- readRDS(file.path("data", "NDVI_regional_data",
                                 "ndvi_optim_and_proportion.rds"))
slope_sd <- ndvi_params$slope_term_sd

# constants for fire spread simulation
upper_limit <- 1
n_veg <- 5
veg_names <- c("wet", "subalpine", "dry", "shrubland", "grassland")
n_terrain <- 2
terrain_names <- c("slope", "wind")
par_names <- c(veg_names, terrain_names, "steps")
n_coef <- length(par_names)

# number of fires to simulate by particle
n_rep <- 20

ext_alpha <- 50
ext_beta <- 30

params_lower <- c(rep(-ext_alpha, n_veg), rep(0, n_terrain), 5)
params_upper <- c(rep(ext_alpha, n_veg), ext_beta / slope_sd, ext_beta, NA)

names(params_lower) <- names(params_upper) <- par_names

gam_formula_full <- formula(
  y ~
    # marginal effects
    s(wet, bs = basis, k = k_side) +
    s(subalpine, bs = basis, k = k_side) +
    s(dry, bs = basis, k = k_side) +
    s(shrubland, bs = basis, k = k_side) +
    s(grassland, bs = basis, k = k_side) +

    s(slope, bs = basis, k = k_side) +
    s(wind, bs = basis, k = k_side) +

    s(steps, bs = basis, k = k_side) +

    # interactions (intercepts)
    ti(wet, slope, k = k_int, bs = basis) +
    ti(wet, wind, k = k_int, bs = basis) +

    ti(subalpine, slope, k = k_int, bs = basis) +
    ti(subalpine, wind, k = k_int, bs = basis) +

    ti(dry, slope, k = k_int, bs = basis) +
    ti(dry, wind, k = k_int, bs = basis) +

    ti(shrubland, slope, k = k_int, bs = basis) +
    ti(shrubland, wind, k = k_int, bs = basis) +

    ti(grassland, slope, k = k_int, bs = basis) +
    ti(grassland, wind, k = k_int, bs = basis) +

    # directional interactions
    ti(slope, wind, k = k_int, bs = basis)
)

a <- "s(wet, k = k_side, bs = \"cr\")"
as.formula(paste("y ~", a))

gam_terms <- list(
  # marginals
  wet       = "s(wet, k = k_side, bs = \"cr\")",
  subalpine = "s(subalpine, k = k_side, bs = \"cr\")",
  dry       = "s(dry, k = k_side, bs = \"cr\")",
  shrubland = "s(shrubland, k = k_side, bs = \"cr\")",
  grassland = "s(grassland, k = k_side, bs = \"cr\")",

  slope = "s(slope, k = k_side, bs = \"cr\")",
  wind  = "s(wind, k = k_side, bs = \"cr\")",
  steps = "s(steps, k = k_side, bs = \"cr\")",

  # interactions vegetation_slope
  wet_slope       = "ti(wet, slope, k = k_int, bs = \"cr\")",
  subalpine_slope = "ti(subalpine, slope, k = k_int, bs = \"cr\")",
  dry_slope       = "ti(dry, slope, k = k_int, bs = \"cr\")",
  shrubland_slope = "ti(shrubland, slope, k = k_int, bs = \"cr\")",
  grassland_slope = "ti(grassland, slope, k = k_int, bs = \"cr\")",

  # interactions vegetation_wind
  wet_wind       = "ti(wet, wind, k = k_int, bs = \"cr\")",
  subalpine_wind = "ti(subalpine, wind, k = k_int, bs = \"cr\")",
  dry_wind       = "ti(dry, wind, k = k_int, bs = \"cr\")",
  shrubland_wind = "ti(shrubland, wind, k = k_int, bs = \"cr\")",
  grassland_wind = "ti(grassland, wind, k = k_int, bs = \"cr\")",

  # other interactions
  slope_wind  = "ti(slope, wind, k = k_int, bs = \"cr\")",
  slope_steps = "ti(slope, steps, k = k_int, bs = \"cr\")",
  wind_steps  = "ti(wind, steps, k = k_int, bs = \"cr\")"
)

terms_names <- names(gam_terms)

# Get landscapes size -----------------------------------------------------

# n_fires <- length(filenames)
# size_data <- data.frame(file = filenames,
#                         fire_id = NA,
#                         size_land = NA,
#                         size_burnable = NA)
# for(i in 1:n_fires) {
#   ll <- readRDS(file.path(data_dir, filenames[i]))
#   size_data$fire_id[i] <- ll$fire_id_spread
#   size_data$size_land[i] <- prod(dim(ll$landscape)[1:2])
#   size_data$size_burnable[i] <- sum(ll$cells_by_veg)
# }
# rm(ll); gc()
# size_data <- size_data[order(size_data$size_burnable), ]
# size_data$size_burn_rel <- size_data$size_burnable / max(size_data$size_burnable)
# saveRDS(size_data, file.path("data", "focal fires data", "fire_size_data.rds"))
size_data <- readRDS(file.path("data", "focal fires data", "fire_size_data.rds"))

# rownames(size_data) <- 1:nrow(size_data)
# write.csv(size_data, file.path("data", "focal fires data", "fire_size_data.csv"))

# Functions ---------------------------------------------------------------

normalize <- function(x) x / sum(x)

# Simulate fires and compare then with the observed one using the overlap_spatial
# function. The landscape argument includes all data to simulate fire, and also
# burned and burned_ids layers.
# Returns a matrix with the mean and variance across replicates of:
#   overlap_spatial,
#   size_diff: size differences, as simulated - observed,
#   edge: number of pixels burned at the edge of the landscape.
# The last two are used to rule out bounded fires.
similarity_simulate_particle <- function(particle, fire_data = NULL,
                                         n_sim = n_rep) {

  ## testo
  # particle <- particles_sim(N = 1)
  ## end testo
  ov <- numeric(n_sim)
  cv <- particle[1:n_veg]
  ct <- particle[(n_veg+1):(n_veg+n_terrain)]

  # simulate and compute metrics
  for(i in 1:n_sim) { # simulated fires by particle
    fire_sim <- simulate_fire_compare(
      vegetation = fire_data$landscape[, , 1],
      terrain = fire_data$landscape[, , -1],
      coef_veg = cv,
      coef_terrain = ct,
      ignition_cells = fire_data$ig_rowcol,
      upper_limit = upper_limit,
      steps = particle[n_coef]
    )

    ov[i] <- overlap_spatial(
      fire_sim, fire_data[c("burned_layer", "burned_ids")]
    )
  }

  return(ov)
}

# The same but returning many metrics. Currently used only for steps used.
# It returns a data.frame with averages.
similarity_simulate_particle_metrics <- function(particle, fire_data = NULL,
                                                 n_sim = n_rep) {

  metrics <- matrix(NA, n_sim, 3)
  colnames(metrics) <- c("overlap", "size_diff", "steps_used")
  cv <- particle[1:n_veg]
  ct <- particle[(n_veg+1):(n_veg+n_terrain)]

  # simulate and compute metrics
  for(i in 1:n_sim) { # simulated fires by particle
    fire_sim <- simulate_fire_compare(
      vegetation = fire_data$landscape[, , 1],
      terrain = fire_data$landscape[, , -1],
      coef_veg = cv,
      coef_terrain = ct,
      ignition_cells = fire_data$ig_rowcol,
      upper_limit = upper_limit,
      steps = particle[n_coef]
    )

    metrics[i, "overlap"] <- overlap_spatial(
      fire_sim, fire_data[c("burned_layer", "burned_ids")]
    )

    metrics[i, "size_diff"] <- ncol(fire_sim$burned_ids) -
      ncol(fire_data$burned_ids)

    metrics[i, "steps_used"] <- fire_sim$steps_used
  }

  mo <- mean(metrics[, "overlap"])
  # extract values
  ll_summ <- c(
    overlap = mo,
    overlap_log = log(mo),
    overlap_logit = qlogis(mo),
    overlap_var = var(metrics[, "overlap"]),       # variance in original scale
    ll = qlogis(mo),
    var = var(metrics[, "overlap"]), # same as overlap_var, but kept for compatibility
    size_diff = mean(metrics[, "size_diff"]),
    steps_used = mean(metrics[, "steps_used"])
  )

  return(ll_summ)
}

# Emulate the loglik over a list of particles in parallel. Returns the metrics
# in rows.
# It requires the matrix of particles to simulate (coefficient values).
similarity_simulate_parallel <- function(particles_mat = NULL,
                                         fire_data = NULL,
                                         n_sim = n_rep) {

  # turn particle matrix into list for parallel evaluation
  particles_list <- lapply(1:nrow(particles_mat),
                           function(x) particles_mat[x, ])

  # simulate in parallel
  result <- foreach(pp = particles_list) %dopar% {
    similarity_simulate_particle(pp, fire_data, n_sim)
  }

  # rbind list result
  overlap_mat <- do.call("rbind", result) %>% as.matrix()

  # return data.frame
  res <- data.frame(
    wave = NA,
    overlap = rowMeans(overlap_mat)
  )
  res$par_values <- particles_mat
  res$ov_values <- overlap_mat


  return(res)
}

# get_bounds: computes the metrics for the largest and smallest fires possible
get_bounds <- function(fire_data, n_sim = n_rep) {

  coef_burn_all <- c(1e6, rep(0, n_coef - 1))
  coef_burn_none <- c(-1e6, rep(0, n_coef - 1))

  small_fire <- similarity_simulate_particle_metrics(coef_burn_none,
                                                     fire_data = fire_data,
                                                     n_sim = n_sim)

  large_fire <- similarity_simulate_particle_metrics(coef_burn_all,
                                                     fire_data = fire_data,
                                                     n_sim = n_sim)

  sim_bounds <- rbind(small_fire, large_fire)
  rownames(sim_bounds) <- c("smallest", "largest")

  return(sim_bounds)
}

wave_plot <- function(data, response = "overlap", alpha = 0.3, best = NULL,
                      x = "par_values", thres = NULL, bin = FALSE,
                      tit = NULL, rc = c(3, 3)) {

  if(!is.null(best)) {
    data <- data[order(data[, response], decreasing = TRUE), ]
  }

  if(bin) {
    data$bin <- jitter(as.numeric(data$overlap >= thres), factor = 0.2)
    response <- "bin"
  }

  yy <- range(data[, response])
  yy[2] <- yy[2] * 1.05

  # title for first plot
  if(is.null(tit)) {
    tit <- deparse(substitute(data))
  }

  par(mfrow = c(rc[1], rc[2]))
  for(i in 1:n_coef) {
    mm <- ifelse(i == 1, tit, NA)
    plot(data[, response] ~ data[, x][, i], ylab = response,
         xlab = par_names[i], ylim = yy, main = mm,
         pch = 19, col = rgb(0, 0, 0, alpha))

    if(!is.null(best)) {
      points(data[1:best, response] ~ data[, x][1:best, i],
             pch = 19, col = rgb(0, 0, 1, alpha))
    }

    if(!is.null(thres) & response == "overlap") {
      abline(h = thres, col = "red")
    }
  }
  par(mfrow = c(1, 1))
}

wave_plot_pairs <- function(data, alpha = 0.05, thres, support = NULL,
                            tit = NA, rc = c(3, 5)) {

  ## TEST
  # data = wlist$like_sim
  # thres = wlist$thres
  # alpha = 0.1
  # tit = "1999_2140469994_r"
  # rc = c(3, 5)
  ##

  if(is.null(support)) {
    ranges <- apply(data$par_values, 2, range)
  } else {
    ranges <- support
  }

  xvars <- list(
    rep("slope", n_veg),
    rep("wind", n_veg),
    c("wind", "steps", "steps")
  )

  yvars <- list(
    veg_names,
    veg_names,
    c("slope", "slope", "wind")
  )

  dd <- data$par_values[data$overlap >= thres, ]

  n_cols <- sapply(yvars, length)

  par(mfrow = c(rc[1], rc[2]))
  for(r in 1:3) {
    for(c in 1:n_cols[r]) {
      # r = 1; c = 1

      tit <- ifelse(r == 1 & c == 1, tit, NA)

      yname <- yvars[[r]][c]
      xname <- xvars[[r]][c]

      plot(dd[, yname] ~ dd[, xname], xlab = xname, ylab = yname,
           xlim = ranges[, xname], ylim = ranges[, yname],
           pch = 19, col = rgb(0, 0, 0, alpha),
           main = tit)
    }
  }
  par(mfrow = c(1, 1))
}

# The same, to be used only with particles (matrix)
wave_plot_pairs2 <- function(x, alpha = 0.05, support = NULL,
                            tit = NA, rc = c(3, 5)) {

  if(is.null(support)) {
    ranges <- apply(x, 2, range)
  } else {
    ranges <- support
  }

  xvars <- list(
    rep("slope", n_veg),
    rep("wind", n_veg),
    c("wind", "steps", "steps")
  )

  yvars <- list(
    veg_names,
    veg_names,
    c("slope", "slope", "wind")
  )

  dd <- x

  n_cols <- sapply(yvars, length)

  par(mfrow = c(rc[1], rc[2]))
  for(r in 1:3) {
    for(c in 1:n_cols[r]) {
      # r = 1; c = 1

      tit <- ifelse(r == 1 & c == 1, tit, NA)

      yname <- yvars[[r]][c]
      xname <- xvars[[r]][c]

      plot(dd[, yname] ~ dd[, xname], xlab = xname, ylab = yname,
           xlim = ranges[, xname], ylim = ranges[, yname],
           pch = 19, col = rgb(0, 0, 0, alpha),
           main = tit)
    }
  }
  par(mfrow = c(1, 1))
}



# https://rpubs.com/binhho660/922614
#' rmvn by hand, Covariance parameterization, or maybe from a sobol sequence
#' @param n number of samples to draw
#' @param mu p-vector mean or p-columns matrix of means.
#' @param Sigma covariance matrix (p x p)
#' @param sobol logical, indicate whether to use or not a sobol sequence.
#' @return matrix of dimension p x n of samples
rmvn_sobol <- function(n, mu, Sigma, sobol = F){
  if(is.null(dim(mu))) {
    mu <- matrix(mu, nrow = 1)
  }

  p <- ncol(mu)

  if(sobol) {
    if(nrow(mu) > 1) stop("Sobol sequence not recommended with more than one center.")
    P <- sobol(n, p, init = F)
    Z <- apply(P, 2, qnorm) %>% t()
  } else {
    Z <- matrix(rnorm(p*n), p, n)
  }

  L <- t(chol(Sigma)) # By default R's chol function returns upper cholesky factor
  X <- L %*% Z

  if(nrow(mu) == 1) { # global centre
    mu <- as.numeric(mu)
    X <- sweep(X, 1, mu, FUN = `+`)
  } else {            # focal centre
    X <- X + t(mu)
  }
  return(X)
}


# positive-definiteness check made by makeSECdistr {sn}. It differs from
# is.positive.definite(), from {matrixcalc}
is_positive_definite <- function(m) {
  min(eigen(m, symmetric = TRUE, only.values = TRUE)$values) > 0
}

# make_positive_definite makes a vcov positive definite. It's supposed that if
# the correlation matrix is pd, then the vcov contructed with any set of
# scales would be. But it isn't so in practice.
# https://stats.stackexchange.com/questions/23317/create-positive-definite-3x3-covariance-matrix-given-specified-correlation-value

# If the positive-definiteness fails after max_iter trials, the result is NA.
# This is useful for optim(), because otherwise the iteration number for finding
# a positive-definite matrix took too long. It's safer to just ignore problematic
# parameter vector
make_positive_definite <- function(m, corr = FALSE, max_iter = 1e4) {
  if(is_positive_definite(m)) {
    return(m)
  } else {
    temp <- tryNULL(nearPD(m, corr = corr, keepDiag = TRUE, maxit = max_iter))
    if(is.null(temp)) return(NA) else {
      mat <- as.matrix(temp$mat)
      if(!is_positive_definite(mat)) {
        return(NA)
      } else return(mat)
    }
  }
}

# mvn_fitted <- gam(
#   list(intercept ~ 1, vfi ~ 1, tfi ~ 1, slope ~ 1, wind ~ 1, steps ~ 1),
#   family = mvn(n_coef), data = as.data.frame(xun_seeds)
# )
# mus <- coef(mvn_fitted)[1:n_coef]
# V <- solve(crossprod(mvn_fitted$family$data$R))

# # Fit MVN by hand
# V <- cov(xun_seeds)
# mus <- apply(xun_seeds, 2, mean)
#
# # Ensure positive definiteness of enlarged fitted covariance matrix
# Vlarge <- V * var_factor
# Vlarge <- make_positive_definite(as.matrix(forceSymmetric(Vlarge)))
#
# # if we did not find a positive definite matrix, use a diagonal vcov.
# if(is.null(dim(Vlarge))) { # returns "NA" when could not find it
#   Vlarge <- diag(apply(xun_seeds, 2, var) * var_factor)
#   # again force numerical PDness
#   diag(Vlarge) <- diag(Vlarge) + 1e-12
# }
#
# # sample new particles using enlarged covariance
# xun_new <- apply(xun_seeds, 1, FUN = function(x) {
#   mgcv::rmvn(n = 1, mu = mus, V = Vlarge)
# }) %>% t




# invlogit_scaled: translates unconstrained variables (logit) to a flat compatc
#   support, by applying invlogit and then scaling.
# the support matrix has dimensions [c(lower, upper), var_names]
invlogit_scaled <- function(xun, support) {
  x <- xun
  for(j in 1:ncol(xun)) {
    x[, j] <- plogis(xun[, j]) * (support[2, j] - support[1, j]) + support[1, j]
  }
  return(x)
}

invlogit_scaled_vec <- function(xun, support) {
  dmax <- support[2, ] - support[1, ]
  x <- plogis(xun) * dmax + support[1, ]
  return(x)
}

# logit_scaled: translates variables with compact support to the logit space,
#   by previously scaling them between 0 and 1.
# the support matrix has dimensions [c(lower, upper), var_names]
logit_scaled <- function(x, support) {
  xun <- x
  for(j in 1:n_coef) {
    xun[, j] <- qlogis((x[, j] - support[1, j]) / (support[2, j] - support[1, j]))
  }
  return(xun)
}

logit_scaled_vec <- function(x, support) {
  dobs <- x - support[1, ]
  dmax <- support[2, ] - support[1, ]
  x01 <- dobs / dmax
  return(qlogis(x01))
}

# function to scale [0, 1] particles to a given support
scale_params <- function(x, support) {
  xs <- x
  for(i in 1:n_coef) {
    xs[, i] <- x[, i] * (support[2, i] - support[1, i]) + support[1, i]
  }
  return(xs)
}

# Function to transform from original (constrained) support to the
# unconstrained one, using scaled logit and log (for steps)
unconstrain <- function(x, support) {
  xun <- x

  names_logit <- colnames(x)[colnames(x) != "steps"]
  for(j in names_logit) {
    xun[, j] <- qlogis((x[, j] - support[1, j]) / (support[2, j] - support[1, j]))
  }

  xun[, "steps"] <- log(x[, "steps"])

  return(xun)
}

# The inverse of unconstrain
constrain <- function(xun, support) {
  xc <- xun

  names_logit <- colnames(xun)[colnames(xun) != "steps"]
  for(j in names_logit) {
    xc[, j] <- plogis(xun[, j]) * (support[2, j] - support[1, j]) + support[1, j]
  }

  xc[, "steps"] <- exp(xun[, "steps"])

  return(xc)
}



# explore_likelihood searches particles in a space with compact support.
# The algorithm starts with a dataset of simulated likelihood over a sobol
# sequence of particles, and then it resamples the best ones according to the
# likelihood. Good particles are reproduced with a MVN kernel.
# data: dataset with previous simulated particles.
# n: number of new particles to simulate.
# var_factor: multiplier of the covariance matrix computed from the best
#   particle. Ideally, between 1 and 2.
# n_best: number of best particles to reproduce. if "all", all data available
#   are resampled.
# p_best: the proportion of best particles to resample. Overrides n_best.
# accept_thres: overlap threshold below which particles will not be resampled.
#   overrides both n_best and p_best.
# prob_fn: function to transform overlap into probability. If not NULL,
#   the probability is used to resample particles.
# support: matrix with lower limits (first row) and upper limits (second row)
#   of all variables.
# spread_data: data to simulate fires.
# centre: "local" or "global". New particles are simulated from a MVN centred
#   at each resampled particle (focal), or at the global mean (global). If global,
#   sobol may be set to TRUE, if focal, sobol is set to FALSE.
explore_likelihood <- function(data, n = 600, var_factor = 1, n_best = "all",
                               p_best = 0.1, accept_thres = NULL,
                               prob_fn = NULL,
                               support, spread_data,
                               centre = "local", sobol = TRUE,
                               phase = NULL) {

  ### Testo
  # data = wave1; n = 500; p_best = 0.15; prob_fn = NULL;
  # support = sup; spread_data = spread_data;
  # centre = "global"; sobol = TRUE;
  # pow = 1; var_factor = 1; accept_thres = NULL
  ### endo testo

  # if threshold is used, ignore p_best
  if(!is.null(accept_thres)) {
    p_best <- NULL
    n_best <- sum(data$overlap >= accept_thres)
  }

  # which particles will be resampled?
  if(!is.null(p_best)) {
    n_best <- round(nrow(data) * p_best)
    if(n_best < 500) {
      n_best <- 500
    }
  }
  if(n_best == "all") {
    n_best <- nrow(data)
  } else {
    data <- data[order(data$overlap, decreasing = TRUE), ]
  }

  # resample the best particles
  weight <- data$overlap

  # sometimes the best particles are too few, and the computation of the vcov
  # fails. In these cases we need to smooth a bit the weights, so more than 20
  # different particles are resampled.
  l <- 0
  k <- 0 # smoothing power (root)
  j <- 0
  while((l < 30) & (j < 100)) {
    ids_rep <- sample(1:n_best, size = n, replace = T,
                      prob = weight[1:n_best] ^ (1 / (1 + k)))
    l <- length(unique(ids_rep))
    k <- k + 0.5
    j <- j + 1
  }

  if(j == 100) {
    dexp <- data[1:n_best, ]
  } else {
    dexp <- data[ids_rep, ]
  }

  # transform to unconstrained scale
  dexp$par_values_raw <- logit_scaled(dexp$par_values, support)

  # MVN kernel
  # first, remove infinite values and redefine n
  finite <- apply(dexp$par_values_raw, 1, function(x) all(is.finite(x))) %>% t
  dexp <- dexp[finite, ]
  n <- nrow(dexp)

  V <- make_positive_definite(cov(dexp$par_values_raw))
  if(anyNA(V)) V <- diag(apply(dexp$par_values_raw, 2, var))
  Vlarge <- V * var_factor

  if(centre == "local") {
    mus <- dexp$par_values_raw # centred at good particles
    sobol <- FALSE
  }
  if(centre == "global") {
    mus <- apply(dexp$par_values_raw, 2, mean) # centred at global mean
  }

  candidates_raw <- rmvn_sobol(n, mus, Vlarge, sobol) %>% t
  # transform to original scale
  colnames(candidates_raw) <- par_names
  candidates <- invlogit_scaled(candidates_raw, support)
  colnames(candidates) <- par_names

  # message("Simulating fires")
  sim_result <- similarity_simulate_parallel(particles_mat = candidates,
                                             fire_data = spread_data)

  # define wave
  sim_result$wave <- max(data$wave) + 1
  sim_result$phase <- phase
  colnames(sim_result$par_values) <- par_names

  # merge old and new datasets
  res <- rbind(
    data[, colnames(sim_result)],
    sim_result
  )

  return(res)
}

# iterate waves of explore_likelihood()
explore_likelihood_iterate <- function(
    n_waves = 10,
    data, n = 500, var_factor = 1.5,
    n_best = "all", p_best = NULL, accept_thres = NULL,
    support, spread_data,
    centre = "global", sobol = TRUE,
    phase = NULL,
    write = FALSE, fire_name = NULL
) {

  last_wave <- max(data$wave)
  m <- paste("Search starts at max overlap = ",
             round(max(data$overlap), 4), sep = "")
  message(m)

  # initialize res as previous data
  res <- data

  # create var_factor sequence for the variable case
  var_factors <- rep(var_factor, ceiling(n_waves / length(var_factor)))

  # # create centre sequence for the variable case
  # each <- ceiling(length(var_factor) / 2) # var_factors are for each centre
  # centres_dup <- rep(centre, each = each)
  # centres <- rep(centres_dup, ceiling(n_waves / length(var_factor)))

  for(i in 1:n_waves) {
    res <- explore_likelihood(data = res, n = n, var_factor = var_factors[i],
                              n_best = n_best, p_best = p_best,
                              accept_thres = accept_thres,
                              support = support, spread_data = spread_data,
                              centre = centre, sobol = sobol,
                              phase = phase)

    m <- paste("Search wave ", last_wave + i, "; max overlap = ",
               round(max(res$overlap), 4), sep = "")

    if(write) {
      nn <- paste(fire_name, "-wave-", last_wave + i, ".rds", sep = "")
      saveRDS(res[res$wave == last_wave + i, ], file.path(target_dir, nn))
      # save only latest
    }

    message(m)
  }

  return(res)
}

# Function to perform adaptive MCMC in parallel by hand. It's just a wrapper
# for adaptMCMC::MCMC, because MCMC.parallel does not allow to fix the inits.
MCMC_parallel <- function(fun, n, adapt, scale, init_list, acc.rate = 0.234,
                          n_chains, n_cores, ...) {

  registerDoMC(n_cores)

  runs <- foreach(pp = init_list) %dopar% {
    MCMC(p = fun, n = n, adapt = adapt, scale = scale, init = pp,
         acc.rate = acc.rate, ...)
  }

  return(runs)
}


# make prob = 0 in particles with high uncertainty, measured as the width of the
# 95 % confidence interval (unc).
drop_uncertain <- function(unc, p, centre = 0.7, slope = -90) {
  return(p * plogis(slope * (unc - centre)))
}
# curve(plogis(-90 * (x - 0.6)), n = 300)

# sample from the GAM taking independent draws (simulation)
# unc_thres [0, 1] determines the maximum uncertainty allowed in the model.
#   if higher, the particle is rejected.
rejection_sample <- function(iter, model, support,
                             centre = 0.7, slope = -90) {

  ntry <- iter * 10
  got <- 0
  wave <- 0
  samples <- matrix(NA, ntry, n_coef)
  colnames(samples) <- par_names

  while(got < iter) {
    unif_draws <- matrix(runif(ntry * n_coef), ntry, n_coef)
    particles <- scale_params(unif_draws, support)
    colnames(particles) <- par_names
    part_df <- as.data.frame(particles)

    # compute probability
    prob <- predict(model, part_df, type = "response")

    # only evaluate uncertainty at particles with high probability.
    # (This step is time consuming)
    unc <- numeric(length(prob))
    high_prob_ids <- which(prob > 0.02)
    if(length(high_prob_ids) > 0) {
      pp <- predict(model, part_df[high_prob_ids, ], se.fit = T)
      lll <- plogis(pp$fit - qnorm(0.975) * pp$se.fit)
      uuu <- plogis(pp$fit + qnorm(0.975) * pp$se.fit)
      unc[high_prob_ids] <- uuu - lll
    }

    # drop probability if highly uncertain
    prob2 <- drop_uncertain(unc, prob, centre, slope)
    stay <- which(runif(ntry) <= prob2)

    l <- length(stay)
    if(l > 0) {
      samples[(got+1):(got+l), ] <- particles[stay, ]
    }

    wave <- wave + 1
    got <- got + l
    # print(paste("wave ", wave, ", got ", got, sep = ""))
  }

  samples <- samples[1:iter, ]
  return(samples)
}

# the same in parallel
rejection_sample_parallel <- function(iter, model, support,
                                      centre = 0.7, slope = -90,
                                      cores = 15) {
  registerDoMC(cores)
  ii <- as.list(1:cores)
  foreach(chain = ii) %dopar% {
    rejection_sample(iter, model, support, centre, slope)
  }
}

# sample from a region where the density is higher than a threshold
rejection_sample_region <- function(iter, kde_model, support, ll_thres) {

  ntry <- iter * 10
  got <- 0
  wave <- 0
  samples <- matrix(NA, ntry, n_coef)
  colnames(samples) <- par_names

  while(got < iter) {
    unif_draws <- matrix(runif(ntry * n_coef), ntry, n_coef)
    particles <- scale_params(unif_draws, support)
    colnames(particles) <- par_names

    # compute probability
    dens <- dkdevine(particles, kde_model)
    log_dens <- log(dens)

    high_prob_ids <- which(log_dens > ll_thres)
    l <- length(high_prob_ids)
    if(l > 0) {
      samples[(got+1):(got+l), ] <- particles[high_prob_ids, ]
    }

    wave <- wave + 1
    got <- got + l
    print(paste("wave ", wave, ", got ", got, sep = ""))
  }

  samples <- samples[1:iter, ]
  return(samples)
}

# sample from a region where the density is higher than a threshold
rejection_sample_region_parallel <- function(iter, kde_model, support, ll_thres,
                                             cores = 15) {
  registerDoMC(cores)
  ii <- as.list(1:cores)
  foreach(chain = ii) %dopar% {
    rejection_sample_region(iter, kde_model, support, ll_thres)
  }
}

# extract samples from a MCMC.parallel run and return a draws_array
# (posterior class). (Also useful for MCMC_parallel result.)
tidy_samples <- function(samples, adapt, support) {
  nc <- length(samples)
  rr <- nrow(samples[[1]]$samples) - adapt
  cc <- ncol(samples[[1]]$samples)
  arr0 <- array(NA, dim = c(rr, cc, nc))
  for(c in 1:nc) {
    arr0[, , c] <- invlogit_scaled(samples[[c]]$samples[-(1:adapt), ], support)
  }
  arr <- aperm(arr0, c(1, 3, 2))
  draws_arr <- as_draws_array(arr)

  draws_arr2 <- rename_variables(draws_arr,
                                 "forest" = "...1",
                                 "shrubland" = "...2",
                                 "grassland" = "...3",
                                 "ndvi" = "...4",
                                 "north" = "...5",
                                 "elev" = "...6",
                                 "slope" = "...7",
                                 "wind" = "...8",
                                 "steps" = "...9")

  return(draws_arr2)
}

# the same as tidy samples to be used with the result from
# rejection_sample_parallel
tidy_samples_ind <- function(samples_list) {
  nc <- length(samples_list)
  rr <- nrow(samples_list[[1]])
  cc <- ncol(samples_list[[1]])
  arr0 <- array(NA, dim = c(rr, cc, nc))
  for(c in 1:nc) {
    arr0[, , c] <- samples_list[[c]]
  }
  arr <- aperm(arr0, c(1, 3, 2))
  draws_arr <- as_draws_array(arr)

  draws_arr2 <- rename_variables(draws_arr,
                                 "forest" = "...1",
                                 "shrubland" = "...2",
                                 "grassland" = "...3",
                                 "ndvi" = "...4",
                                 "north" = "...5",
                                 "elev" = "...6",
                                 "slope" = "...7",
                                 "wind" = "...8",
                                 "steps" = "...9")

  return(draws_arr2)
}

# function to evaluate ABC-posterior, based on an overlap threshold
fn_like_sim_bin <- function(x, support = NULL, fire_data = NULL, thres = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  ov <- similarity_simulate_particle(xx, n_sim = n_rep, fire_data = fire_data) %>% mean
  if(ov >= thres) {
    return(sum(dlogis(x, log = TRUE)))
  } else {
    return(-Inf)
  }
}

# function to evaluate ABC-posterior, based on an overlap threshold
fn_like_sim_binom <- function(x, support = NULL, fire_data = NULL, thres = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  ov <- similarity_simulate_particle(xx, n_sim = n_rep, fire_data = fire_data)
  sum_above <- sum(ov >= thres)
  if(sum_above > 0) {
    log_post <- log(sum_above / n_rep) + sum(dlogis(x, log = TRUE))
    return(log_post)
  } else {
    return(-Inf)
  }
}

# function to reduce support, so the rejection sampling is faster.
reduce_support <- function(x, support, prop = 0.75) {
  ranges <- apply(x, 2, range)
  widths <- apply(ranges, 2, diff)

  rnew <- ranges

  rnew[1, ] <- ranges[1, ] - widths * prop
  rnew[2, ] <- ranges[2, ] + widths * prop

  # constrain to the original support
  rnew[1, ] <- ifelse(rnew[1, ] < support[1, ], support[1, ], rnew[1, ])
  rnew[2, ] <- ifelse(rnew[2, ] > support[2, ], support[2, ], rnew[2, ])

  return(rnew)
}

# Function to subset data to a given support. Returns the rows ids of the
# passing observations.
within_support <- function(x, support) {
  ## test
  # x <- wlist$like_sim$par_values
  # support <- support_reduced

  keep_mat <- matrix(NA, nrow(x), ncol(x))
  for(i in 1:ncol(x)) {
    keep_mat[, i] <- x[, i] >= support[1, i] & x[, i] <= support[2, i]
  }
  okrows <- apply(keep_mat, 1, all)
  return(which(okrows))
}

# Function to create formula based on the terms to include (names). See the
# gam_terms object above.
make_gam_formula <- function(terms) {
  # ## test
  # terms <- c("wet", "steps", "wind_steps")
  rhs <- paste(gam_terms[terms], collapse = " + ")
  ff <- as.formula(paste("y", rhs, sep = " ~ "))
  return(ff)
}

make_seml_formula <- function(terms) {
  ## test
  # terms <- c("wet", "steps", "wind")
  lhs_in <- paste(terms, collapse = ", ")
  lhs <- paste("cbind(", lhs_in, ")", sep = "")
  ff <- as.formula(paste(lhs, " ~ 1", sep = ""))
  return(ff)
}

within_support <- function(x, support) {
  evals <- sapply(1:ncol(x), function(j) {
    res <- x[, j] >= support[1, j] & x[, j] <= support[2, j]
  })
  out <- apply(evals, 1, all)
  return(out)
}

# Function to sample particles from a GAM in a compact support and retain only
# those with similarity above a given threshold, after running the simulator
# on them.
abc_gam <- function(gam, support_sampling, threshold, spread_data,
                    n = 600, n_cores_gam = 15, n_cores_fire = 15) {

  ## TESTO
  # gam = wlist$gam_bern; support_sampling =  support_sample
  # threshold = wlist$thres; spread_data = spread_data
  # n = 600; n_cores_gam = 15; n_cores_fire = n_cores_fire
  ## ENDO TESTO
  n_batch <- ceiling(n / n_cores_gam)

  # Simulate particles from the GAM
  samples_try <- rejection_sample_parallel(
    iter = n_batch, model = gam,
    support = support_sampling,
    centre = 0.5, slope = -90,
    cores = n_cores_gam
  )

  samples_try <- as.matrix(do.call("rbind", samples_try))

  registerDoMC(n_cores_fire)
  sim <- similarity_simulate_parallel(particles_mat = samples_try,
                                      fire_data = spread_data)

  # filter particles above threshold
  use <- sim$overlap >= threshold
  colnames(sim$par_values) <- par_names

  samples_out <- sim$par_values[use, ]
  return(samples_out)
}

# Simulation, 40000 particles ----------------------------------------------

# First, simulate the 48 smallest fires
last_fire <- size_data[size_data$size_burn_rel < 0.1, ] %>% nrow() # 48
which(size_data$fire_id == "1999_27j_N") # problems at fire 31

# (tested with fire 1, error at fire 31)
for(f in (last_fire+1):n_fires) {

  # f = 1

  fire_file <- size_data$file[f]
  fire_name <- size_data$fire_id[f]
  print(paste("Fire:", fire_name))
  full_data <- readRDS(file.path(data_dir, fire_file))
  # subset data needed for spread (to be cloned across workers)
  spread_data <- full_data[c("landscape", "ig_rowcol",
                             "burned_layer", "burned_ids")]

  bb <- get_bounds(fire_data = spread_data)
  sup <- rbind(params_lower, params_upper)
  sup[2, "steps"] <- bb["largest", "steps_used"]

  # Explore likelihood
  # Phase 1: sobol
  if(fire_name == "2015_50") {
    registerDoMC(13)
  } else {
    registerDoMC(15)
  }

  ss <- sobol(n = 10000, dim = n_coef)
  particles <- scale_params(ss, sup)
  waves_1 <- rep(1:20, each = 500)
  # loop to save
  wave1 <- NULL
  print(paste("Phase: sobol.", fire_name, sep = ""))
  for(w in 1:20) {
    # w = 1
    wlocal <- similarity_simulate_parallel(particles[waves_1 == w, ],
                                           spread_data)
    wlocal$wave <- w
    wlocal$phase <- "sobol"
    colnames(wlocal$par_values) <- par_names

    nn <- paste(fire_name, "-wave-", w, ".rds", sep = "")
    saveRDS(wlocal, file.path(target_dir, nn))

    wave1 <- rbind(wave1, wlocal)
    rr <- round(max(wave1$overlap), 4)
    print(paste("wave ", w, ", overlap max: ", rr, sep = ""))
  }

  # wave 2: reproduce 15 % best (500 * 20 = 10000)
  print(paste("Phase: search higher. Fire: ", fire_name, sep = ""))
  wave2 <- explore_likelihood_iterate(
    n_waves = 20, data = wave1, n = 500, p_best = 0.15,
    var_factor = c(1.5, 2, 1),
    support = sup, spread_data = spread_data,
    centre = "global", sobol = TRUE, phase = "search_higher",
    write = T, fire_name = fire_name
  )

  # 10000 iter reproducing the 100 best, to reach the max
  print(paste("Phase: search maximum. Fire: ", fire_name, sep = ""))
  wave3 <- explore_likelihood_iterate(
    n_waves = 100, data = wave2, n = 100, n_best = 100,
    var_factor = c(1.5, 2, 1),
    support = sup, spread_data = spread_data,
    centre = "global", sobol = TRUE, phase = "search_maximum",
    write = T, fire_name = fire_name
  )

  thres <- max(wave3$overlap) - 0.05
  # wave_plot(wave4, alpha = 0.1, thres = thres)

  # 10000 iter above thres
  print(paste("Phase: above threshold. Fire: ", fire_name, sep = ""))
  wave4 <- explore_likelihood_iterate(
    n_waves = 20, data = wave3, n = 500, accept_thres = thres,
    var_factor = c(1.5, 2, 1),
    support = sup, spread_data = spread_data,
    centre = "global", sobol = TRUE, phase = "above_threshold",
    write = T, fire_name = fire_name
  )

  like_sim <- wave4

  nn <- paste(fire_name, "-waveplot.png", sep = "")
  png(filename = file.path(target_dir, nn),
      width = 20, height = 11, units = "cm", res = 300)
  wave_plot(like_sim, alpha = 0.05,
            thres = thres, tit = fire_name, rc = c(2, 4))
  dev.off()

  # Fit gam

  print(paste("Fitting GAM. Fire: ", fire_name, sep = ""))

  data_gam_bern <- as.data.frame(cbind(like_sim$par_values,
                                       y = as.numeric(like_sim$overlap >= thres)))

  k_side <- 15
  k_int <- 6
  basis <- "cr"

  formula_use <- gam_formula_full

  gam_bern <- bam(
    formula_use, data = data_gam_bern, family = "binomial",
    method = "fREML", discrete = T, nthreads = 8
  )

  # Sample GAM-approximated posterior

  print(paste("Sampling GAM. Fire: ", fire_name, sep = ""))

  sup_reduced <- reduce_support(like_sim$par_values[like_sim$overlap >= thres, ],
                                sup, 0.5)
  r_gam <- rejection_sample_parallel(700, gam_bern, sup_reduced,
                                     centre = 0.7,
                                     cores = 15)
  draws <- do.call("rbind", r_gam) %>% as.data.frame

  # Check GAM

  ids <- sample(1:nrow(draws), size = 2000, replace = F)
  ppmat <- as.matrix(draws[ids, ])

  print(paste("Checking GAM. Fire: ", fire_name, sep = ""))

  overlap_check <- similarity_simulate_parallel(particles = ppmat,
                                                fire_data = spread_data)

  nn <- paste(fire_name, "-gamcheck.png", sep = "")
  png(filename = file.path(target_dir, nn),
      width = 10, height = 10, units = "cm", res = 300)
  hist(overlap_check$overlap, xlim = c(0, 1), main = fire_name,
       xlab = "Overlap", ylab = "Relative frequency", freq = F)
  abline(v = thres, col = 2, lwd = 2)
  abline(v = mean(overlap_check$overlap), col = 4, lwd = 2)
  dev.off()

  # save result
  res <- list(
    fire_name = fire_name,
    veg_props = round(normalize(full_data$cells_by_veg) * 100, 2),
    like_sim = like_sim,             # first 40000 simulations
    thres = thres,                   # overlap threshold
    gam_bern = gam_bern,             # bernoulli GAM
    draws = draws,                   # draws sampled from the GAM
    overlap_check = overlap_check    # subsample from draws for which the overlap
                                     #   was computed
  )

  nn <- paste(fire_name, "-simulations_list.rds", sep = "")
  saveRDS(res, file.path(target_dir, nn))

  # remove temporal files
  all_files_saved <- list.files(target_dir, pattern = "-wave-")
  invisible(lapply(all_files_saved, function(f) unlink(file.path(target_dir, f))))

  # clean
  remove(full_data, spread_data, wave1, wave2, wave3, wave4, like_sim)
  gc()
}



# Fire 31 = 1999_27j_N (Falso Granito) had problems, it got stuck twice when
# sampling above threshold. It has a very peaky function, and the GAM has serious
# trouble to fit it.

# Took ~ 1.3 days to run the 48 smallest fires, 40000 + 2000 particles
# Intermediate file saving was disabled for small fires, but enabled for large
# ones.
# Start run of large fires at 18/07/2024 00:56

# In total it took ~ 4.3 days

# Export pair-plots --------------------------------------------------------

# Visualize plots to decide which terms to include in GAMs. Write that in the
# gam_terms_to_include.ods file.

simfiles <- list.files(target_dir, pattern = "-simulations_list.rds")

xvars <- list(
  rep("slope", n_veg),
  rep("wind", n_veg),
  c("wind", "steps", "steps")
)

yvars <- list(
  veg_names,
  veg_names,
  c("slope", "slope", "wind")
)

n_cols <- sapply(yvars, length)

for(f in 1:n_fires) {
  # f = 51
  wlist <- readRDS(file.path(target_dir, simfiles[f]))

  parvalues <- rbind(wlist$like_sim$par_values,
                     wlist$overlap_check$par_values)
  ovv <- c(wlist$like_sim$overlap, wlist$overlap_check$overlap)

  ranges <- apply(parvalues, 2, range)
  dd <- parvalues[ovv >= wlist$thres, ]

  # waveplot including test data
  nn <- paste(wlist$fire_name, "-waveplot_2.png", sep = "")
  png(filename = file.path(target_dir, nn),
      width = 22, height = 11, units = "cm", res = 300)

  wave_plot(wlist$like_sim, alpha = 0.05, thres = wlist$thres,
            rc = c(2, 5), tit = wlist$fire_name)
  dev.off()

  # pairs plot
  nn <- paste(wlist$fire_name, "-waveplot_pairs.png", sep = "")
  png(filename = file.path(target_dir, nn),
      width = 22, height = 15, units = "cm", res = 300)

  # start plot
  par(mfrow = c(3, 5))
  for(r in 1:3) {
    for(c in 1:n_cols[r]) {
      # r = 1; c = 1

      tit <- ifelse(r == 1 & c == 1, wlist$fire_name, NA)

      yname <- yvars[[r]][c]
      xname <- xvars[[r]][c]

      plot(dd[, yname] ~ dd[, xname], xlab = xname, ylab = yname,
           xlim = ranges[, xname], ylim = ranges[, yname],
           pch = 19, col = rgb(0, 0, 0, alpha),
           main = tit)
    }
  }
  par(mfrow = c(1, 1))
  # end plot

  dev.off()
}

# The problem was mostly related to low data availability, not so much
# to GAMs with excessive complexity.
# It was not low k either; GAMs with low k are able to fit very sharp surfaces.

# Get more particles in the interesting area ------------------------------

# run 10000 sobol particles in the square of high overlap based on the data
# obtained up to the moment.
# Then, fit a GAM, and sample 10000 particles. Evaluate and refit.
# Sample 1000 more and check.

# evaluate which require improvement
drefit <- readxl::read_xls(file.path("data", "focal fires data",
                           "gam_terms_to_include.xls"),
                           sheet = 3)

# Data to simplify GAMs
dterms <- readxl::read_xls(file.path("data", "focal fires data",
                                     "gam_terms_to_include.xls"),
                           sheet = 2) %>% as.data.frame()
rownames(dterms) <- dterms$fire_id
dterms <- as.matrix(dterms[, -1])
dterms[is.na(dterms)] <- 0

# barplot(table(drefit$refit))
# View(drefit)

# probar con algunos:
# fires_try <- c("1999_25j", "1999_28", "2014_1", "2014_3",  "2015_11")
# ids_try <- which(size_data$fire_id %in% fires_try)

fires_improve <- drefit$fire_id[drefit$refit > 0]
ids_improve <- which(size_data$fire_id %in% fires_improve)

for(f in ids_improve[-1]) { #1:n_fires) { # al 1 lo hice por separado.

  # f = 2
  fire_file <- size_data$file[f]
  fire_name <- size_data$fire_id[f]
  print(paste("Fire:", fire_name))
  full_data <- readRDS(file.path(data_dir, fire_file))
  # subset data needed for spread (to be cloned across workers)
  spread_data <- full_data[c("landscape", "ig_rowcol",
                             "burned_layer", "burned_ids")]

  bb <- get_bounds(fire_data = spread_data)
  sup <- rbind(params_lower, params_upper)
  sup[2, "steps"] <- bb["largest", "steps_used"]

  # import previous simulations
  fn <- paste(fire_name, "-simulations_list.rds", sep = "")
  wlist <- readRDS(file.path(target_dir, fn))
  thres <- wlist$thres

  # merge all previous simulations
  wlist$overlap_check$phase <- "gam_check_1"

  like_sim2 <- rbind(
    wlist$like_sim,
    wlist$overlap_check
  )
  above <- like_sim2$overlap >= thres

  # reduce support
  sup_reduced <- reduce_support(like_sim2$par_values[above, ],
                                sup, prop = 0.2)

  # Throw particles in the reduced support.
  if(fire_name == "2015_50") {
    registerDoMC(13)
  } else {
    registerDoMC(15)
  }

  ss <- sobol(n = 10000, dim = n_coef)
  particles <- scale_params(ss, sup_reduced)
  waves_1 <- rep(1:10, each = 1000)
  # loop to save
  wave1 <- NULL
  print(paste("Phase: sobol. Fire: ", fire_name, sep = ""))
  for(w in 1:10) {
    # w = 1
    wlocal <- similarity_simulate_parallel(particles[waves_1 == w, ],
                                           spread_data)
    wlocal$wave <- w + max(like_sim2$wave)
    wlocal$phase <- "sobol_2"
    colnames(wlocal$par_values) <- par_names

    nn <- paste(fire_name, "-wave-", w, ".rds", sep = "")
    saveRDS(wlocal, file.path(target_dir, nn))

    wave1 <- rbind(wave1, wlocal)
    rr <- round(max(wave1$overlap), 4)
    print(paste("wave ", w, ", overlap max: ", rr, sep = ""))
  }

  # Merge old and new data (including data used for testing)
  like_sim3 <- rbind(
    like_sim2,
    wave1[, colnames(like_sim2)]
  )

  # Fit gam
  print(paste("Fitting GAM. Fire: ", fire_name, sep = ""))

  data_gam_bern <- as.data.frame(
    cbind(
      like_sim3$par_values,
      y = as.numeric(like_sim3$overlap >= thres)
    )
  )

  k_side <- 15
  k_int <- 7
  basis <- "cr"

  terms_include <- terms_names[dterms[fire_name, ] == 1]
  formula_use <- make_gam_formula(terms_include)

  gam_bern <- bam(
    formula_use, data = data_gam_bern, family = "binomial",
    method = "fREML", discrete = T, nthreads = 8
  )

  # Sample GAM-approximated posterior

  print(paste("Sampling GAM. Fire: ", fire_name, sep = ""))

  sup_gam <- reduce_support(like_sim3$par_values[like_sim3$overlap >= thres, ],
                            sup, 0.5)
  r_gam <- rejection_sample_parallel(670, gam_bern, sup_gam,
                                     centre = 0.5,
                                     cores = 15)
  draws <- do.call("rbind", r_gam) %>% as.data.frame

  # Simulate GAM-based samples
  ppmat <- as.matrix(draws)

  print(paste("Simulating GAM samples. Fire: ", fire_name, sep = ""))

  wave_ids <- rep(1:15, each = 670)[1:nrow(draws)]
  # loop to save
  wave_gam <- NULL
  print(paste("Phase: sim_gam. Fire: ", fire_name, sep = ""))

  if(fire_name == "2015_50") {
    registerDoMC(13)
  } else {
    registerDoMC(15)
  }

  for(w in 1:15) {
    # w = 1
    wlocal <- similarity_simulate_parallel(ppmat[wave_ids == w, ],
                                           spread_data)
    wlocal$wave <- w + max(like_sim3$wave)
    wlocal$phase <- "sim_gam"
    colnames(wlocal$par_values) <- par_names

    nn <- paste(fire_name, "-wave-", w, ".rds", sep = "")
    saveRDS(wlocal, file.path(target_dir, nn))

    wave_gam <- rbind(wave_gam, wlocal)
    print(w)
  }

  # Merge data
  like_sim4 <- rbind(
    like_sim3,
    wave_gam[, colnames(like_sim3)]
  )

  # Refit GAM
  print(paste("Re-fitting GAM. Fire: ", fire_name, sep = ""))

  data_gam_bern <- as.data.frame(
    cbind(
      like_sim4$par_values,
      y = as.numeric(like_sim4$overlap >= thres)
    )
  )

  gam_bern <- bam(
    formula_use, data = data_gam_bern, family = "binomial",
    method = "fREML", discrete = T, nthreads = 8
  )

  # Sample GAM-approximated posterior again (test)
  print(paste("Testing GAM. Fire: ", fire_name, sep = ""))
  sup_gam <- reduce_support(like_sim4$par_values[like_sim4$overlap >= thres, ],
                            sup, 0.5)
  r_gam <- rejection_sample_parallel(100, gam_bern, sup_gam,
                                     centre = 0.5,
                                     cores = 15)
  draws <- do.call("rbind", r_gam) %>% as.data.frame

  if(fire_name == "2015_50") {
    registerDoMC(13)
  } else {
    registerDoMC(15)
  }

  wgam <- similarity_simulate_parallel(particles = as.matrix(draws),
                                       fire_data = spread_data)
  wgam$wave <- max(like_sim4$wave) + 1
  wgam$phase <- "gam_check_2"

  nn <- paste(fire_name, "-gamcheck2.png", sep = "")
  png(filename = file.path(target_dir, nn),
      width = 10, height = 10, units = "cm", res = 300)
  hist(wgam$overlap, xlim = c(0, 1), main = fire_name,
       xlab = "Overlap", ylab = "Relative frequency", freq = F)
  abline(v = thres, col = 2, lwd = 2)
  abline(v = mean(wgam$overlap), col = 4, lwd = 2)
  dev.off()

  # save result
  res <- list(
    fire_name = fire_name,
    veg_props = round(normalize(full_data$cells_by_veg) * 100, 2),
    like_sim = like_sim4,            # simulations
    thres = thres,                   # overlap threshold
    gam_bern = gam_bern,             # bernoulli GAM
    draws = draws,                   # draws sampled from the GAM
    overlap_check = wgam,            # tested GAM samples
    support = sup
  )

  nn <- paste(fire_name, "-simulations_list2.rds", sep = "")
  saveRDS(res, file.path(target_dir, nn))

  # remove temporal files
  all_files_saved <- list.files(target_dir, pattern = "-wave-")
  invisible(lapply(all_files_saved, function(f) unlink(file.path(target_dir, f))))

  # clean
  remove(full_data, spread_data, like_sim4)
  gc()
}



# Assessing bad GAMs ------------------------------------------------------

dterms <- readxl::read_xls(file.path("data", "focal fires data",
                                     "gam_terms_to_include.xls"),
                           sheet = 2) %>% as.data.frame()
rownames(dterms) <- dterms$fire_id
dterms <- as.matrix(dterms[, -1])
dterms[is.na(dterms)] <- 0
terms_names

# Table indicating which to refit
dcheck <- readxl::read_xls(file.path("data", "focal fires data",
                                     "gam_terms_to_include.xls"),
                           sheet = 3) %>% as.data.frame()
rownames(dcheck) <- dcheck$fire_id
dcheck <- as.matrix(dcheck[, -1])
dcheck[is.na(dcheck)] <- 0
dcheck[dcheck[, "refit2"] == 1, ]



# Fixing GAM 1999_27j_N ---------------------------------------------------

fire_name <- "1999_27j_N"
wlist <- readRDS(file.path(target_dir,
                           paste(fire_name, "-simulations_list2.rds", sep = "")))
fire_file <- size_data$file[size_data$fire_id == fire_name]
full_data <- readRDS(file.path(data_dir, fire_file))
spread_data <- full_data[c("landscape", "ig_rowcol",
                           "burned_layer", "burned_ids")]

x11()
wave_plot_pairs(wlist$like_sim, thres = wlist$thres)

dterms[fire_name, ]

# refitted GAM with higher k but (didn't work) (marg = 20, int = 15)
# refit in a smaller space: reduce_support(prop = 1)

rows_above <- wlist$like_sim$overlap >= wlist$thres
support_fit <- reduce_support(wlist$like_sim$par_values[rows_above, ],
                              support = wlist$support,
                              prop = 1)
rows_in_fit <- within_support(wlist$like_sim$par_values, support_fit)

formula_simple <- make_gam_formula(terms_names[dterms[fire_name, ] == 1])
k_side <- 15; k_int <- 7
basis <- "cr"

data_gam_bern <- as.data.frame(
  cbind(
    wlist$like_sim$par_values[rows_in_fit, ],
    y = as.numeric(wlist$like_sim$overlap[rows_in_fit] >= wlist$thres)
  )
)

gam_bern <- bam(
  formula_simple, data = data_gam_bern, family = "binomial",
  method = "fREML", discrete = T, nthreads = 8
)

# Sample GAM-approximated posterior again (test)
fff <- wlist$like_sim$overlap >= wlist$thres
sup_gam <- reduce_support(wlist$like_sim$par_values[fff, ],
                          wlist$support, 0.5)
r_gam <- rejection_sample_parallel(100, gam_bern, sup_gam,
                                   centre = 0.5,
                                   cores = 15)
draws <- do.call("rbind", r_gam) %>% as.data.frame

registerDoMC(15)
wgam <- similarity_simulate_parallel(particles = as.matrix(draws),
                                     fire_data = spread_data)
wgam$wave <- max(wlist$like_sim$wave) + 1
wgam$phase <- "gam_check_2"

nn <- paste(fire_name, "-gamcheck3.png", sep = "")
png(filename = file.path(target_dir, nn),
    width = 10, height = 10, units = "cm", res = 300)
hist(wgam$overlap, xlim = c(0, 1), main = fire_name,
     xlab = "Overlap", ylab = "Relative frequency", freq = F)
abline(v = wlist$thres, col = 2, lwd = 2)
abline(v = mean(wgam$overlap), col = 4, lwd = 2)
dev.off()
sum(wgam$overlap >= wlist$thres) / nrow(wgam)


x11()
wave_plot_pairs(wgam, thres = wlist$thres)
# Could not fix it, use first gam.

# Fixing GAM 2021_2146405150_E ---------------------------------------------------

fire_name <- "2021_2146405150_E"
wlist <- readRDS(file.path(target_dir,
                           paste(fire_name, "-simulations_list2.rds", sep = "")))
fire_file <- size_data$file[size_data$fire_id == fire_name]
full_data <- readRDS(file.path(data_dir, fire_file))
spread_data <- full_data[c("landscape", "ig_rowcol",
                           "burned_layer", "burned_ids")]

x11()
wave_plot_pairs(wlist$like_sim, thres = wlist$thres)

dterms[fire_name, ]

# remove some interactions

terms_less <- c("dry", "shrubland", "grassland",
                "shrubland_wind", "grassland_wind",
                "slope_wind")
formula_simple <- make_gam_formula(terms_less)
k_side <- 15; k_int <- 7
basis <- "cr"

data_gam_bern <- as.data.frame(
  cbind(
    wlist$like_sim$par_values,
    y = as.numeric(wlist$like_sim$overlap >= wlist$thres)
  )
)

gam_bern <- bam(
  formula_simple, data = data_gam_bern, family = "binomial",
  method = "fREML", discrete = T, nthreads = 8
)

# Sample GAM-approximated posterior again (test)
fff <- wlist$like_sim$overlap >= wlist$thres
sup_gam <- reduce_support(wlist$like_sim$par_values[fff, ],
                          wlist$support, 0.5)
r_gam <- rejection_sample_parallel(100, gam_bern, sup_gam,
                                   centre = 0.5,
                                   cores = 15)
draws <- do.call("rbind", r_gam) %>% as.data.frame

registerDoMC(15)
wgam <- similarity_simulate_parallel(particles = as.matrix(draws),
                                     fire_data = spread_data)
wgam$wave <- max(wlist$like_sim$wave) + 1
wgam$phase <- "gam_check_2"

nn <- paste(fire_name, "-gamcheck3.png", sep = "")
png(filename = file.path(target_dir, nn),
    width = 10, height = 10, units = "cm", res = 300)
hist(wgam$overlap, xlim = c(0, 1), main = fire_name,
     xlab = "Overlap", ylab = "Relative frequency", freq = F)
abline(v = wlist$thres, col = 2, lwd = 2)
abline(v = mean(wgam$overlap), col = 4, lwd = 2)
dev.off()
sum(wgam$overlap >= wlist$thres) / nrow(wgam)


x11()
wave_plot_pairs(wgam, thres = wlist$thres)

# Can't fix it, use the first one.

# Fixing GAM 2021_865 ---------------------------------------------------

fire_name <- "2021_865"
wlist1 <- readRDS(file.path(target_dir,
                            paste(fire_name, "-simulations_list.rds", sep = "")))
wlist2 <- readRDS(file.path(target_dir,
                            paste(fire_name, "-simulations_list2.rds", sep = "")))


sum(wlist1$overlap_check$overlap >= wlist1$thres) / nrow(wlist1$overlap_check)
sum(wlist2$overlap_check$overlap >= wlist2$thres) / nrow(wlist2$overlap_check)



# Proportion of good particles sampled from GAM ---------------------------

dcheck <- readxl::read_xls(file.path("data", "focal fires data",
                                     "gam_terms_to_include.xls"),
                           sheet = 3) %>% as.data.frame()
dcheck[is.na(dcheck)] <- 0
dcheck$prop_above <- NA

for(f in 1:nrow(df_check)) {
  # f = 55
  fire_name <- rownames(df_check)[f]
  data_type <- ifelse(df_check$refit[f] == 1 & df_check$use_first[f] == 0,
                      "-simulations_list2.rds",
                      "-simulations_list.rds")

  wlist <- readRDS(file.path(
    target_dir, paste(fire_name, data_type, sep = "")
  ))

  pp <- sum(wlist$overlap_check$overlap >= wlist$thres) / nrow(wlist$overlap_check)
  dcheck$prop_above[f] <- pp
}
# hist(dcheck$prop_above)

# Sample GAM-based posteriors ---------------------------------------------

# Load check table to know which GAM to use.
dcheck <- readxl::read_xls(file.path("data", "focal fires data",
                                     "gam_terms_to_include.xls"),
                           sheet = 3) %>% as.data.frame()
dcheck[is.na(dcheck)] <- 0

for(f in 1:nrow(size_data)) {
  # f = 31

  fire_file <- size_data$file[f]
  fire_name <- size_data$fire_id[f]
  full_data <- readRDS(file.path(data_dir, fire_file))
  # subset data needed for spread (to be cloned across workers)
  spread_data <- full_data[c("landscape", "ig_rowcol",
                             "burned_layer", "burned_ids")]

  # load list with GAM, simulations, and threshold
  data_type <- ifelse(dcheck$refit[dcheck$fire_id == fire_name] != 0 &
                        dcheck$use_first[dcheck$fire_id == fire_name] == 0,
                      "-simulations_list2.rds",
                      "-simulations_list.rds")

  wlist <- readRDS(file.path(
    target_dir, paste(fire_name, data_type, sep = "")
  ))

  # proportion of accepted particles
  pp <- sum(wlist$overlap_check$overlap >= wlist$thres) / nrow(wlist$overlap_check)

  mm <- paste("Fire: ", fire_name, "; accept probability = ",
              round(pp * 100, 2), " %", sep = "")
  message(mm)

  # Get support to sample
  above <- wlist$like_sim$overlap >= wlist$thres
  support_original <- rbind(params_lower, params_upper)
  support_original[2, "steps"] <- max(wlist$like_sim$par_values[, "steps"])
  support_sample <- reduce_support(wlist$like_sim$par_values[above, ],
                                   support_original, prop = 0.25)

  n_cores_fire <- ifelse(fire_name == "2015_50", 13, 15)

  # Sample in batches of 600 particles to try
  samples_out <- NULL
  got <- 0
  wave <- 0

  samples_out <- readRDS(file.path(target_dir, "1999_27j_N -posterior_samples.rds"))
  got <- nrow(samples_out)
  wave <- 51

  while(got < 5000) {

    samples_i <- abc_gam(
      gam = wlist$gam_bern, support_sampling =  support_sample,
      threshold = wlist$thres, spread_data = spread_data,
      n = 600, n_cores_gam = 15, n_cores_fire = n_cores_fire
    )

    wave <- wave + 1

    if(!is.null(samples_i)) {
      samples_out <- rbind(samples_out, samples_i)
      got <- got + nrow(samples_i)

      temp_file_name <- paste(fire_name, "-sampling_temp_wave_", wave, ".rds",
                              sep = "")
      saveRDS(samples_i, file.path(target_dir, temp_file_name))
    }

    mm <- paste("Wave ", wave, "; got ", got, sep = "")
    message(mm)
  }

  # add fire_id, support_original and threshold as attributes
  attr(samples_out, "fire_id") <- fire_name
  attr(samples_out, "support") <- support_original
  attr(samples_out, "threshold") <- wlist$thres

  # write final matrix of samples
  fn <- paste(fire_name, "-posterior_samples.rds")
  saveRDS(samples_out, file.path(target_dir, fn))

  # remove temporal files
  all_files_saved <- list.files(target_dir, pattern = "-sampling_temp_wave_")
  invisible(lapply(all_files_saved, function(f) unlink(file.path(target_dir, f))))

  # clean
  remove(full_data, spread_data, wlist)
  gc()
}
# done?
