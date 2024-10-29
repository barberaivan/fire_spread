# This code inherits from
# <fire_wise_posteriors_FI.R>.

# I don't trust the previous exploration, as it used the abc_prob instead of
# overlap, which is very peaky and problematic if an expquad kernel is to be
# used later. From the overlap POV, many waveplots are quite rare.
# Perhaps it is safer to explore the parameter space again, as follows:

# 10000 sobol
# 10000 reproducing 15 % best
# 10000 reproducing 1000 best
# 10000 reproducing 100 best

# with n_rep = 1, it should take ~ 8 h


# (the fire simulator is edited from here on, so it doesn't compute
# the non-directional layer before simulation, as that only makes sense if the
# simulation is to be repeated many times for the same particle)

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
target_dir <- file.path("files", "posterior_samples_stage1_exploration")

# load file with constants to standardize
fi_params <- readRDS(file.path("data", "flammability indices",
                               "flammability_indices.rds"))
slope_sd <- fi_params$slope_term_sd

# constants for fire spread simulation
upper_limit <- 1
n_veg <- 5
veg_names <- c("wet", "subalpine", "dry", "shrubland", "grassland")
n_terrain <- 2
terrain_names <- c("slope", "wind")
terrain_variables <- c("elevation", "wdir", "wspeed")
n_nd <- n_fi <- 2        # flammability indices
nd_variables <- c("vfi", "tfi")

par_names <- c("intercept", nd_variables, terrain_names, "steps")
n_coef <- length(par_names)

# number of fires to simulate by particle
n_rep <- 1

ext_alpha <- 50
ext_beta <- 30

params_lower <- c(-ext_alpha, rep(0, n_coef-2), 5)
params_upper <- c(ext_alpha, rep(ext_beta, n_coef-2), NA)
names(params_lower) <- names(params_upper) <- par_names
params_upper["slope"] <- ext_beta / slope_sd

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

## CAREFUL: below another variable was added to size_data

# Functions ---------------------------------------------------------------

normalize <- function(x) x / sum(x)

# compute scale for a gaussian kernel where y = p * pmax at delta_crit
get_scale <- function(delta_crit, p) {
  ss <- delta_crit / sqrt(-log(p))
  return(ss)
}

xcrit <- 0.5; p <- 0.5
kscale <- get_scale(xcrit, p)

# Compute acceptance probability with an exponential quadratic function,
# based on scaled overlap.
kernel_expquad <- function(overlap, overlap_max,
                           ov_scale = 0.05, log = F) {

  # Scale overlap to have maximum = 1
  ov_unit <- overlap / overlap_max
  ov_unit[ov_unit > 1] <- 1
  ov_dist <- 1 - ov_unit

  out <- dnorm(ov_dist, sd = ov_scale, log = log)

  return(out)
}

# Simpler exponential quadratic kernel (not normalized)
expquad <- function(x, sigma, log = F) {
  ylog <- -(x / sigma) ^ 2
  if(log) return(ylog) else return(exp(ylog))
}

expquad_ov <- function(overlap, overlap_max, sigma, log = F) {
  ov_unit <- overlap / overlap_max
  ov_unit[ov_unit > 1] <- 1
  x <- 1 - ov_unit # turn into distance
  ylog <- -(x / sigma) ^ 2

  if(log) return(ylog) else return(exp(ylog))
}

expquad_flat <- function(overlap, overlap_high,
                         sigma = 0.01, log = F) {

  ov_unit <- overlap / overlap_high
  ov_unit[ov_unit > 1] <- 1
  x <- 1 - ov_unit # turn into distance
  ylog <- -(x / sigma) ^ 2

  if(log) return(ylog) else return(exp(ylog))
}


# Simulate fires and compare then with the observed one using the overlap_spatial
# function. The landscape argument includes all data to simulate fire, and also
# burned and burned_ids layers.
# Returns a matrix with the mean and variance across replicates of:
#   overlap_spatial,
#   size_diff: size differences, as simulated - observed,
#   edge: number of pixels burned at the edge of the landscape.
# The last two are used to rule out bounded fires.
similarity_simulate_particle <- function(particle, fire_data = NULL) {

  coef_intercepts <- rep(particle[1], n_veg)
  coef_nd <- particle[2:3]
  coef_terrain <- particle[4:(n_coef-1)]
  steps <- particle[n_coef]

  fire_sim <- simulate_fire_compare(
    layer_vegetation = fire_data$landscape[, , "veg"],
    layer_nd = fire_data$landscape[, , nd_variables],
    layer_terrain = fire_data$landscape[, , terrain_variables],
    coef_intercepts = coef_intercepts,
    coef_nd = coef_nd,
    coef_terrain = coef_terrain,
    ignition_cells = fire_data$ig_rowcol,
    upper_limit = upper_limit,
    steps = steps
  )

  ov <- overlap_spatial(
    fire_sim[c("burned_layer", "burned_ids")],
    fire_data[c("burned_layer", "burned_ids")]
  )

  return(ov)
}

# The same but returning many metrics. Currently used only for steps used.
# It returns a data.frame with averages.
similarity_simulate_particle_metrics <- function(particle, fire_data = NULL) {

  metrics <- numeric(3)
  names(metrics) <- c("overlap", "size_diff", "steps_used")

  coef_intercepts <- rep(particle[1], n_veg)
  coef_nd <- particle[2:3]
  coef_terrain <- particle[4:(n_coef-1)]
  steps <- particle[n_coef]

  fire_sim <- simulate_fire_compare(
    layer_vegetation = fire_data$landscape[, , "veg"],
    layer_nd = fire_data$landscape[, , nd_variables],
    layer_terrain = fire_data$landscape[, , terrain_variables],
    coef_intercepts = coef_intercepts,
    coef_nd = coef_nd,
    coef_terrain = coef_terrain,
    ignition_cells = fire_data$ig_rowcol,
    upper_limit = upper_limit,
    steps = steps
  )

  metrics["overlap"] <- overlap_spatial(
    fire_sim, fire_data[c("burned_layer", "burned_ids")]
  )

  metrics["size_diff"] <-
    ncol(fire_sim$burned_ids) -
    ncol(fire_data$burned_ids)

  metrics["steps_used"] <- fire_sim$steps_used

  return(metrics)
}

# Emulate the loglik over a list of particles in parallel. Returns the metrics
# in rows.
# It requires the matrix of particles to simulate (coefficient values).
similarity_simulate_parallel <- function(particles_mat = NULL,
                                         fire_data = NULL) {

  # turn particle matrix into list for parallel evaluation
  particles_list <- lapply(1:nrow(particles_mat),
                           function(x) particles_mat[x, ])

  # simulate in parallel
  result <- foreach(pp = particles_list) %dopar% {
    similarity_simulate_particle(pp, fire_data)
  }

  # rbind list result
  ov <- unlist(result)

  # return data.frame
  res <- data.frame(
    wave = NA,
    overlap = ov
  )
  res$par_values <- particles_mat
  res$phase <- NA
  return(res)
}

# get_bounds: computes the metrics for the largest and smallest fires possible
get_bounds <- function(fire_data) {

  coef_burn_all <- c(1e6, rep(0, n_coef - 1)) # steps = 0 means infinite
  coef_burn_none <- c(-1e6, rep(0, n_coef - 2), 1)

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
                      tit = NULL, rc = c(2, 3)) {

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
  # check NA or inf before definiteness
  if(anyNA(m) | !all(is.finite(m))) {
    return(NA)
  }

  if(is_positive_definite(m)) {
    return(m)
  } else {
    temp <- tryNULL(nearPD(m, corr = corr, keepDiag = TRUE, maxit = max_iter))
    if(is.null(temp)) return(NA) else {
      mat <- as.matrix(temp$mat)

      # check NA or inf before definiteness
      if(anyNA(mat) | !all(is.finite(mat))) {
        return(NA)
      }

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
                    n = 600, n_cores_gam = 15, n_cores_fire = 15,
                    n_rep = 20) {

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
                                      fire_data = spread_data,
                                      n_sim = n_rep)

  # filter particles above threshold
  use <- sim$overlap >= threshold
  colnames(sim$par_values) <- par_names

  samples_out <- sim$par_values[use, ]
  return(samples_out)
}

selm_robust <- function(data){
  tryCatch(
    # try to fit smn
    {
      snfit <- selm(
        formula = selm_form,
        data = data,
        family = "SN"
      )

      dp <- coef(snfit, "DP", vector = FALSE)
      names(dp) <- c("xi", "Omega", "alpha")
      dp$xi <- as.numeric(dp$xi)
      names(dp$xi) <- par_names

      return(dp)
    },
    # return empirical MVN if fails
    error = function(e) {
      message("SN fit failed")

      V0 <- cov(data)
      V0[lower.tri(V0)] = t(V0)[lower.tri(V0)]

      V <- make_positive_definite(V0)
      if(anyNA(V)) V <- diag(apply(data, 2, var))

      dp = list(
        xi = colMeans(data),
        Omega = V,
        alpha = rep(0, n_coef)
      )
      return(dp)
    }
  )
}

# Get steps bounds for all fires ------------------------------------------

# size_data$steps_upper <- NA
# for(f in 1:n_fires) {
#   print(f)
#   fire_file <- size_data$file[f]
#   fire_name <- size_data$fire_id[f]
#   full_data <- readRDS(file.path(data_dir, fire_file))
#   # subset data needed for spread (to be cloned across workers)
#   spread_data <- full_data[c("landscape", "ig_rowcol",
#                              "burned_layer", "burned_ids",
#                              "counts_veg")]
#   bb <- get_bounds(fire_data = spread_data)
#   size_data$steps_upper[f] <- bb["largest", "steps_used"]
# }
# saveRDS(size_data, file.path("data", "focal fires data", "fire_size_data.rds"))
size_data <- readRDS(file.path("data", "focal fires data", "fire_size_data.rds"))

# Simulate fires ----------------------------------------------------------

for(f in 2:n_fires) {
  write_file <- ifelse(f > 30, TRUE, FALSE)

  fire_file <- size_data$file[f]
  fire_name <- size_data$fire_id[f]
  print(paste("Fire:", fire_name))
  full_data <- readRDS(file.path(data_dir, fire_file))
  # subset data needed for spread (to be cloned across workers)
  spread_data <- full_data[c("landscape", "ig_rowcol",
                             "burned_layer", "burned_ids",
                             "counts_veg")]

  bb <- get_bounds(fire_data = spread_data)
  sup <- rbind(params_lower, params_upper)
  sup[2, "steps"] <- bb["largest", "steps_used"]

  if(fire_name %in% c("2015_50", "1999_2140469994_r")) {
    registerDoMC(10)
  } else {
    registerDoMC(15)
  }

  # Explore likelihood

  # Phase 1: sobol
  ss <- sobol(n = 10000, dim = n_coef)
  particles <- scale_params(ss, sup)
  wave1 <- similarity_simulate_parallel(particles, spread_data)
  wave1$wave <- 1
  wave1$phase <- "sobol"

  if(write_file) {
    nn <- paste(fire_name, "-wave1.rds", sep = "")
    saveRDS(wave1, file.path(target_dir, nn))
  }

  rr <- round(max(wave1$overlap), 4)
  print(paste("wave 1; overlap max: ", rr, sep = ""))

  # wave 2: reproduce 15 % best (500 * 20 = 10000)
  print(paste("Phase: search higher. Fire: ", fire_name, sep = ""))
  wave2 <- explore_likelihood_iterate(
    n_waves = 20, data = wave1, n = 500, p_best = 0.15,
    var_factor = c(1.5, 2, 1),
    support = sup, spread_data = spread_data,
    centre = "global", sobol = TRUE, phase = "search_higher",
    write = F, fire_name = fire_name
  )

  if(write_file) {
    nn <- paste(fire_name, "-wave2.rds", sep = "")
    saveRDS(wave2, file.path(target_dir, nn))
  }
  rr <- round(max(wave2$overlap), 4)
  print(paste("wave 2; overlap max: ", rr, sep = ""))

  # wave3: 10000 iter reproducing the 1000 best, to reach the max
  print(paste("Phase: search maximum (1000 best). Fire: ", fire_name, sep = ""))
  wave3 <- explore_likelihood_iterate(
    n_waves = 20, data = wave2, n = 500, n_best = 1000,
    var_factor = c(1.5, 2, 1),
    support = sup, spread_data = spread_data,
    centre = "global", sobol = TRUE, phase = "search_max_1000",
    write = F, fire_name = fire_name
  )

  if(write_file) {
    nn <- paste(fire_name, "-wave3.rds", sep = "")
    saveRDS(wave3, file.path(target_dir, nn))
  }
  rr <- round(max(wave3$overlap), 4)
  print(paste("wave 3; overlap max: ", rr, sep = ""))

  # wave4: 10000 iter reproducing 100 best
  print(paste("Phase: search maximum (100 best). Fire: ", fire_name, sep = ""))
  wave4 <- explore_likelihood_iterate(
    n_waves = 40, data = wave3, n = 500, n_best = 100,
    var_factor = c(1.5, 2, 1),
    support = sup, spread_data = spread_data,
    centre = "global", sobol = TRUE, phase = "search_max_100",
    write = F, fire_name = fire_name
  )

  if(write_file) {
    nn <- paste(fire_name, "-wave4.rds", sep = "")
    saveRDS(wave4, file.path(target_dir, nn))
  }
  rr <- round(max(wave4$overlap), 4)
  print(paste("wave 4; overlap max: ", rr, sep = ""))

  # merge and export
  explore <- rbind(wave1, wave2, wave3, wave4)
  # this was a mistake, because explore_.._iterate is cumulative!!!
  # wave4 already had all the waves.

  nn <- paste(fire_name, "-simulations_list.rds", sep = "")
  saveRDS(explore, file.path(target_dir, nn))

  # waveplot including test data
  nn <- paste(fire_name, "-waveplot_overlap.png", sep = "")
  png(filename = file.path(target_dir, nn),
      width = 22, height = 11, units = "cm", res = 300)
  wave_plot(explore, "overlap", alpha = 0.01,
            rc = c(2, 3), tit = fire_name)
  dev.off()

  # remove temporal files
  all_files_saved <- list.files(target_dir,
                                pattern = "-wave1|-wave2|-wave3|-wave4")
  if(length(all_files_saved) > 0) {
    lapply(all_files_saved, function(f) unlink(file.path(target_dir, f)))
  }
  # clean
  remove(full_data, spread_data, wave1, wave2, wave3, wave4, explore)
  gc()
}
# started at 2024-08-28; 00:58 h

# # Amend a mistake. Waves were cumulative, but I merged them.
# nw <- c(10000, 10000, 10000, 20000)
# nc <- cumsum(nw)
# ncum <- cumsum(nc)
# from <- ncum[3]+1
# to <- max(ncum)
#
# files <- list.files(target_dir, "-simulations_list.rds")
# for(f in 1:length(files)) {
#   print(f)
#   simslong <- readRDS(file.path(target_dir, files[f]))
#   to_local <- nrow(simslong)
#   print(to_local)
#   sims <- simslong[from:to_local, ]
#   sims$phase[sims$phase == "above_threshold"] <- "search_max_100"
#   sims$phase[sims$phase == "search_maximum"] <- "search_max_1000"
#
#   firename <- strsplit(files[f], "-simulations_list.rds")[[1]]
#   saveRDS(
#     sims,
#     file.path(target_dir, paste(firename, "-simulations.rds", sep = ""))
#   )
# }

# Sample fire-wise posteriors ------------------------------------------------

# For every fire, define the kernel based on the highest overlap values, to put
# the quadratic part of the expquad kernel where there are many points.

# In a few fires with low overlap, the overlap is very noisy, showing very large
# values that occurred with very low frequency. If those are just noise,
# setting the kernel based on those particles will bring problems. Perhaps the
# kernel width might be forcet to include the 1000 best particles.

# Constants
n_best <- 500 # number of best particles; the min overlap will be the threshold
sigma_kernel <- 0.025 # below the threshold, kernel for abc_prob
sd_factor_prior <- 5  # widen prior_sd will be sample_sd * sd_factor_prior
neff_get <- 10000    # target neff
nwave <- 20000        # minimum size of simulation wave
selm_form <- cbind(intercept, vfi, tfi, slope, wind, steps) ~ 1

target_dir_samples <- file.path("files", "posterior_samples_stage1")

# table to store info
sims_data <- size_data
sims_data$efficiency <- NA
sims_data$neff <- NA
sims_data$nsim <- NA
sims_data$overlap_high <- NA
sims_data$overlap_mean <- NA

for(f in 1:n_fires) {

  fire_id <- size_data$fire_id[f]
  mm <- paste(f, ": Fire ", fire_id,
              " --------------------------------------------------",
              sep = "")
  message(mm)

  # Import landscape
  landscape_file <- size_data$file[f]
  full_data <- readRDS(file.path(data_dir, landscape_file))
  spread_data <- full_data[c("landscape", "ig_rowcol",
                             "burned_layer", "burned_ids",
                             "counts_veg")]

  # Set multicore
  if(fire_id %in% c("2015_50", "1999_2140469994_r")) {
    registerDoMC(10)
  } else {
    registerDoMC(15)
  }

  # Load and prepare previous simulations
  simulations_file <- paste(fire_id, "-simulations.rds", sep = "")
  sims <- readRDS(file.path(target_dir, simulations_file))
  sup <- rbind(params_lower, params_upper)
  sup["params_upper", "steps"] <- size_data$steps_upper[size_data$fire_id == fire_id]
  # unconstrain parameters
  sims$par_values_unc <- logit_scaled(sims$par_values, sup)
  finite <- apply(sims$par_values_unc, 1, function(x) all(is.finite(x))) |>  t()
  sims <- sims[finite, ]

  # Define acceptance kernel and compute abc-prob
  perc_best <- 1 - n_best / nrow(sims)
  ov_high <- quantile(sims$overlap, probs = perc_best, method = 8) |> unname()
  sims$prob <- expquad_flat(sims$overlap, ov_high, sigma = sigma_kernel)

  # Resample particles to fit proposal distribution
  nfit <- 10000
  sims_sample <- sims[sims$prob > 0.01, ]
  ids <- sample(1:nrow(sims_sample), nfit, replace = T,
                prob = sims_sample$prob)
  sims_fit_prop <- sims_sample[ids, ]
  prop_df <- as.data.frame(sims_fit_prop$par_values_unc)
  colnames(prop_df) <- par_names

  # Fit proposal density (multivariate skew-normal, or simply mvn if it fails)
  dp_prop <- selm_robust(prop_df)
  sn_prop_dist <- makeSECdistr(dp_prop, "SN", "proposal", par_names) #just to plot

  # Define prior
  dp_prior <- list(
    "xi" = colMeans(prop_df),
    "Omega" = cov(prop_df) * sd_factor_prior ^ 2,
    "alpha" = rep(0, n_coef)
  )
  sn_prior_dist <- makeSECdistr(dp_prior, "SN", "prior", par_names)

  # Sample
  wave <- 1
  out <- NULL
  neff <- 0
  sim_count <- 0

  while(neff < neff_get1) {

    # Search particles to simulate
    got <- 0
    ntry <- nwave * 5
    candidates <- NULL

    while(got < nwave) {
      candidates_try <- rmsn(ntry, dp = dp_prop)
      log_prior <- dmsn(candidates_try, dp = dp_prior, log = T)
      log_prop <- dmsn(candidates_try, dp = dp_prop, log = T)
      prob_sim <- exp(log_prior - log_prop)
      ids_sim <- runif(ntry) < prob_sim
      new <- sum(ids_sim)
      if(new > 0) {
        candidates <- rbind(candidates, candidates_try[ids_sim, ])
        got <- got + new

        # optional message:
        mm <- paste("Wave ", wave, "; got ", got, sep = "")
        message(m)
       }
    }

    # Simulate particles
    candidates_const <- invlogit_scaled(candidates, sup)
    colnames(candidates_const) <- par_names
    out_tmp <- similarity_simulate_parallel(candidates_const,
                                             fire_data = spread_data)
    out_tmp$prob <- expquad_flat(out_tmp$overlap, ov_high, sigma_kernel)
    out <- rbind(out, out_tmp)

    # process results to plot advance
    w <- out$prob / sum(out$prob)
    neff <- 1 / sum(w ^ 2)
    sim_count <- nrow(out)
    rate <- neff / sim_count * 100
    ov_mean <- sum(w * out$overlap)

    # Print progress
    mm <- paste("Wave ", wave, " (", fire_id, ") *****************", sep = "")
    message(mm)
    mm <- paste("Efficiency =", round(rate, 2), "%")
    message(mm)
    mm <- paste("Neff =", round(neff, 2))
    message(mm)
    mm <- paste("Simulations =", sim_count)
    message(mm)
    mm <- paste("Overlap high =", round(ov_high, 3))
    message(mm)
    mm <- paste("Overlap mean =", round(ov_mean, 3))
    message(mm)

    # save temporal result if large
    if(f > 40) {
      nn <- paste(fire_id, "-samples_temporal_wave_", wave, ".rds", sep = "")
      saveRDS(out_tmp, file.path(target_dir_samples, nn))
    }
    wave <- wave + 1
  }

  result <- list(
    samples = out,
    dp_prop = dp_prop,
    dp_prior = dp_prior,
    neff = neff,
    sim_count = sim_count,
    fire_id = fire_id,
    overlap_high = ov_high,
    overlap_mean = ov_mean
  )

  nn <- paste(fire_id, "-samples.rds", sep = "")
  saveRDS(result, file.path(target_dir_samples, nn))

  # remove temporal files
  patt <- paste(fire_id, "-samples_temporal_wave_", sep = "")
  all_files_saved <- list.files(target_dir_samples,
                                pattern = patt)
  if(length(all_files_saved) > 0) {
    lapply(all_files_saved, function(f) unlink(file.path(target_dir_samples, f)))
  }
  gc()

  # save data
  sims_data$efficiency[f] <- rate
  sims_data$neff[f] <- neff
  sims_data$nsim[f] <- sim_count
  sims_data$overlap_high[f] <- ov_high
  sims_data$overlap_mean[f] <- ov_mean
}

write.csv(sims_data, file.path(target_dir_samples, "simulations_data.csv"))

# start at 30/08/2024, 02:35 h
# ended at 30/08/2024, 10:30 h


# Merge samples in a single array -----------------------------------------

# To run mcmc on the hyperparameters of the hierarchical model, resample all
# posteriors to 100000 samples, and save array.
# Also, export list with priors.

target_dir <- file.path("files", "posterior_samples_stage1")
samples_files <- list.files(target_dir, pattern = "-samples.rds")
J1 <- length(samples_files)

# import a list with matrix-samples
samples_list <- lapply(samples_files, function(s) {
  rrr <- readRDS(file.path(target_dir, s))
  return(rrr)
})

fire_ids <- lapply(samples_list, function(x) x$fire_id) |> unlist()
names(samples_list) <- fire_ids

# tidy data with steps bounds
rownames(size_data) <- size_data$fire_id
size_data <- size_data[fire_ids, ]

sup <- rbind(params_lower, params_upper)

# resample to 100000 particles
samples_list_boot <- vector("list", J1)
names(samples_list_boot) <- fire_ids
for(j in 1:J1) {
  # j = 1
  sup_local <- sup
  sup_local["params_upper", "steps"] <- size_data$steps_upper[j]
  xx <- samples_list[[j]]$samples

  pp <- xx$par_values
  # remove inf, checked at logit scale
  pp_logit <- logit_scaled(pp, sup_local)
  keep <- apply(pp_logit, 1, function(x) all(is.finite(x)))

  # turn into logit-log scale
  pp_unc <- unconstrain(pp[keep, ], sup_local)

  # resample
  ids <- sample(1:nrow(pp_unc), size = 100000, prob = xx$prob[keep], replace = T)
  samples_list_boot[[j]] <- pp_unc[ids, ]
}

samples_arr <- abind::abind(samples_list_boot, along = 3)
dimnames(samples_arr) <- list(
  iter = 1:100000,
  param = par_names,
  fire_id = fire_ids
)

saveRDS(samples_arr, file.path(target_dir, "samples_boot_all_fires.rds"))

samples_list[[1]]$dp_prior

mu0 <- sapply(1:J1, function(i) samples_list[[i]]$dp_prior$xi)
colnames(mu0) <- fire_ids
S0 <- sapply(1:J1, function(i) {samples_list[[i]]$dp_prior$Omega},
             simplify = "array")
dimnames(S0)[[3]] <- fire_ids

str(mu0); str(S0)
priors <- list(mu0 = mu0, S0 = S0)
saveRDS(priors, file.path(target_dir, "priors_stage1.rds"))

# BELOW: trials that did not work ----------------------------------------

# Multivariate-skew-normal approximation to the posterior ----------------

target_dir_samples <- file.path("files", "posterior_samples_stage1")
selm_form <- cbind(intercept, vfi, tfi, slope, wind, steps) ~ 1

for(f in 1:n_fires) {
  # f = 56
  fire_id <- size_data$fire_id[f]
  mm <- paste(f, ": Fire ", fire_id,
              " --------------------------------------------------",
              sep = "")
  message(mm)

  # Import landscape
  landscape_file <- size_data$file[f]
  full_data <- readRDS(file.path(data_dir, landscape_file))
  spread_data <- full_data[c("landscape", "ig_rowcol",
                             "burned_layer", "burned_ids",
                             "counts_veg")]
  sup <- rbind(params_lower, params_upper)
  sup["params_upper", "steps"] <- size_data$steps_upper[size_data$fire_id == fire_id]

  # Set multicore
  if(fire_id %in% c("2015_50", "1999_2140469994_r")) {
    registerDoMC(10)
  } else {
    registerDoMC(15)
  }

  # Load samples
  posterior_file <- paste(fire_id, "-samples.rds", sep = "")
  posterior_list <- readRDS(file.path(target_dir_samples, posterior_file))
  samples <- posterior_list$samples
  # unconstrain parameters
  samples$par_values_unc <- logit_scaled(samples$par_values, sup)
  finite <- apply(samples$par_values_unc, 1, function(x) all(is.finite(x)))
  samples <- samples[finite, ]

  # resample to get N = 10000
  set.seed(453459)
  ids <- sample(1:nrow(samples), size = 10000, prob = samples$prob, replace = T)
  samples_boot <- samples[ids, ]

  # Fit skew-normal approximation
  samples_boot_df <- as.data.frame(samples_boot$par_values_unc)
  colnames(samples_boot_df) <- par_names
  # Fit proposal density (multivariate skew-normal, or simply mvn if it fails)
  dp_posterior <- selm_robust(samples_boot_df)
  dist_posterior <- makeSECdistr(dp_posterior, "SN", "posterior", par_names)

  # Simulate overlap from the approximate posterior
  params_test <- rmsn(1000, dp = dp_posterior)
  params_test_const <- invlogit_scaled(params_test, sup)
  sims_test <- similarity_simulate_parallel(params_test_const,
                                            fire_data = spread_data)

  # Plot density to assess how good is the approximation
  dobs <- density(samples_boot$overlap, from = 0, to = 1, n = 2^11)
  dsim <- density(sims_test$overlap, from = 0, to = 1, n = 2^11)
  ym <- max(c(dobs$y, dsim$y))

  nn <- paste(fire_id, "-overlap_check_msn.png", sep = "")
  png(filename = file.path(target_dir_samples, nn),
      width = 10, height = 9, units = "cm", res = 300)
  plot(dobs, ylim = c(0, ym), xlab = "overlap", main = fire_id)
  lines(dsim, col = "red")
  dev.off()

  # Save result with approximated msn
  posterior_list$dp_posterior <- dp_posterior
  nn <- paste(fire_id, "-samples_msn_approx.rds", sep = "")
  saveRDS(posterior_list, file.path(target_dir_samples, nn))
}


# Multivariate normal mixture approximation to the posterior ----------------

library(MGMM)
target_dir_samples <- file.path("files", "posterior_samples_stage1")
selm_form <- cbind(intercept, vfi, tfi, slope, wind, steps) ~ 1

for(f in 1:n_fires) {
  # f = 25
  fire_id <- size_data$fire_id[f]
  mm <- paste(f, ": Fire ", fire_id,
              " --------------------------------------------------",
              sep = "")
  message(mm)

  # Import landscape
  landscape_file <- size_data$file[f]
  full_data <- readRDS(file.path(data_dir, landscape_file))
  spread_data <- full_data[c("landscape", "ig_rowcol",
                             "burned_layer", "burned_ids",
                             "counts_veg")]
  sup <- rbind(params_lower, params_upper)
  sup["params_upper", "steps"] <- size_data$steps_upper[size_data$fire_id == fire_id]

  # Set multicore
  if(fire_id %in% c("2015_50", "1999_2140469994_r")) {
    registerDoMC(10)
  } else {
    registerDoMC(15)
  }

  # Load samples
  posterior_file <- paste(fire_id, "-samples.rds", sep = "")
  posterior_list <- readRDS(file.path(target_dir_samples, posterior_file))
  samples <- posterior_list$samples
  # unconstrain parameters
  samples$par_values_unc <- logit_scaled(samples$par_values, sup)
  finite <- apply(samples$par_values_unc, 1, function(x) all(is.finite(x)))
  samples <- samples[finite, ]

  # resample to get N = 10000
  set.seed(453459)
  ids <- sample(1:nrow(samples), size = 10000, prob = samples$prob, replace = T)
  samples_boot <- samples[ids, ]

  # Fit mvnm approximation, choosing number of components
  samples_boot_df <- as.data.frame(samples_boot$par_values_unc)
  colnames(samples_boot_df) <- par_names
  # Fit posterior density (multivariate skew-normal, or simply mvn if it fails)

  # ?MGMM::ChooseK()
  # ?MGMM::FitGMM
  # getk <- ChooseK(
  #   samples_boot$par_values_unc,
  #   k0 = 2,
  #   k1 = 6,
  #   boot = 0#,
  #   # init_means = NULL,
  #   # fix_means = FALSE,
  #   # init_covs = NULL,
  #   # init_props = NULL,
  #   # maxit = 10,
  #   # eps = 1e-04,
  #   # report = TRUE
  # )

  kk <- 3
  ff <- FitGMM(
    samples_boot$par_values_unc,
    k = kk, maxit = 300
  )

  # Simulate overlap from the approximate posterior
  params_test <- rGMM(n = 1e3, d = 6, k = kk,
                      means = ff@Means, covs = ff@Covariances,
                      pi = ff@Proportions)
  params_test_const <- invlogit_scaled(params_test, sup)
  sims_test <- similarity_simulate_parallel(params_test_const,
                                            fire_data = spread_data)

  # Plot density to assess how good is the approximation
  dobs <- density(samples_boot$overlap, from = 0, to = 1, n = 2^11)
  dsim <- density(sims_test$overlap, from = 0, to = 1, n = 2^11)
  ym <- max(c(dobs$y, dsim$y))

  # nn <- paste(fire_id, "-overlap_check_msn.png", sep = "")
  # png(filename = file.path(target_dir_samples, nn),
  #     width = 10, height = 9, units = "cm", res = 300)
  plot(dobs, ylim = c(0, ym), xlab = "overlap", main = fire_id)
  lines(dsim, col = "red")
  # dev.off()

  # # Save result with approximated msn
  # posterior_list$dp_posterior <- dp_posterior
  # nn <- paste(fire_id, "-samples_msn_approx.rds", sep = "")
  # saveRDS(posterior_list, file.path(target_dir_samples, nn))
}


## It didn't work using normal mixtures


# MCMC --------------------------------------------------------------------

library(adaptMCMC)

posterior_ll <- function(x, spread_data, ov_high, sigma, support) {
  xx <- invlogit_scaled_vec(x, support)
  ov <- similarity_simulate_particle(xx, fire_data = spread_data)
  loglike <- expquad_flat(ov, ov_high, sigma = sigma, log = T)
  logprior <- sum(dlogis(xx, log = T))
  return(loglike + logprior)
}

target_dir_samples <- file.path("files", "posterior_samples_stage1")

f = 1
fire_id <- size_data$fire_id[f]
mm <- paste(f, ": Fire ", fire_id,
            " --------------------------------------------------",
            sep = "")
message(mm)

# Import landscape
landscape_file <- size_data$file[f]
full_data <- readRDS(file.path(data_dir, landscape_file))
spread_data <- full_data[c("landscape", "ig_rowcol",
                           "burned_layer", "burned_ids",
                           "counts_veg")]
sup <- rbind(params_lower, params_upper)
sup["params_upper", "steps"] <- size_data$steps_upper[size_data$fire_id == fire_id]

# Set multicore
if(fire_id %in% c("2015_50", "1999_2140469994_r")) {
  registerDoMC(10)
} else {
  registerDoMC(15)
}

# Load samples
posterior_file <- paste(fire_id, "-samples.rds", sep = "")
posterior_list <- readRDS(file.path(target_dir_samples, posterior_file))
samples <- posterior_list$samples

# unconstrain parameters
samples$par_values_unc <- logit_scaled(samples$par_values, sup)
finite <- apply(samples$par_values_unc, 1, function(x) all(is.finite(x)))
samples <- samples[finite, ]

# resample to get N = 10000
set.seed(453459)
ids <- sample(1:nrow(samples), size = 10000, prob = samples$prob, replace = T)
samples_boot <- samples[ids, ]

V <- cov(samples_boot$par_values_unc)
start <- samples_boot$par_values_unc[which.max(samples_boot$overlap), ]

m0 <- MCMC(p = posterior_ll, n = 50000, adapt = 20000,
           scale = V, init = start,
           acc.rate = 0.243,
           spread_data = spread_data,
           ov_high = 1, posterior_list$overlap_high,
           sigma = 0.1, # 0.025 before
           support = sup)

ss <- m0$samples[-(1:20000), ]
for(i in 1:6) {acf(ss[, i], lag.max = 2000)}
plot(ss[, 5], type = "l")
# Muuuuuuy ineficiente
m0$log.p
m0$cov.jump # No salta nada!!
m0$cov.jump |> diag() |> sqrt() # pequeñísimo el kernel de movida.
V |> diag() |> sqrt()

pairs(samples_boot$par_values_unc, col = rgb(0, 0, 0, 0.01), pch = 19)



# Assessing the prior effect ----------------------------------------------

lim <- 5
a <- 10
target_samples <- rsn(1000, alpha = a)
prop_mu <- mean(target_samples)
prop_sd <- sd(target_samples)



prior_factor <- 5
prior_mu <- prop_mu
prior_sd <- prop_sd * prior_factor


# visualize prior/proposal
curve(dsn(x, alpha = a), from = -lim, to = lim, n = 300)
curve(dnorm(x, prop_mu, prop_sd), add = T, col = "green")
curve(dnorm(x, prior_mu, prior_sd), add = T, col = "blue", lty = 2)
curve(dnorm(x, prior_mu, prior_sd) / dnorm(x, prop_mu, prop_sd),
      add = T, col = "red")

# Explore varying sd_factor:

lim <- 5
a <- 0.5

par(mfrow = c(3, 3))
for(prior_factor in 2:10) {
  target_samples <- rsn(1000, alpha = a)
  prop_mu <- mean(target_samples)
  prop_sd <- sd(target_samples)
  # prior_factor <- 5
  prior_mu <- prop_mu
  prior_sd <- prop_sd * prior_factor

  xseq <- seq(-lim, lim, length.out = 2000)
  xd <- diff(xseq)[1]
  dpost_un <- dnorm(xseq, prior_mu, prior_sd) * dsn(xseq, alpha = a)
  z <- sum(dpost_un * xd)
  dpost <- dpost_un / z
  dprior <- dnorm(xseq, prior_mu, prior_sd)

  plot(dpost ~ xseq, type = "l", col = "red", ylim = c(0, max(dpost) *1.2),
       main = paste("prior factor =", prior_factor))
  curve(dsn(x, alpha = a), lwd = 1, add = T, n = 2000)
  curve(dnorm(x, prior_mu, prior_sd), col = "blue", add = T, n = 2000)
}
par(mfrow = c(1, 1))



####
# uniform-mixture target

lim <- 8
target_samples <- c(runif(700, -3, -2), runif(300, 1, 2))
prop_mu <- mean(target_samples)
prop_sd <- sd(target_samples)

prior_factor <- 5
prior_mu <- prop_mu
prior_sd <- prop_sd * prior_factor


# visualize prior/proposal
curve(dunif(x, -3, -2) * 0.7 + dunif(x, 1, 2) * 0.3,
      from = -lim, to = lim, n = 300, ylab = NA)
curve(dnorm(x, prop_mu, prop_sd), add = T, col = "green")
curve(dnorm(x, prior_mu, prior_sd), add = T, col = "blue", lty = 2)
curve(dnorm(x, prior_mu, prior_sd) / dnorm(x, prop_mu, prop_sd),
      add = T, col = "red")

# Explore varying sd_factor:

lim <- 10

par(mfrow = c(3, 3))
for(prior_factor in 2:10) {
  target_samples <- c(runif(700, -3, -2), runif(300, 1, 2))
  prop_mu <- mean(target_samples)
  prop_sd <- sd(target_samples)
  # prior_factor <- 5
  prior_mu <- prop_mu
  prior_sd <- prop_sd * prior_factor

  xseq <- seq(-lim, lim, length.out = 2000)
  xd <- diff(xseq)[1]
  dpost_un <-
    dnorm(xseq, prior_mu, prior_sd) *
    (dunif(xseq, -3, -2) * 0.7 + dunif(xseq, 1, 2) * 0.3)
  z <- sum(dpost_un * xd)
  dpost <- dpost_un / z
  dprior <- dnorm(xseq, prior_mu, prior_sd)

  plot(dpost ~ xseq, type = "l", col = "red", ylim = c(0, max(dpost) *1.2),
       main = paste("prior factor =", prior_factor))
  curve(dunif(x, -3, -2) * 0.7 + dunif(x, 1, 2) * 0.3,
        col = "black", add = T, n = 2000)
  curve(dnorm(x, prior_mu, prior_sd), col = "blue", add = T, n = 2000)
}
par(mfrow = c(1, 1))