# This code inherits from
# <overlap_function_estimation-many_intercepts.R>.

# Here I run the simulation of small fires, to see how GAMs work.
# Evaluate 40000 particles, fit GAM, test gam on 2000 particles.

# I do the same as before, but with a model based on flammability indices
# (FI), which has 6 data-level parameters:

# intercept, vfi, tfi, slope, wind, steps

# I will try metrics besides the overlap

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
library(mvtnorm)
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

data_dir <- file.path("data", "focal fires data", "landscapes_FI")
filenames <- list.files(data_dir)
n_fires <- length(filenames)

# dir to save output
target_dir <- file.path("files", "overlaps_FI")

# load file with constants to standardize
fi_params <- readRDS(file.path("data", "NDVI_regional_data",
                               "flammability_indices.rds"))
slope_sd <- fi_params$slope_term_sd

# constants for fire spread simulation
upper_limit <- 1
n_veg <- 5
veg_names <- c("wet", "subalpine", "dry", "shrubland", "grassland")
n_terrain <- 2
terrain_names <- c("slope", "wind")
n_fi <- 2        # flammability indices
n_nd <- n_fi + 1 # non-directional terms in the linear predictor
par_names <- c("intercept", "vfi", "tfi", terrain_names, "steps")
n_coef <- length(par_names)

# number of fires to simulate by particle
n_rep <- 20

ext_alpha <- 50
ext_beta <- 30

params_lower <- c(-ext_alpha, rep(0, n_coef-2), 5)
params_upper <- c(ext_alpha, rep(ext_beta, n_coef-2), NA)
names(params_lower) <- names(params_upper) <- par_names
params_upper["slope"] <- ext_beta / slope_sd

gam_formula <- formula(
  y ~

    # marginal effects
    s(intercept, bs = basis, k = k_side, fx = fixed) +
    s(vfi, bs = basis, k = k_side, fx = fixed) +
    s(tfi, bs = basis, k = k_side, fx = fixed) +
    s(slope, bs = basis, k = k_side, fx = fixed) +
    s(wind, bs = basis, k = k_side, fx = fixed) +
    s(steps, bs = basis, k = k_side, fx = fixed) +

    # interactions
    ti(intercept, vfi, k = k_int, bs = basis, fx = fixed) +
    ti(intercept, tfi, k = k_int, bs = basis, fx = fixed) +
    ti(intercept, slope, k = k_int, bs = basis, fx = fixed) +
    ti(intercept, wind, k = k_int, bs = basis, fx = fixed) +
    ti(intercept, steps, k = k_int, bs = basis, fx = fixed) +

    ti(vfi, tfi, k = k_int, bs = basis, fx = fixed) +
    ti(vfi, slope, k = k_int, bs = basis, fx = fixed) +
    ti(vfi, wind, k = k_int, bs = basis, fx = fixed) +
    ti(vfi, steps, k = k_int, bs = basis, fx = fixed) +

    ti(tfi, slope, k = k_int, bs = basis, fx = fixed) +
    ti(tfi, wind, k = k_int, bs = basis, fx = fixed) +
    ti(tfi, steps, k = k_int, bs = basis, fx = fixed) +

    ti(slope, wind, k = k_int, bs = basis, fx = fixed) +
    ti(slope, steps, k = k_int, bs = basis, fx = fixed) +

    ti(wind, steps, k = k_int, bs = basis, fx = fixed)
)

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

# Function to compute size-based metrics. It takes the vectors of burned pixels
# by vegetation type, and returns a vector of the same length with the relative
# absolute differences.
compare_size_veg <- function(counts_obs, counts_sim) {
  size_real <- sum(counts_obs)
  dr <- abs(counts_obs - counts_sim) / size_real
  return(dr)
}

# compute scale for a gaussian kernel where y = p * pmax at xcrit
get_scale <- function(xcrit, p) {
  ss <- xcrit / sqrt(-log(p))
  return(ss)
}

xcrit <- 0.5; p <- 0.5
kscale <- get_scale(xcrit, p)

# NOT USED
# compare_size <- function(fire1, fire2,
#                          qcrit = 1.5, p = 0.5,
#                          scale = get_scale(qcrit, p)) {
#   s1 <- sum(fire1$counts_veg); s2 <- sum(fire2$counts_veg)
#   sq <- 1 - max(c(s1, s2)) / min(c(s1, s2))
#   y <- exp(-(sq / scale) ^ 2)
#   return(y)
# }

# s1 <- c(100, 200, 300, 400, 500)
# s2 <- c(10, 210, 900, 200, 1000)
# compare_size_veg(s1, s2) |> sum()
# compare_size_veg(s1*2, s2*2) |> sum()
# compare_size_veg(s1, s2) |> prod()
# compare_size_veg(s1*2, s2*2) |> prod()

# Compute acceptance probability with an exponential quadratic function,
# based on overlap (similarity) and size difference by veg type.
kernel_expquad <- function(overlap, sizediff, scale = kscale) {
  # ## test:
  # overlap <- 0.5
  # sizediff <- compare_size_veg(s1, s2)
  # scale = kscale
  # ##
  lp <-
    ((1-overlap) / scale) ^ 2 +
    sum((sizediff / scale) ^ 2)
  return(exp(-lp))
}

# Similar as above, but scaling the overlap to max = 1, so
# a similar kernel is used in all fires.
kernel_expquad2 <- function(metrics, overlap_max,
                            ov_scale = 0.05, diff_scale = 0.1, log = F) {

  if(is.null(dim(metrics))) {
    metrics <- matrix(metrics, nrow = 1)
  }
  # ## test:
  # overlap <- runif(1)
  # sizediff <- matrix(runif(1 * 5), 1, 5)
  # ov_scale = 0.05; diff_scale = 0.1; overlap_max = 1; log = F
  # ##

  # Scale overlap to have maximum = 1
  ov_unit <- metrics[1] / overlap_max
  ov_unit[ov_unit > 1] <- 1
  ov_dist <- 1 - ov_unit
  metrics[1] <- ov_dist

  Sigma <- diag(c(ov_scale, rep(diff_scale, ncol(metrics)-1)) ^ 2)
  out <- dmvnorm(metrics, rep(0, ncol(metrics)), Sigma,
                 log = log, checkSymmetry = F)

  return(out)
}

kernel_expquad3 <- function(overlap, overlap_max,
                            ov_scale = 0.05, log = F) {

  # Scale overlap to have maximum = 1
  ov_unit <- overlap / overlap_max
  ov_unit[ov_unit > 1] <- 1
  ov_dist <- 1 - ov_unit

  out <- dnorm(ov_dist, sd = ov_scale, log = log)

  return(out)
}


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
  sdiff <- matrix(NA, n_sim, n_veg)

  # terrain coefficients
  ct <- particle[(n_nd+1):(n_coef-1)]

  # non-directional effects (intercepf + fis)
  coef_nd <- particle[1:n_nd]
  nd_layer <- make_nd_layer(fire_data$landscape[, , c("vfi", "tfi")], coef_nd)

  # simulate and compute metrics
  for(i in 1:n_sim) { # simulated fires by particle
    fire_sim <- simulate_fire_compare_veg(
      nd_terms = nd_layer,
      vegetation = fire_data$landscape[, , 1],
      terrain = fire_data$landscape[, , c("elevation", "wdir", "wspeed")],
      coef_terrain = ct,
      ignition_cells = fire_data$ig_rowcol,
      upper_limit = upper_limit,
      steps = particle[n_coef],
      n_veg = n_veg
    )

    ov[i] <- overlap_spatial(
      fire_sim[c("burned_layer", "burned_ids")],
      fire_data[c("burned_layer", "burned_ids")]
    )

    sdiff[i, ] <- compare_size_veg(fire_data$counts_veg,
                                   fire_sim$counts_veg)

  }

  out <- cbind(ov, sdiff)
  colnames(out) <- c("ov", paste("sdiff", 1:n_veg, sep = ""))

  return(out)
}

# The same but returning many metrics. Currently used only for steps used.
# It returns a data.frame with averages.
similarity_simulate_particle_metrics <- function(particle, fire_data = NULL,
                                                 n_sim = n_rep) {

  metrics <- matrix(NA, n_sim, 3)
  colnames(metrics) <- c("overlap", "size_diff", "steps_used")

  # terrain coefficients
  ct <- particle[(n_nd+1):(n_coef-1)]

  # non-directional effects (intercepf + fis)
  coef_nd <- particle[1:n_nd]
  nd_layer <- make_nd_layer(fire_data$landscape[, , c("vfi", "tfi")], coef_nd)


  # simulate and compute metrics
  for(i in 1:n_sim) { # simulated fires by particle
    fire_sim <- simulate_fire_compare(
      nd_terms = nd_layer,
      vegetation = fire_data$landscape[, , 1],
      terrain = fire_data$landscape[, , c("elevation", "wdir", "wspeed")],
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

  # ### test:
  # particles_mat <- particles[waves_1 == w, ]
  # fire_data <- spread_data
  # n_sim <- 20
  # ###

  # turn particle matrix into list for parallel evaluation
  particles_list <- lapply(1:nrow(particles_mat),
                           function(x) particles_mat[x, ])

  # simulate in parallel
  result <- foreach(pp = particles_list) %dopar% {
    similarity_simulate_particle(pp, fire_data, n_sim)
  }

  # rbind list result
  metrics_arr <- abind::abind(result, along = 3)

  # for each particle, get the average abc_prob, computed before averaging
  # replicates
  abc_prob_mat <- apply(metrics_arr, c(1, 3), function(x) {
    kernel_expquad(x[1], x[2:(n_veg+1)], scale = kscale)
  })
  abc_prob_summ <- colMeans(abc_prob_mat)

  metrics_summ <- apply(metrics_arr, 3, colMeans) |> t()

  # return data.frame
  res <- data.frame(
    wave = NA,
    abc_prob = abc_prob_summ
  )
  res$metrics <- metrics_summ
  res$par_values <- particles_mat
  return(res)
}

# get_bounds: computes the metrics for the largest and smallest fires possible
get_bounds <- function(fire_data, n_sim = 1) {

  coef_burn_all <- c(1e6, rep(0, n_coef - 1))
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

wave_plot <- function(data, response = "abc_prob", alpha = 0.3, best = NULL,
                      x = "par_values", thres = NULL, bin = FALSE,
                      tit = NULL, rc = c(2, 3)) {

  if(!is.null(best)) {
    data <- data[order(data[, response], decreasing = TRUE), ]
  }

  if(bin) {
    data$bin <- jitter(as.numeric(data$abc_prob >= thres), factor = 0.2)
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

    if(!is.null(thres) & response == "abc_prob") {
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

  dd <- data$par_values[data$abc_prob >= thres, ]

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

  names_logit <- colnames(support)[colnames(support) != "steps"]
  for(j in names_logit) {
    xun[, j] <- qlogis((x[, j] - support[1, j]) / (support[2, j] - support[1, j]))
  }

  xun[, "steps"] <- log(x[, "steps", drop = F])

  return(xun)
}

# The same as above, but with all variables bounded by logit link
unconstrain2 <- function(x, support) {
  xun <- x
  for(j in 1:ncol(x)) {
    xun[, j] <- qlogis((x[, j] - support[1, j]) / (support[2, j] - support[1, j]))
  }
  return(xun)
}

unconstrain_vec <- function(x, support) {
  xun <- x

  llgg <- colnames(support) != "steps"
  xun[llgg] <- qlogis(
    (x[llgg] - support[1, llgg]) /
      (support[2, llgg] - support[1, llgg])
  )
  xun[!llgg] <- log(x[!llgg])

  return(xun)
}

unconstrain2_vec <- function(x, support) {
  xun <- x
  xun <- qlogis(
    (x - support[1, ]) /
    (support[2, ] - support[1, ])
  )
  return(xun)
}

# The inverse of unconstrain
constrain <- function(xun, support) {
  xc <- xun

  names_logit <- colnames(support)[colnames(support) != "steps"]
  for(j in names_logit) {
    xc[, j] <- plogis(xun[, j]) * (support[2, j] - support[1, j]) + support[1, j]
  }

  xc[, "steps"] <- exp(xun[, "steps", drop = F])

  return(xc)
}

# The inverse of unconstrain2
constrain2 <- function(xun, support) {
  xc <- xun
  for(j in 1:ncol(xun)) {
    xc[, j] <- plogis(xun[, j]) * (support[2, j] - support[1, j]) + support[1, j]
  }
  return(xc)
}

constrain_vec <- function(xun, support) {
  xc <- xun

  llgg <- colnames(support) != "steps" # logit positions
  xc[llgg] <- plogis(xun[llgg]) * (support[2, llgg] - support[1, llgg]) + support[1, llgg]

  xc[!llgg] <- exp(xun[!llgg])

  return(xc)
}

constrain2_vec <- function(xun, support) {
  xc <- xun
  xc <- plogis(xun) * (support[2, ] - support[1, ]) + support[1, ]
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
# accept_thres: abc_prob threshold below which particles will not be resampled.
#   overrides both n_best and p_best.
# prob_fn: function to transform abc_prob into probability. If not NULL,
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
    n_best <- sum(data$abc_prob >= accept_thres)
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
    data <- data[order(data$abc_prob, decreasing = TRUE), ]
  }

  # resample the best particles
  weight <- data$abc_prob

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
  m <- paste("Search starts at max abc_prob = ",
             round(max(data$abc_prob), 4), sep = "")
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

    m <- paste("Search wave ", last_wave + i, "; max abc_prob = ",
               round(max(res$abc_prob), 4), sep = "")

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

# function to evaluate ABC-posterior, based on an abc_prob threshold
fn_like_sim_bin <- function(x, support = NULL, fire_data = NULL, thres = NULL) {
  xx <- invlogit_scaled_vec(x, support)
  ov <- similarity_simulate_particle(xx, n_sim = n_rep, fire_data = fire_data) %>% mean
  if(ov >= thres) {
    return(sum(dlogis(x, log = TRUE)))
  } else {
    return(-Inf)
  }
}

# function to evaluate ABC-posterior, based on an abc_prob threshold
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
  use <- sim$abc_prob >= threshold
  colnames(sim$par_values) <- par_names

  samples_out <- sim$par_values[use, ]
  return(samples_out)
}


# Simulate fires ----------------------------------------------------------

# last_fire <- size_data[size_data$size_burn_rel < 0.1, ] %>% nrow() # 48
# which(size_data$fire_id == "1999_25j")

# (tested with fire 1, error at fire 26)

for(f in 26:n_fires) {

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
  print(paste("Phase: sobol. ", fire_name, sep = ""))
  for(w in 1:20) {
    # w = 1
    wlocal <- similarity_simulate_parallel(particles[waves_1 == w, ],
                                           spread_data)
    wlocal$wave <- w
    wlocal$phase <- "sobol"
    colnames(wlocal$par_values) <- par_names

    nn <- paste(fire_name, "-wave-", w, ".rds", sep = "")

    if(write_file) {
      saveRDS(wlocal, file.path(target_dir, nn))
    }

    wave1 <- rbind(wave1, wlocal)
    rr <- round(max(wave1$abc_prob), 10)
    print(paste("wave ", w, ", abc_prob max: ", rr, sep = ""))
  }

  # wave 2: reproduce 15 % best (500 * 20 = 10000)
  print(paste("Phase: search higher. Fire: ", fire_name, sep = ""))
  wave2 <- explore_likelihood_iterate(
    n_waves = 20, data = wave1, n = 500, p_best = 0.15,
    var_factor = c(1.5, 2, 1),
    support = sup, spread_data = spread_data,
    centre = "global", sobol = TRUE, phase = "search_higher",
    write = write_file, fire_name = fire_name
  )

  # 10000 iter reproducing the 100 best, to reach the max
  print(paste("Phase: search maximum. Fire: ", fire_name, sep = ""))
  wave3 <- explore_likelihood_iterate(
    n_waves = 100, data = wave2, n = 100, n_best = 100,
    var_factor = c(1.5, 2, 1),
    support = sup, spread_data = spread_data,
    centre = "global", sobol = TRUE, phase = "search_maximum",
    write = write_file, fire_name = fire_name
  )

  thres <- max(wave3$abc_prob) * 0.95

  # 10000 iter above thres
  print(paste("Phase: above threshold. Fire: ", fire_name, sep = ""))
  wave4 <- explore_likelihood_iterate(
    n_waves = 20, data = wave3, n = 500, accept_thres = thres,
    var_factor = c(1.5, 2, 1),
    support = sup, spread_data = spread_data,
    centre = "global", sobol = TRUE, phase = "above_threshold",
    write = write_file, fire_name = fire_name
  )

  # 5000 sobol particles in the rectangle of high prob.
  above <- wave4$abc_prob >= thres

  # reduce support
  sup_reduced <- reduce_support(wave4$par_values[above, ],
                                sup, prop = 1)

  # Throw particles in the reduced support.
  ss <- sobol(n = 5000, dim = n_coef)
  particles <- scale_params(ss, sup_reduced)
  waves_1 <- rep(1:10, each = 500)

  # loop to save
  wave1 <- NULL
  print(paste("Phase: sobol at high prob. Fire: ", fire_name, sep = ""))
  for(w in 1:10) {
    # w = 1
    wlocal <- similarity_simulate_parallel(particles[waves_1 == w, ],
                                           spread_data)
    wlocal$wave <- w + max(wave4$wave)
    wlocal$phase <- "sobol_2"
    colnames(wlocal$par_values) <- par_names

    nn <- paste(fire_name, "-wave-", w, ".rds", sep = "")

    if(write_file) {
      saveRDS(wlocal, file.path(target_dir, nn))
    }

    wave1 <- rbind(wave1, wlocal)
    rr <- round(max(wave1$abc_prob), 10)
    print(paste("wave ", w, ", abc_prob max: ", rr, sep = ""))
  }

  # Merge old and new data (including data used for testing)
  like_sim <- rbind(
    wave4,
    wave1[, colnames(wave4)]
  )

  nn <- paste(fire_name, "-waveplot.png", sep = "")
  png(filename = file.path(target_dir, nn),
      width = 20, height = 11, units = "cm", res = 300)
  wave_plot(like_sim, alpha = 0.05,
            thres = thres, tit = fire_name, rc = c(2, 3))
  dev.off()

  # Fit gam

  print(paste("Fitting GAM. Fire: ", fire_name, sep = ""))

  data_gam_bern <- as.data.frame(cbind(like_sim$par_values,
                                       y = as.numeric(like_sim$abc_prob >= thres)))

  k_side <- 15
  k_int <- 6
  basis <- "cr"
  fixed <- F

  gam_bern <- bam(
    gam_formula, data = data_gam_bern, family = "binomial",
    method = "fREML", discrete = T, nthreads = 8
  )

  # Sample GAM-approximated posterior

  print(paste("Sampling GAM. Fire: ", fire_name, sep = ""))

  sup_reduced <- reduce_support(like_sim$par_values[like_sim$abc_prob >= thres, ],
                                sup, 0.5)
  r_gam <- rejection_sample_parallel(700, gam_bern, sup_reduced,
                                     centre = 0.7,
                                     cores = 15)
  draws <- do.call("rbind", r_gam) %>% as.data.frame

  # Check GAM

  ids <- sample(1:nrow(draws), size = 1000, replace = F)
  ppmat <- as.matrix(draws[ids, ])

  print(paste("Checking GAM. Fire: ", fire_name, sep = ""))

  abc_prob_check <- similarity_simulate_parallel(particles = ppmat,
                                                 fire_data = spread_data)

  # acceptance probability from gam
  acp <- sum(abc_prob_check$abc_prob >= thres) / nrow(abc_prob_check) * 100

  nn <- paste(fire_name, "-gamcheck.png", sep = "")
  ttt <- paste(fire_name, "; acceptance prob = ", round(acp, 2), " %", sep = "")
  png(filename = file.path(target_dir, nn),
      width = 10, height = 10, units = "cm", res = 300)
  hist(abc_prob_check$abc_prob, xlim = c(0, 1), main = ttt,
       xlab = "abc_prob", ylab = "Relative frequency", freq = F)
  abline(v = thres, col = 2, lwd = 2)
  abline(v = mean(abc_prob_check$abc_prob), col = 4, lwd = 2)
  dev.off()

  # for(i in 1:6) {
  #   plot(like_sim$abc_prob ~ like_sim$metrics[, i],
  #        col = rgb(0, 0, 0, 0.03), pch = 19,
  #        main = colnames(like_sim$metrics)[i])
  # }

  # save result
  res <- list(
    fire_name = fire_name,
    veg_props = round(normalize(full_data$cells_by_veg) * 100, 2),
    like_sim = like_sim,             # first 40000 simulations
    thres = thres,                   # abc_prob threshold
    gam_bern = gam_bern,             # bernoulli GAM
    draws = draws,                   # draws sampled from the GAM
    abc_prob_check = abc_prob_check, # subsample from draws for which the abc_prob
    #   was computed
    acceptance = acp
  )

  nn <- paste(fire_name, "-simulations_list.rds", sep = "")
  saveRDS(res, file.path(target_dir, nn))

  # remove temporal files
  all_files_saved <- list.files(target_dir, pattern = "-wave-")
  if(length(all_files_saved) > 0) {
    lapply(all_files_saved, function(f) unlink(file.path(target_dir, f)))
  }
  # clean
  remove(full_data, spread_data, wave1, wave2, wave3, wave4, like_sim)
  gc()
}

## start run 2024-08-21, 21:32 h
## restart at 2024-08-22, 01:05 h

## Old comment:
# Start run of large fires at 18/07/2024 00:56
# In total it took ~ 4.3 days





# Exploring importance sampling -------------------------------------------

f <- 1
fire_file <- size_data$file[f]
fire_name <- size_data$fire_id[f]

sims_name <- file.path("files", "overlaps_FI",
                       paste(fire_name, "-simulations_list.rds", sep = ""))
sims <- readRDS(sims_name)
dat <- rbind(sims$like_sim,
             sims$abs_prob_check)

# Overlap will be optimized
full_data <- readRDS(file.path(data_dir, fire_file))
# subset data needed for spread (to be cloned across workers)
spread_data <- full_data[c("landscape", "ig_rowcol",
                           "burned_layer", "burned_ids",
                           "counts_veg")]

bb <- get_bounds(fire_data = spread_data)
sup <- rbind(params_lower, params_upper)
sup[2, "steps"] <- bb["largest", "steps_used"]

registerDoMC(7)

# Optimize overlap
fn <- function(xun) {
  x <- constrain2_vec(xun, support = sup)
  out <- similarity_simulate_particle(as.vector(x), fire_data = spread_data,
                                      n_sim = n_rep)
  out <- colMeans(out)
  ll <- out[1]
  return(log(ll))
}

opt <- optim(
  unconstrain2_vec(dat$par_values[which.max(dat$metrics[, "ov"]), , drop = F],
                   support = sup),
  fn,
  method = "BFGS",
  control = list(fnscale = -1, maxit = 100),
  hessian = TRUE
)
overlap_max <- opt$value |> exp()
ov_scale = 0.025

# Explore continuous kernel evaluations -----------------------------------

dat$ll <- kernel_expquad2(dat$metrics, overlap_max,
                          ov_scale = ov_scale,
                          diff_scale = 0.6)
dat$overlap <- dat$metrics[, "ov"]
dat$like_overlap <- kernel_expquad3(dat$overlap, overlap_max,
                                    ov_scale = ov_scale)

# wave_plot(dat, "ll", alpha = 0.05)
# wave_plot(dat, "overlap", alpha = 0.08)
wave_plot(dat, "like_overlap", alpha = 0.05)    ## usar esto

# El kernel multivariado toma una forma muy complicada


# Build MSN proposal ------------------------------------------------------

# select the half best particles, unconstrain, and fit sn model.
thresh2 <- max(dat$like_overlap) / 2
thresh3 <- 1
xtrain <- dat$par_values[dat$like_overlap >= thresh3, ]
dim(xtrain)
samples_train <- unconstrain2(xtrain, sup)
apply(samples_train, 2, range)

# filter inf
ids_ok <- sapply(1:nrow(samples_train), function(i) {
  return(all(is.finite(samples_train[i, ])) && !anyNA(samples_train[i, ]))
})
samples_train <- samples_train[ids_ok, ]

# Fit skew normal distrib proposal
selm_form <- cbind(intercept, vfi, tfi, slope, wind, steps) ~ 1

snfit <- sn::selm(formula = selm_form,
                  data = as.data.frame(samples_train), family = "SN")
dp <- coef(snfit, "DP", vector = FALSE)
var_factor <- 3
dp$Omega <- dp$Omega * var_factor
sndist <- makeSECdistr(dp, "SN", "proposal", par_names)

# Build sampler avoiding unuseful simulation --------------------------------

nget <- 10000
candidates <- NULL
got <- 0
nbatch <- 5e5

while(got < nget) {
  candidates_try <- rmsn(nbatch, dp = dp)
  log_prior <- rowSums(apply(candidates_try, 2, dlogis, log = T))
  log_prop <- dmsn(candidates_try, dp = dp, log = T)
  prob_sim <- exp(log_prior - log_prop)
  ids_sim <- runif(nbatch) < prob_sim
  new <- sum(ids_sim)
  if(new > 0) {
    candidates <- rbind(candidates, candidates_try[ids_sim, ])
    got <- got + new
    print(got)
  }
}

plot(sndist)
pairs(candidates[1:5000, ], col = rgb(0, 0, 0, 0.1), pch = 19, cex = 0.7)
# simulation prob is ~ 69 / 1e5, biutiful
# Siii, las distrib a simular es muy similar a la propuesta!


# Simulate particles with prob prior / proposal, and define weight as
# k(s) / k(0)

candidates_const <- constrain2(candidates, sup)
ov_evals <- similarity_simulate_parallel(candidates_const,
                                         fire_data = spread_data)
log_kernel <- kernel_expquad3(ov_evals$metrics[, "ov"], overlap_max,
                              ov_scale = 0.025,
                              log = T)
ov_evals$ll <- exp(log_kernel)
# wave_plot(ov_evals, "ll", alpha = 0.1)

### (todo esto es viejo, de importance sampling simple)
# mode_prop <- modeSECdistr(object = sndist)
# prop_max <- dmsn(mode_prop, dp = dp)
#
# # obtain the maximum prior / proposal density ratio
# densq <- function(x) {
#   exp(sum(dlogis(x, log = T)) - dmsn(x, dp = dp, log = T))
# }
# opt2 <- optim(mode_prop, densq, method = "BFGS",
#               control = list(fnscale = -1))
# M <- opt2$value * kernel_expquad3(1, 1, ov_scale = 0.025, log = F)
#
# accept_prob <- exp(log_kernel + log_prior - log(M) - log_prop)
# accept_prob <- exp(log_kernel + log_prior) / (M * exp(log_prop))
#
# hist(accept_prob)
# summary(accept_prob)
# # las probs son bajísimas, y si les hago importance,
# w <- accept_prob / sum(accept_prob)
# 1 / sum(w ^ 2) # el effective sample size es 1


# normalize weights
w <- ov_evals$ll / sum(ov_evals$ll)
1 / sum(w ^ 2)
# 287 / 10000 # Eficiencia del 2.87 %


# Repito pero simulando sólo una réplica por partícula (como debería ser!!)

ov_evals1 <- similarity_simulate_parallel(candidates_const,
                                         fire_data = spread_data,
                                         n_sim = 1)
log_kernel1 <- kernel_expquad3(ov_evals1$metrics[, "ov"], overlap_max,
                               ov_scale = 0.025,
                               log = T)
ov_evals1$ll <- exp(log_kernel1)

w1 <- ov_evals1$ll / sum(ov_evals1$ll)
1 / sum(w1 ^ 2)
# 305 / nrow(candidates) # Eficiencia del 3 %, pero con 1/simulaciones comparado
# con antes.

# Si la eff es del 2 % y queremos 5000 muestras por fuego,
# necesitamos simular ~ 5000 / 0.02 = 250000 partículas por fuego.
# Antes simulábamos 20 fuegos por partícula; si simulamos sólo 1,
# sería equivalente a simular 250000 / 20 = 12500 partículas con el
# esquema anterior (puede tardar 1 o 2 días).

# Esto puede ser peor si la propuesta no se asemeja tanto a la
# posterior, pero no es tan grave.
ids_re <- sample(1:nrow(candidates), size = 5000, prob = ov_evals1$ll,
                 replace = T)
ov_resamp <- ov_evals1[ids_re, ]
wave_plot(ov_resamp, "ll", alpha = 0.1)
pairs(ov_resamp$par_values[1:5000, ], col = rgb(0, 0, 0, 0.1), pch = 19, cex = 0.7)

# Efectividad 3 % implica guardar 1 / 30 muestras en un MCMC... no está
# taaaan mal. Pensá que el MCMC en R guardaba 1 / 1000
# Todo habría sido más barato si de entrada simulábamos un fuego por
# partícula.



# Y si uso un kernel absoluto, que ponga los overlaps bien puntudo --------


log_kernel1 <- kernel_expquad3(ov_evals1$metrics[, "ov"],
                               overlap_max = overlap_max,
                               ov_scale = 0.025,
                               log = T)
ov_evals1$ll <- exp(log_kernel1)

w1 <- ov_evals1$ll / sum(ov_evals1$ll)
1 / sum(w1 ^ 2) # lo mismo

# No no, mejor no. Usar un kernel relativo, o sea, variando el overlap_max
# entre fuegos, hace que la like sea cuadrática en el modo; si no, puede ponerse
# muy exponencial, con forma tipo laplaciana, y eso hace a la distri
# mucho más difícil de muestrear. O sea, el algoritmo se vuelve muy ineficiente.

# El mismo efecto tiene ponerse muy quisquilloso con el ov_scale:
# mientras más astringente, menos eficiente el simulador.



# Cuánto varía el overlap entre réplicas? ---------------------------------

nr <- 10
matsim <- matrix(NA, nrow(ov_resamp), nr)
for(k in 1:nr) {
  print(k)
  rr <- similarity_simulate_parallel(ov_resamp$par_values,
                                     fire_data = spread_data,
                                     n_sim = 1)
  matsim[, k] <- rr$metrics[, "ov"]
}

mat_all <- cbind(ov_resamp$metrics[, "ov"], matsim)
hist(mat_all)
mads <- apply(mat_all, 1, mad)
hist(mads)
hist(apply(mat_all, 1, sd))

# El overlap varía re poco entre réplicas!! Se podría usar un kernel aún
# más picante. Si la proposal density se estima de la misma forma,
# no debería ser un problema.
#
# Podemos agarrar las mejores partículas de la tanda anterior, reproducirlas
# mucho, y calcular este desvío, y en base a ese desvío definir la finitud
# del kernel.
# De la misma forma, usar sólo las partículas con mayor densidad
# para estimar la proposal skewnormal, tipo el percentil más alto.
#
#
# ESTO ES UNA BELLESA

# MCMC again?? ------------------------------------------------------------

# Si con 8 * 4 iter se logra algo decente, es una buena alternativa.
# ´

# Function to perform adaptive MCMC in parallel by hand. It's just a wrapper
# for adaptMCMC::MCMC, because MCMC.parallel does not allow to fix the inits.
MCMC_parallel <- function(fun, n, adapt, scale, init_list, acc_rate = 0.234,
                          n_chains, n_cores, ...) {

  registerDoMC(n_cores)

  runs <- foreach(pp = init_list) %dopar% {
    MCMC(p = fun, n = n, adapt = adapt, scale = scale, init = pp,
         acc.rate = acc_rate, ...)
  }

  return(runs)
}

# simulate fire
lp_fun <- function(x, support = NULL, fire_data = NULL,
                   overlap_max = 1, ov_scale = 0.025) {
  xx <- constrain2_vec(x, support)
  mets <- similarity_simulate_particle(xx, n_sim = 20, fire_data = fire_data)
  ov <- mean(mets[, 1])
  ll <- kernel_expquad3(ov, overlap_max, ov_scale = ov_scale, log = T)
  lp <- ll + sum(dlogis(x, log = TRUE))
  return(lp)
}

sampling_iters <- 1e5
adapt_iters <- 2000
sigma_init <- cov(samples_train)
nc <- 2
start <- lapply(1:nc, function(i) opt$par)


# # sample 5 posteriors
# r0 <- MCMC_parallel(lp_fun,
#                     n = sampling_iters + adapt_iters,
#                     adapt = adapt_iters,
#                     n_chains = nc, n_cores = nc,
#                     scale = sigma_init, init = start, acc_rate = 0.234,
#                     support = sup, fire_data = spread_data,
#                     ov_scale = 0.025, overlap_max = overlap_max)

r0 <- MCMC(
  lp_fun, init = opt$par,
  n = sampling_iters + adapt_iters,
  adapt = adapt_iters,
  scale = cov(samples_train),
  support = sup, fire_data = spread_data,
  ov_scale = 0.05, overlap_max = overlap_max, acc.rate = 0.234
)
xx <- as_draws_matrix(r0$samples[-(1:adapt_iters), ])
View(summarise_draws(xx))
# Nop, hay algunos neff muy muy bajos. (89, en 100000 muestras)
# Hay varios ESS muu bajos, muy ineficiente.


# Eliptical-contrained vainilla-ABC ----------------------------------------

highest_prob <- 0.98
kernel_thres <- dnorm(qnorm((1-highest_prob) / 2, sd = 0.025), sd = 0.025)
# kernel density at the 98 % cutoff.

# select the half best particles, unconstrain, and fit sn model.
rows_use <- dat$like_overlap >= kernel_thres
xtrain <- dat$par_values[rows_use, ]
ll <- dat$like_overlap[rows_use]
dim(xtrain)
samples_train <- unconstrain2(xtrain, sup)

# filter inf
ids_ok <- sapply(1:nrow(samples_train), function(i) {
  return(all(is.finite(samples_train[i, ])) && !anyNA(samples_train[i, ]))
})
samples_train <- samples_train[ids_ok, ]
ll <- ll[ids_ok]

# Fit skew normal distrib proposal
selm_form <- cbind(intercept, vfi, tfi, slope, wind, steps) ~ 1
snfit <- sn::selm(formula = selm_form,
                  data = as.data.frame(samples_train), family = "SN")
dp <- coef(snfit, "DP", vector = FALSE)
var_factor <- 1
dp$Omega <- dp$Omega * var_factor
sndist <- makeSECdistr(dp, "SN", "proposal", par_names)

# plot(sndist)
densities <- dmsn(samples_train, dp = dp, log = T)
higher_prob2 <- 0.99
thres_dens <- quantile(densities, prob = 1-higher_prob2)
non_extreme <- densities > thres_dens
hist(densities[non_extreme])
plot(ll[non_extreme] ~ densities[non_extreme])


# Sample
ntry <- 1e8
prop <- matrix(rlogis(ntry * n_coef), ntry, n_coef)
dens <- dmsn(prop, dp = dp, log = T)
keep <- dens > thres_dens
nsim <- sum(keep)
prop_sim <- prop[keep, ]

sims <- similarity_simulate_parallel(constrain2(prop_sim, sup),
                                     fire_data = spread_data)
sims$ll <- kernel_expquad3(sims$metrics[, "ov"], overlap_max, ov_scale = 0.025)

# compute weights and neff
w <- sims$ll / sum(sims$ll)
neff <- 1 / sum(w ^ 2)

# Efficience ?
cat("eff first step =", nsim / ntry)
cat("eff second step =", neff / nsim)

wave_plot(sims, "ll", alpha = 0.2)
ids_resample <- sample(1:nrow(sims), size = 10000, prob = sims$ll, replace = T)

sims_boot <- sims[ids_resample, ]
wave_plot(sims_boot, "ll", alpha = 0.05)

pairs(sims_boot$par_values, col = rgb(0, 0, 0, 0.1), pch = 19)



# Tareas ------------------------------------------------------------------

# Lo último es lo que va (lo de vainilla abc eliptical constrained)
# La función muestreadora debería hacer lo siguiente:
# 01) juntar 1000 partículas para simular
# 02) estimar el neff. Si eficiencia < 10 %,
#     subir el umbral de densidad de aceptación, de a poco, para tener al menos
#     15 %.
# 03) una vez logrado eso, simular partículas hasta que el neff
#     llegue a 2000.
# 04) reproducir hasta tener 20000 re-muestras bootstrap, ajustar msn
#     y listo, go to make inference.

### NOOOOO, ignorar eso de la elíptica.
### usar importance-rejection sampling para expensive simulators
### (explicado arriba con sampler avoiding unuseful simulations.)

# Habrá que editar el simulador para que multiplique y sume internamente
# los efectos no direccionales, porque ya no será una constante entre
# réplicas.

# lEER tambieńe lo que puse en esta sección, está re interesting lo que pasa.
# El overlap no varía tanto entre réplicas, estuvimos perdiendo el tiempo.
# (al menos en donde realmente importa. Creo que sí varía mucho cuando le
# dejamos steps altos.)
# Cuánto varía el overlap entre réplicas?
