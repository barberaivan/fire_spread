# Code inherited from
# <single fire posterior - surrogate vs simulator - binary gam.R>.

# Trying to approximate single-fire similarity functions through models was too
# hard. The most accurate approximation was the binary GAM, where the probability
# of finding overlap values above some threshold was modelled. However,
# overlap functions were quite flat in most dimensions, so accepting or not based
# on a threshold could lead to unreasonable inferences. In particular, there were
# parameters with high correlation, hard to distinguish in small fires. For example,
# if a small fire burned only forest and at a north slope, the northing was estimated
# negative and the forest, at the highest possible value.
# Then, I decided to abandon the single-fire approach.

# To resolve the problem of steps varying by fires, I will try to take it as a
# known parameters, fixed at the minimum number of steps required to burn all
# the burned pixels. Then, a separate model to predict the steps of future
# fires will be fitted separately. The posterior distribution of the remaining
# parameters will be defined conditioning on the known steps, which vary by fire.

# I will try to use a similarity function defined as sum(log(overlap_f)), treating
# it as a log-likelihood.
# If this metric shows a well-behaved shape, maybe we could approximate it with
# a model, to perform the MCMC later.

# The joint similarity will be explored by reproducing the good particles, in an
# approach similar to the used before.

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

# Multicore settings -----------------------------------------------------

n_cores <- 15
registerDoMC(n_cores)

# Constants --------------------------------------------------------------

data_dir <- file.path("data", "focal fires data", "landscapes")
filenames <- list.files(data_dir)

# dir to save output
target_dir <- file.path("files", "pseudolikelihood_estimation")

# constants for fire spread simulation
upper_limit <- 1
n_coef <- 8
par_names <- c("forest", "shrubland", "grassland",
               "ndvi", "north", "elev", "slope", "wind")

# number of fires to simulate by particle
n_rep <- 20

gam_formula_bern <- formula(
  y ~
    # marginal effects
    s(forest, bs = basis, k = k_side) +
    s(shrubland, bs = basis, k = k_side) +
    s(grassland, bs = basis, k = k_side) +
    s(ndvi, bs = basis, k = k_side) +
    s(north, bs = basis, k = k_side) +
    s(elev, bs = basis, k = k_side) +
    s(slope, bs = basis, k = k_side) +
    s(wind, bs = basis, k = k_side) +

    # interactions (intercepts)
    ti(forest, ndvi, k = k_int, bs = basis) +
    ti(forest, north, k = k_int, bs = basis) +
    ti(forest, elev, k = k_int, bs = basis) +
    ti(forest, slope, k = k_int, bs = basis) +
    ti(forest, wind, k = k_int, bs = basis) +

    ti(shrubland, ndvi, k = k_int, bs = basis) +
    ti(shrubland, north, k = k_int, bs = basis) +
    ti(shrubland, elev, k = k_int, bs = basis) +
    ti(shrubland, slope, k = k_int, bs = basis) +
    ti(shrubland, wind, k = k_int, bs = basis) +

    ti(grassland, ndvi, k = k_int, bs = basis) +
    ti(grassland, north, k = k_int, bs = basis) +
    ti(grassland, elev, k = k_int, bs = basis) +
    ti(grassland, slope, k = k_int, bs = basis) +
    ti(grassland, wind, k = k_int, bs = basis) +


    ti(elev, slope, k = k_int, bs = basis) +
    ti(elev, wind, k = k_int, bs = basis) +
    ti(elev, ndvi, k = k_int, bs = basis) +
    ti(slope, wind, k = k_int, bs = basis) +
    ti(north, ndvi, k = k_int, bs = basis)
)


# Functions ---------------------------------------------------------------

normalize <- function(x) x / sum(x)

# Obtain the minimum number of steps required to burn all the burnable pixels,
# in a setting with fire spread probability = 1.
# Returns the list of fire data with the landscape trimmed to remove areas
# that will not be used because the steps does not allow to burn further.

get_steps <- function(spread_data) {

  ## test
  # spread_data <- fire1

  # duplicate fire data
  spread_data2 <- spread_data

  steps_count <- simulate_fire_animate(
    landscape = spread_data$landscape[, , -1],
    vegetation = spread_data$landscape[, , 1],
    ignition_cells = spread_data$ig_rowcol,
    coef = c(1e9, 1e9, 1e9, 0, 0, 0, 0, 0),
    steps = 0
  )

  # steps used to burn all the burned cells
  steps_masked <- spread_data$burned_layer * steps_count
  steps_max <- max(steps_masked)

  # find rows-cols range used
  lte_max <- steps_count <= steps_max

  cols_seq <- 1:ncol(spread_data$landscape)
  cols_used <- colSums(lte_max)
  cols_range <- range(cols_seq[cols_used > 0])

  rows_seq <- 1:nrow(spread_data$landscape)
  rows_used <- rowSums(lte_max)
  rows_range <- range(rows_seq[rows_used > 0])

  # add variables to duplicated fire data
  spread_data2$steps <- steps_max
  spread_data2$rows_range <- rows_range
  spread_data2$cols_range <- cols_range

  return(spread_data2)
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
  # simulate and compute metrics
  for(i in 1:n_sim) { # simulated fires by particle
    fire_sim <- simulate_fire_compare(
      landscape = fire_data$landscape[, , -1],
      vegetation = fire_data$landscape[, , 1],
      ignition_cells = fire_data$ig_rowcol,
      coef = particle[1:(n_coef-1)],
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

  # simulate and compute metrics
  for(i in 1:n_sim) { # simulated fires by particle
    fire_sim <- simulate_fire_compare(
      landscape = fire_data$landscape[, , -1],
      vegetation = fire_data$landscape[, , 1],
      ignition_cells = fire_data$ig_rowcol,
      coef = particle[1:(n_coef-1)],
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
                      x = "par_values", thres = NULL, bin = FALSE, title = NULL,
                      rc = c(3, 3)) {

  par_names <- colnames(data$par_values)
  n_coef <- length(par_names)

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
  if(is.null(title)) {
    title <- deparse(substitute(data))
  }

  par(mfrow = c(rc[1], rc[2]))
  for(i in 1:n_coef) {
    mm <- ifelse(i == 1, title, NA)
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
  # data = sim1; n = 500; p_best = 0.15; prob_fn = NULL;
  # support = sup; spread_data = spread_data;
  # centre = "global"; sobol = TRUE;
  # pow = 1; var_factor = 1
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
  while(l < 30) {
    ids_rep <- sample(1:n_best, size = n, replace = T,
                      prob = weight[1:n_best] ^ (1 / (1 + k)))
    l <- length(unique(ids_rep))
    k <- k + 0.5
  }

  dexp <- data[ids_rep, ]

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
drop_uncertain <- function(unc, p) {
  return(p * plogis(73 - 90 * unc))
}

# sample from the GAM taking independent draws (simulation)
# unc_thres [0, 1] determines the maximum uncertainty allowed in the model.
#   if higher, the particle is rejected.
rejection_sample <- function(iter, model, support) {

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
    prob <- predict(gam_bern, part_df, type = "response")

    # only evaluate uncertainty at particles with high probability.
    # (This step is time consuming)
    unc <- numeric(length(prob))
    high_prob_ids <- which(prob > 0.02)
    if(length(high_prob_ids) > 0) {
      pp <- predict(gam_bern, part_df[high_prob_ids, ], se.fit = T)
      lll <- plogis(pp$fit - qnorm(0.975) * pp$se.fit)
      uuu <- plogis(pp$fit + qnorm(0.975) * pp$se.fit)
      unc[high_prob_ids] <- uuu - lll
    }

    # drop probability if highly uncertain
    prob2 <- drop_uncertain(unc, prob)
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
rejection_sample_parallel <- function(iter, model, support, cores = 15) {
  registerDoMC(cores)
  ii <- as.list(1:cores)
  foreach(chain = ii) %dopar% {
    rejection_sample(iter, model, support)
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


# Compare fixed steps with previous optimum -------------------------------

# fire_name <- strsplit(fire_file, ".rds")[[1]]

spread_data_list <- lapply(1:length(filenames), function(f) {
  dd <- readRDS(file.path(data_dir,  filenames[f]))
  dd2 <- get_steps(dd)
  return(dd2)
})

tabb <- data.frame(
  fire_id = lapply(spread_data_list, function(x) x[["fire_id_spread"]]) %>% unlist,
  steps = lapply(spread_data_list, function(x) x[["steps"]]) %>% unlist,
  steps_mean = 0,
  steps_min = 0,
  steps_max = 0,
  overlap_max = 0,
  overlap_match = 0
)

# import simulations to get optimal steps and other stuff

simulation_files <- list.files(file.path("files", "pseudolikelihoods_OLD"),
                               pattern = "-simulations_data")

simulations <- lapply(simulation_files, function(f) {
  readRDS(file.path("files", "pseudolikelihoods_OLD", f))
})

names(simulations) <- sapply(simulation_files, function(x) {
  strsplit(x, "-simulations_data.rds")[[1]][1]
})

for(i in 1:nrow(tabb)) {
  # i = 1
  simdata <- simulations[[tabb$fire_id[i]]]

  # choose 100 better steps
  simdata <- simdata[order(simdata$overlap, decreasing = T), ]

  # get fitted steps
  filt <- simdata$overlap >= max(simdata$overlap) - 0.05

  good_steps <- simdata$par_values[filt, "steps"]

  tabb$steps_mean[i] <- mean(good_steps)
  tabb$steps_min[i] <- min(good_steps)
  tabb$steps_max[i] <- max(good_steps)
  tabb$overlap_max[i] <- max(simdata$overlap)

  # get maximum overlap in a step sequence, to get an estimate of the
  # optimal
  simdata <- simdata[order(simdata$par_values[, "steps"]), ]
  simdata2 <- data.frame(steps = simdata$par_values[, "steps"],
                         overlap = simdata$overlap)
  ng <- 100
  nrow(simdata2) / ng
  simdata2$steps_class <- rep(1:100, each = nrow(simdata2) / ng)
  simdata_agg <- aggregate(overlap ~ steps_class, simdata2, max)
  simdata_agg$steps <- aggregate(steps ~ steps_class, simdata2, mean)[, "steps"]

  m1 <- gam(overlap ~ s(steps, k = ng / 2, bs = "cr"), data = simdata_agg)

  png(filename = paste("pseudolikelihood estimation/steps-matching-figures/",
                       tabb$fire_id[i], ".png", sep = ""),
      width = 10, height = 8, units = "cm", res = 300)
  plot(overlap ~ steps, simdata_agg, main = tabb$fire_id[i],
       pch = 19, col = rgb(0, 0, 0, 0.8),
       ylim = c(0, max(simdata$overlap) * 1.05))
  lines(fitted(m1) ~ steps, simdata_agg, col = "blue")
  abline(v = tabb$steps[i], col = "red")
  dev.off()

  tabb$overlap_match[i] <- predict(m1, data.frame(steps = tabb$steps[i]))
}

rr <- range(c(tabb$steps_max, tabb$steps_min, tabb$steps))

ggplot(tabb, aes(steps, steps_mean, ymin = steps_min, ymax = steps_max)) +
  geom_linerange(alpha = 0.8) +
  geom_point(alpha = 0.8, size = 2) +
  geom_abline(slope = 1, intercept = 0, color = viridis(1, begin = 0.2)) +
  coord_fixed() +
  scale_y_continuous(limits = c(0, rr[2]), expand = c(0.01, 0.01)) +
  scale_x_continuous(limits = c(0, rr[2]), expand = c(0.01, 0.01)) +
  ylab("Optimal steps") +
  xlab("Fixed steps")

ggsave("pseudolikelihood estimation/steps_compare_max-005.png",
       width = 12, height = 12, units = "cm")

ggplot(tabb, aes(overlap_match, overlap_max)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_abline(slope = 1, intercept = 0, color = viridis(1, begin = 0.2)) +
  coord_fixed() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  ylab("Overlap at optimal steps") +
  xlab("Overlap at fixed steps")

ggsave("pseudolikelihood estimation/steps_compare_max-005_overlaps.png",
       width = 12, height = 12, units = "cm")

tabb$steps_diff <- tabb$steps - tabb$steps_mean
tabb$overlap_diff <- tabb$overlap_match - tabb$overlap_max

ggplot(tabb, aes(steps_diff, overlap_diff)) +
  geom_hline(yintercept = 0, color = viridis(1, begin = 0.5, option = "B"),
             linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = viridis(1, begin = 0.5, option = "B"),
             linetype = "dashed", linewidth = 0.5) +
  geom_point(alpha = 0.8, size = 2.5) +
  xlab("Steps difference\n(fixed - optimal)") +
  ylab("Overlap difference\n(fixed steps - optimal steps)")

ggsave("pseudolikelihood estimation/steps_compare_diff-diff.png",
       width = 14, height = 12, units = "cm")


# Y hacemos los wave_plots de ejemplo

png(filename = paste("pseudolikelihood estimation/wave-plots-old-model/",
                     "1999_2140469994_r", ".png", sep = ""),
    width = 15, height = 11, units = "cm", res = 300)

wave_plot(simulations[["1999_2140469994_r"]], alpha = 0.1,
          title = "1999_2140469994_r", rc = c(2, 3))
dev.off()

png(filename = paste("pseudolikelihood estimation/wave-plots-old-model/",
                     "2015_53", ".png", sep = ""),
    width = 15, height = 11, units = "cm", res = 300)

wave_plot(simulations[["2015_53"]], alpha = 0.1, title = "2015_53",
          rc = c(2, 3))
dev.off()

png(filename = paste("pseudolikelihood estimation/wave-plots-old-model/",
                     "2015_50", ".png", sep = ""),
    width = 15, height = 11, units = "cm", res = 300)

wave_plot(simulations[["2015_50"]], alpha = 0.1, title = "2015_50",
          rc = c(2, 3))
dev.off()



# OLD CODE BELOW ----------------------------------------------------------
# Simulation waves ---------------------------------------------------------

ext_alpha <- 30
ext_beta <- 30
ext_ndvi <- 200


fire_file <- "2008_3.rds" #"2000_8.rds"#filenames[f] # #"1999_25j.rds"#"2015_53.rds"#
fire_name <- strsplit(fire_file, ".rds")[[1]]
print(paste("Fire:", fire_name))
full_data <- readRDS(file.path(data_dir, fire_file))
# subset data needed for spread (to be cloned across workers)
spread_data <- full_data[c("landscape", "ig_rowcol",
                           "burned_layer", "burned_ids")]

bb <- get_bounds(fire_data = spread_data)

# params_lower <- c("forest" = -ext_alpha,
#                   "shrubland" = -ext_alpha,
#                   "grassland" = -ext_alpha,
#                   "ndvi" = -ext_ndvi,
#                   "north" = 0,
#                   "elev" = -ext_beta,
#                   "slope" = 0,
#                   "wind" = 0,
#                   "steps" = 5)
#
# params_upper <- c("forest" = ext_alpha,
#                   "shrubland" = ext_alpha,
#                   "grassland" = ext_alpha,
#                   "ndvi" = 0,
#                   "north" = ext_beta * 2,
#                   "elev" = 0,
#                   "slope" = ext_beta * 2,
#                   "wind" = ext_beta,
#                   "steps" = bb["largest", "steps_used"])

## free sign support:

params_lower <- c("forest" = -ext_alpha,
                  "shrubland" = -ext_alpha,
                  "grassland" = -ext_alpha,
                  "ndvi" = -ext_ndvi,
                  "north" = -ext_beta * 2,
                  "elev" = -ext_beta,
                  "slope" = -ext_beta * 2,
                  "wind" = -ext_beta,
                  "steps" = 5)

params_upper <- c("forest" = ext_alpha,
                  "shrubland" = ext_alpha,
                  "grassland" = ext_alpha,
                  "ndvi" = ext_ndvi,
                  "north" = ext_beta * 2,
                  "elev" = ext_beta,
                  "slope" = ext_beta * 2,
                  "wind" = ext_beta,
                  "steps" = bb["largest", "steps_used"])

sup <- rbind(params_lower, params_upper)


# Explore likelihood
# Phase 1: sobol
registerDoMC(15)

ss <- sobol(n = 10000, dim = n_coef)
particles <- scale_params(ss, sup)
waves_1 <- rep(1:20, each = 500)
# loop to save
wave1 <- NULL
print(paste("Phase: sobol. Fire: ", fire_name, sep = ""))
for(w in 1:20) {
  # w = 1
  wlocal <- similarity_simulate_parallel(particles[waves_1 == w, ],
                                         spread_data)
  wlocal$wave <- w
  wlocal$phase <- "sobol"
  colnames(wlocal$par_values) <- par_names

  # nn <- paste(fire_name, "-wave-", w, ".rds", sep = "")
  # saveRDS(wlocal, file.path(target_dir, nn))

  wave1 <- rbind(wave1, wlocal)
  rr <- round(max(wave1$overlap), 4)
  print(paste("wave ", w, ", overlap max: ", rr, sep = ""))
}

# wave 2: reproduce 15 % best (500 * 20)
print(paste("Phase: search higher. Fire: ", fire_name, sep = ""))
wave2 <- explore_likelihood_iterate(
  n_waves = 20, data = wave1, n = 500, p_best = 0.15,
  var_factor = c(1.5, 2, 1),
  support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE, phase = "search_higher",
  write = F, fire_name = fire_name
)
# wave_plot(wave2)

# 5000 iter reproducing the 500 best, only to reach the max
print(paste("Phase: search maximum. Fire: ", fire_name, sep = ""))
wave3 <- explore_likelihood_iterate(
  n_waves = 10, data = wave2, n = 500, n_best = 500,
  var_factor = c(1.5, 2, 1),
  support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE, phase = "search_maximum",
  write = F, fire_name = fire_name
)

# 15000 iter reproducing the 200 best, only to reach the max
print(paste("Phase: search maximum. Fire: ", fire_name, sep = ""))
wave4 <- explore_likelihood_iterate(
  n_waves = 30, data = wave3, n = 500, n_best = 100,
  var_factor = c(1.5, 2, 1),
  support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE, phase = "search_maximum",
  write = F, fire_name = fire_name
)

thres <- max(wave4$overlap) - 0.05
# wave_plot(wave4, alpha = 0.1, thres = thres)

# 1000 iter above thres
print(paste("Phase: above threshold. Fire: ", fire_name, sep = ""))
wave5 <- explore_likelihood_iterate(
  n_waves = 20, data = wave4, n = 500, accept_thres = thres,
  var_factor = c(1.5, 2, 1),
  support = sup, spread_data = spread_data,
  centre = "global", sobol = TRUE, phase = "above_threshold",
  write = F, fire_name = fire_name
)


like_sim <- wave5

wave_plot(like_sim[like_sim$overlap > 0.5, ], alpha = 0.1,
          thres = thres)


# Fit gam -----------------------------------------------------------------

data_gam_bern <- as.data.frame(cbind(like_sim$par_values,
                                     y = as.numeric(like_sim$overlap >= thres)))

k_side <- 15
k_int <- 6
basis <- "cr"

gam_bern <- bam(
  gam_formula_bern, data = data_gam_bern, family = "binomial",
  method = "fREML", discrete = T, nthreads = 8
)
summary(gam_bern)

fbern <- fitted(gam_bern)

# Sample GAM-approximated posterior

# print(paste("Sampling posterior. Fire: ", fire_name, sep = ""))
sup_reduced <- reduce_support(like_sim$par_values[like_sim$overlap >= thres, ],
                              sup, 0.5)
r_gam <- rejection_sample_parallel(1000, gam_bern, sup_reduced, cores = 15)
draws <- do.call("rbind", r_gam) %>% as.data.frame

par(mfrow = c(3, 3))
for(i in 1:ncol(draws)) {
  dd <- density(draws[, i], from = sup[1, i], to = sup[2, i])
  plot(dd, main = names(draws)[i],
       ylim = c(0, max(dd$y) * 1.05),
       xlim = sup[, i], ylab = NA, xlab = NA)
}
par(mfrow = c(1, 1))

# MCMC simulating fire ----------------------------------------------------

# initial values (fixed, so all chains start equal)
nc <- 20
best_100 <- order(like_sim$overlap, decreasing = TRUE)[1:200]
set.seed(23432)
best_ids <- sample(best_100, nc, replace = F)
best_parts <- logit_scaled(like_sim$par_values[best_ids, ],
                           sup)
init_list <- lapply(1:nc, function(i) best_parts[i, ])
vv <- logit_scaled(like_sim$par_values[like_sim$overlap >= thres, ],
                   support = sup)
# remove Inf from vv
finites <- apply(vv, 1, function(x) all(is.finite(x)))
vv <- vv[finites, ]
sigma_init <- cov(vv)

sampling_iters <- 8000
adapt_iters <- 2000

r_sim <- MCMC_parallel(fun = fn_like_sim_bin,
                       n = sampling_iters + adapt_iters,
                       adapt = adapt_iters,
                       n_chains = nc,
                       n_cores = 15,
                       scale = sigma_init, init_list = init_list,
                       support = sup, fire_data = spread_data, thres = thres)
gc()
# con 2008_3 tardó mucho y lo cancelé


# Comparo posteriores -----------------------------------------------------

# postprocessing
post_names <- c("simulation", "gam") # "gam_bin"
arr_all <- list(tidy_samples(r_sim, adapt_iters, sup),
                tidy_samples_ind(r_gam))
names(arr_all) <- post_names

# Compare with densities
dflong <- do.call("rbind", lapply(arr_all, function(a) {
  # tt <- thin_draws(a, 10)
  return(as_draws_df(a, .nchains = nc))
}))
names(dflong) <- par_names
dflong$sampling <- factor(
  rep(post_names, each = nrow(dflong) / length(post_names)),
  levels = post_names
)

dflong <- dflong[, c(par_names, "sampling")]
dflonger <- pivot_longer(dflong, all_of(1:n_coef), values_to = "par_value",
                         names_to = "par_names")
dflonger$par_names <- factor(dflonger$par_names, levels = par_names)

ggplot(dflonger, aes(x = par_value, fill = sampling, color = sampling)) +
  geom_density(alpha = 0, adjust = 2) +
  facet_wrap(vars(par_names), scales = "free") +
  theme(panel.grid.minor = element_blank()) +
  viridis::scale_color_viridis(discrete = TRUE, end = 0.7, option = "D")

# Miramos por cadena el MCMC
mcmc_dens_overlay(arr_all[["simulation"]]) +
  scale_color_viridis(discrete = TRUE, option = "A", end = 0.8) +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("simulation")

# posteriores de a pares

combs <- as.data.frame(combn(par_names, 2) %>% t)
names(combs) <- c("x", "y")
ncomb <- nrow(combs)
plist <- vector("list", ncomb)
dflong2 <- dflong
dflong2$sampling <- factor(dflong$sampling, levels = levels(dflong$sampling),
                           labels = c("simulation", "gam"))

for(i in 1:ncomb) {
  # i = 1
  take <- combs[i, ] %>% as.character()
  dlocal <- dflong2[, c(take, "sampling")]

  p <-
  ggplot(dlocal, aes_string(x = take[1], y = take[2])) +
    geom_hdr(probs = seq(0.05, 0.95, by = 0.1),
             method = method_kde(adjust = 2)) +
    facet_wrap(vars(sampling), nrow = 1) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          strip.background = element_rect(fill = "white", color = "white"),
          strip.text = element_text(size = 10))


  if(!(i %in% 1:6)) {
    p <- p +
      theme(strip.background = element_blank(),
            strip.text = element_blank())
  }

  # if(i == 9) {
  #   p <- p + theme(legend.position = "right")
  # }

  plist[[i]] <- p
}

x11()
pairs_compare <- egg::ggarrange(plots = plist, ncol = 6,
                                draw = FALSE, newpage = FALSE)
ggsave(file.path(target_dir, paste(fire_name, "-plot_pairs_mcmc-gam.png", sep = "")),
       plot = pairs_compare, height = 26, width = 34, units = "cm")
dev.off()


# uniform product ---------------------------------------------------------

up <- sapply(1:10000, function(i) {
  prod(runif(57))
})
plot(density(up))

up <- sapply(1:10000, function(i) {
  sum(log(runif(57)))
})

plot(density(up))
plot(density(exp(up)))
log(0.2) * 57

sum(log(0.3) * 57)
sum(log(0.8) * 57)

ll_ov <- function(x = 0.5, n = 57) log(x) * n

curve(ll_ov(x, 57), from = 0.01, to = 1)
curve(log(x) * 57, from = 0.01, to = 1)

# Ideas para definir steps independiente de lo demás ----------------------

firesim_count <- simulate_fire_animate(
  landscape = spread_data$landscape[, , -1],
  vegetation = spread_data$landscape[, , 1],
  ignition_cells = spread_data$ig_rowcol,
  coef = c(1e9, 1e9, 1e9, 0, 0, 0, 0, 0),
  steps = 0
)
range(firesim_count)

# en qué step todo lo quemado llega a quemarse?
sss <- sapply(1:max(firesim_count), function (s) {
  # s = 1
  spread_data$burned_layer
  mstep <- firesim_count <= s & firesim_count > 0
  mmatch <- sum(mstep * spread_data$burned_layer)

  return(mmatch)
})

mean(like_sim$par_values[best_100, "steps"])
range(like_sim$par_values[best_100, "steps"])

plot(sss)
abline(v = range(like_sim$par_values[like_sim$overlap >= thres, "steps"]))

firesim_count <= 1


