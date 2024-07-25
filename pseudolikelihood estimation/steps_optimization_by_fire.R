# This code inherits from
# <pseudolikelihood_estimation_ig-known.R>.

# Here I simply try to explore the overlap function for each fire to get a
# point estimate of the steps parameter. The overlap function exploration scheme
# is modified to spend less time in areas with low overlap. Trials were made in
# <single fire posterior - optimizing steps by fire.R>

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())

library(Rcpp)
library(terra)
library(abind)         # manage arrays

library(randtoolbox)   # sobol sequences
library(MGMM)          # fit MVN distribution
library(Matrix)        # nearPD, searches the closest PD matrix.
library(stutils)       # tryNULL, for search of nearPD

library(foreach)       # parallelization
library(doMC)          # doRNG had problems with cpp functions

library(FireSpread)    # spread and similarity functions

source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for rast_from_mat and a few constants

# Constants --------------------------------------------------------------

data_dir <- file.path("data", "focal fires data", "landscapes")
filenames <- list.files(data_dir)

# dir to save output
target_dir <- file.path("files", "steps_optimization")

# load file with constants to standardize
ndvi_params <- readRDS(file.path("data", "NDVI_regional_data",
                                 "ndvi_optim_and_proportion.rds"))

slope_sd <- ndvi_params$slope_term_sd

# constants for fire spread simulation
upper_limit <- 1
n_coef <- 9
par_names <- c("forest", "shrubland", "grassland",
               "ndvi", "north", "elev", "slope", "wind",
               "steps")

n_rep <- 20 # fires simulated by particle

ext_alpha <- 50
ext_beta <- 30

params_lower <- c("forest" = -ext_alpha,
                  "shrubland" = -ext_alpha,
                  "grassland" = -ext_alpha,
                  "ndvi" = -ext_beta,
                  "north" = 0,
                  "elev" = -ext_beta,
                  "slope" = 0,
                  "wind" = 0,
                  "steps" = 5)

params_upper <- c("forest" = ext_alpha,
                  "shrubland" = ext_alpha,
                  "grassland" = ext_alpha,
                  "ndvi" = 0,
                  "north" = ext_beta,
                  "elev" = 0,
                  "slope" = ext_beta / slope_sd, # because it is not standardized
                  "wind" = ext_beta,
                  "steps" = NA)       # to be filled later, varies by fire

# Parameter space exploration settings

# Sobol phase settings
sobol_n <- 1000
kw <- 4
npw_sobol <- sobol_n / kw

# Maximization settings
npw_maxim <- 100
maxim_n <- 7000   # get 8000 particles in total, 4 + 70 waves
n_best <- 50
vf <- c(1.5, 2, 1, 3) # variance factors for explore_likelihood

# Multicore settings
n_cores <- 15

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

  for(i in 1:n_sim) {
    # simulate and compute metrics
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

  return(colMeans(metrics))
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
get_bounds <- function(fire_data, n_sim = 1) {

  coef_burn_all <- numeric(n_coef) # intercepts, slopes, steps
  coef_burn_all[1:3] <- 1e9  # intercepts
  coef_burn_all[n_coef] <- 0 # steps (0 means no limit)

  coef_burn_none <- numeric(n_coef) # intercepts, slopes, steps
  coef_burn_none[1:3] <- 1e-9  # intercepts
  coef_burn_none[n_coef] <- 1  # steps

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
  # fails. In these cases we need to smooth a bit the weights, so >= 30
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

# Loop over fires -------------------------------------------------------

for(f in 1:length(filenames)) {

  fire_file <- filenames[f]
  fire_name <- strsplit(fire_file, ".rds")[[1]]
  print(paste("Fire:", fire_name))

  if(fire_name == "2015_50") {
    registerDoMC(12) # to avoid errors in some cores, which are likely due to ram issues
  } else {
    registerDoMC(n_cores)
  }

  full_data <- readRDS(file.path(data_dir, fire_file))
  # subset data needed for spread (to be cloned across workers)
  spread_data <- full_data[c("landscape", "ig_rowcol",
                             "burned_layer", "burned_ids")]

  bb <- get_bounds(fire_data = spread_data)
  params_upper["steps"] <- bb["largest", "steps_used"]
  sup <- rbind(params_lower, params_upper)

  # Sobol phase
  ss <- sobol(n = sobol_n, dim = n_coef)
  particles <- scale_params(ss, sup)
  waves_1 <- rep(1:kw, each = npw_sobol)

  # loop to save intermediate results
  wave1 <- NULL
  print("Sobol phase.")
  for(w in 1:kw) {
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
    # print(paste("wave ", w, ", overlap max: ", rr, sep = ""))

    mm <- paste("wave ", w, ", overlap max: ", rr, sep = "")
    message(mm)
  }

  print("Maximum search.")

  wave2 <- explore_likelihood_iterate(
    n_waves = maxim_n / npw_maxim, n = npw_maxim, n_best = n_best,
    var_factor = vf, data = wave1,
    support = sup, spread_data = spread_data,
    centre = "global", sobol = TRUE, phase = "search_maximum",
    write = T, fire_name = fire_name
  )

  # save final result
  nn <- paste(fire_name, "-steps_optimization_simulations.rds", sep = "")
  saveRDS(wave2, file.path(target_dir, nn))

  # save wave plot
  nn <- paste(fire_name, "-steps_optimization_plot.png", sep = "")
  png(filename = file.path(target_dir, nn),
      width = 15, height = 16, units = "cm", res = 300)
  wave_plot(wave2, title = fire_name, alpha = 0.2)
  dev.off()

  # remove intermediate files
  all_files_saved <- list.files(target_dir, pattern = "-wave-")
  invisible(lapply(all_files_saved, function(f) unlink(file.path(target_dir, f))))

  # clean
  remove(full_data, spread_data, wave1, wave2)
  gc()
}
#endline



# Exploring results -------------------------------------------------------

# obtain the mean overlap in the best 50 particles of each fire.

sim_files <- list.files(target_dir, pattern = "-steps_optimization_simulations.rds")

ov_mean <- sapply(sim_files, function(x) {
  d <- readRDS(file.path(target_dir, x))
  ids <- order(d$overlap, decreasing = T)[1:50]
  return(mean(d$overlap[ids]))
})

hist(ov_mean, xlim = c(0, 1), breaks = 7)
abline(v = mean(ov_mean), lty = 2)
