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
target_dir <- file.path("files", "joint_overlap_estimation")

# load file with constants to standardize
ndvi_params <- readRDS(file.path("data", "NDVI_regional_data",
                                 "ndvi_optim_and_proportion.rds"))

slope_sd <- ndvi_params$slope_term_sd

# constants for fire spread simulation
upper_limit <- 1
n_coef <- 9
n_coef_est <- n_coef - 1
par_names <- c("forest", "shrubland", "grassland",
               "ndvi", "north", "elev", "slope", "wind")

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
                  "wind" = 0)

params_upper <- c("forest" = ext_alpha,
                  "shrubland" = ext_alpha,
                  "grassland" = ext_alpha,
                  "ndvi" = 0,
                  "north" = ext_beta,
                  "elev" = 0,
                  "slope" = ext_beta / slope_sd, # because it is not standardized
                  "wind" = ext_beta)

sup <- rbind(params_lower, params_upper)
rownames(sup) <- c("lower", "upper")

# Parameter space exploration settings

# Sobol phase settings
sobol_n <- 1000
kw <- 4
npw_sobol <- sobol_n / kw

# Maximization settings
npw_maxim <- 100
maxim_n <- 7000   # get 8000 particles in total, 4 + 70 waves
n_best <- 50
vf <- c(1.5, 2, 1) # variance factors for explore_likelihood

# Multicore settings
n_cores <- 15


# Get optimal steps by fire -----------------------------------------------

# sim_files <- list.files(file.path("files", "steps_optimization"),
#                         pattern = "_simulations.rds")
#
# fire_names <- sapply(sim_files, function(x) {
#   strsplit(x, "-steps_optimization_simulations.rds")[[1]][1]
# }) %>% unname
# n_fires <- length(fire_names)
#
# sims <- lapply(sim_files, function(x) {
#   ff <- file.path("files", "steps_optimization", x)
#   return(readRDS(ff))
# })
# names(sims) <- fire_names
#
# steps_optim <- data.frame(fire_id = fire_names,
#                           steps = NA)

# for(i in 1:n_fires) {
#   # i = 1
#   d <- sims[[i]]
#   row_max <- which.max(d$overlap)
#   steps_max <- d$par_values[row_max, "steps"]
#
#   # get ovelap at max steps
#   row_steps_max <- which.max(d$par_values[, "steps"])
#   ov_steps_max <- d$overlap[row_steps_max]
#   ov_max <- d$overlap[row_max]
#
#   ov_diff <- ov_max - ov_steps_max
#
#   if(ov_diff >= 0.01) {
#     steps_optim$steps[i] <- round(d$par_values[row_max, "steps"])
#   } else {
#     rows_good <- which(d$overlap >= (ov_max - 0.01))
#     steps_min <- min(d$par_values[rows_good, "steps"])
#     steps_optim$steps[i] <- round(steps_min)
#   }
#
#   plot(d$overlap ~ d$par_values[, "steps"], ylab = "overlap", xlab = "steps",
#        pch = 19, col = rgb(0, 0 ,0, 0.1), main = fire_names[i])
#   abline(v = steps_optim$steps[i], col = "blue")
#
# } # NOISY ESTIMATION, USE GAM


# for(i in 1:n_fires) {
#   # i = 2
#   d <- sims[[i]]
#
#   # the cpp function rounds below continuous values if integer is required
#   d$par_values[, "steps"] <- floor(d$par_values[, "steps"])
#
#   s_uni <- unique(d$par_values[, "steps"])
#
#   # Get maximum overlaps in a fine sequence of steps
#   cut_points <- seq(min(s_uni) - 0.5, max(s_uni) + 0.5,
#                     length.out = min(length(s_uni), 201))
#   d$xcat <- cut(d$par_values[, "steps"], breaks = cut_points,
#                 include_lowest = T)
#
#   dfsub <- data.frame(overlap = d$overlap,
#                       steps = d$par_values[, "steps"],
#                       xcat = d$xcat)
#   dfagg <- aggregate(cbind(overlap, steps) ~ xcat, dfsub, max)
#
#   # Fit unidimensional GAM by fire
#   nk <- 20
#   kk <- seq(min(s_uni), max(s_uni), length.out = nk)
#
#   mm <- gam(overlap ~
#               s(steps, k = nk, bs = "cr"),
#             data = dfagg, method = "REML", knots = list(parvalues = kk))
#
#   # data to predict
#   dpred <- data.frame(steps = seq(min(s_uni), max(s_uni)))
#   dpred$overlap <- predict(mm, dpred) %>% as.vector
#
#   ov_thres <- max(dpred$overlap) - 0.01
#   # get minimum steps at overlap >= (max - 0.02)
#   dpred_good <- dpred[dpred$overlap >= ov_thres, ]
#   steps_choose <- min(dpred_good$steps)
#
#   ov_thres2 <- max(dpred$overlap) - 0.005
#   # get minimum steps at overlap >= (max - 0.02)
#   dpred_good2 <- dpred[dpred$overlap >= ov_thres2, ]
#   steps_choose2 <- min(dpred_good2$steps)
#
#   steps_choose3 <- dpred$steps[which.max(dpred$overlap)]
#
#   steps_optim$steps[i] <- steps_choose3 # use the fitted max
#
#   # plot(overlap ~ steps, dpred, type = "l", main = fire_names[i])
#   # abline(v = steps_choose, col = "blue")
#   # abline(v = steps_choose2, col = "red")
#   # abline(v = steps_choose3, col = "forestgreen")
# }

# write.csv(steps_optim, "files/steps_optimization/optimal_steps.csv",
#           row.names = F)
steps_optim <- read.csv("files/steps_optimization/optimal_steps.csv")

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

# The same but returning many metrics. Currently used only for steps used.
# It returns a data.frame with averages.
similarity_simulate_particle_size <- function(particle, fire_data = NULL) {

  fire_sim <- simulate_fire_compare(
    landscape = fire_data$landscape[, , -1],
    vegetation = fire_data$landscape[, , 1],
    ignition_cells = fire_data$ig_rowcol,
    coef = particle[1:(n_coef-1)],
    upper_limit = upper_limit,
    steps = particle[n_coef]
  )

  ov <- overlap_spatial(
    fire_sim, fire_data[c("burned_layer", "burned_ids")]
  )
  size <- ncol(fire_sim$burned_ids)

  return(c(ov, size))
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

# Simulate similarity in parallel in all fires. It returns an array with the
# simulated overlaps by particle (row), replicate (col), and fire (slice).
# The wave id is used to give a name to the intermediate results that are written.
similarity_simulate_parallel_af <- function(particles_mat = NULL,
                                            wave = NULL,
                                            fire_ids = 1:length(filenames),
                                            n_sim = n_rep) {

  # ### TESTO
  # np <- 40
  # particles_mat <- matrix(rnorm(np, n_coef_est), np, n_coef_est)
  # wave = 1
  # fire_ids = 1:length(filenames)
  # n_sim = 10
  # ### ENDO TESTO

  # array to be filled with simulations
  n_fires <- length(fire_ids)

  overlap_arr <- array(
    NA, dim = c(nrow(particles_mat), n_sim, n_fires),
    dimnames = list(particle = paste("wave", wave, "_",
                                     1:nrow(particles_mat), sep = ""),
                    replicate = 1:n_sim,
                    fire = 1:n_fires)
  )

  # Loop over fires, saving the result at the end of every one.
  tt <- paste("Wave", wave)
  message(tt)

  for(f in 1:n_fires) {
    full_data <- readRDS(file.path(data_dir, filenames[fire_ids[f]]))
    spread_data <- full_data[c("landscape", "ig_rowcol",
                               "burned_layer", "burned_ids")]

    fire_name <- full_data$fire_id_spread
    message(fire_name)
    # add steps parameter to the particles matrix
    steps_focal <- steps_optim$steps[steps_optim$fire_id == fire_name]
    particles_mat2 <- cbind(particles_mat, steps_focal)

    # turn particle matrix into list for parallel evaluation
    particles_list <- lapply(1:nrow(particles_mat2),
                             function(x) particles_mat2[x, ])

    # simulate fires in parallel
    if(fire_name == "2015_50") {
      registerDoMC(12) # to avoid errors in some cores, which are likely due to ram issues
    } else {
      registerDoMC(n_cores)
    }
    result <- foreach(pp = particles_list) %dopar% {
      similarity_simulate_particle(pp, spread_data, n_sim)
    }

    # rbind list result
    overlap_arr[, , f] <- do.call("rbind", result)
    dimnames(overlap_arr)[[3]][f] <- fire_name

    # save temporary result
    nn <- paste(fire_name, "-simulations-wave-", wave, ".rds", sep = "")
    saveRDS(overlap_arr[, , f], file.path(target_dir, nn))

    rm(full_data, spread_data); gc()
  }

  nn <- paste("all_fires_simulations-wave", wave, ".rds", sep = "")
  saveRDS(overlap_arr, file.path(target_dir, nn))

  # remove intermediate files
  all_files_saved <- list.files(target_dir, pattern = "-simulations-")
  invisible(lapply(all_files_saved, function(f) unlink(file.path(target_dir, f))))

  return(overlap_arr)
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

wave_plot <- function(data, response = "overlap_log_sum", alpha = 0.3, best = NULL,
                      x = "par_values", thres = NULL, bin = FALSE, title = NULL,
                      rc = c(4, 2)) {

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
         xlab = par_names[i], main = mm, #ylim = yy,
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
#' @param sobol logical, indicate whether to initialize the sequence.
#' @return matrix of dimension p x n of samples
rmvn_sobol <- function(n, mu, Sigma, sobol = F, sobol_init = F){
  if(is.null(dim(mu))) {
    mu <- matrix(mu, nrow = 1)
  }

  p <- ncol(mu)

  if(sobol) {
    if(nrow(mu) > 1) stop("Sobol sequence not recommended with more than one center.")
    P <- sobol(n, p, init = sobol_init)
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
logit_scaled <- function(x, support, nc = n_coef_est) {
  xun <- x
  for(j in 1:nc) {
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
scale_params <- function(x, support, nc = n_coef_est) {
  xs <- x
  for(i in 1:nc) {
    xs[, i] <- x[, i] * (support[2, i] - support[1, i]) + support[1, i]
  }
  return(xs)
}


# Function to make data.frame with overlap summarized across replicates and
# fires. It takes the particles matrix and the overlap array.
summarise_overlap <- function(overlap_array, particles) {

  # overlap_array = sobol_arr
  # particles

  # average overlap by fire (aggregate replicates)
  ov_mean <- apply(overlap_array, c(1, 3), mean)

  # data.frame with results
  d <- data.frame(overlap_log_sum = rowSums(log(ov_mean)),
                  overlap_prod = apply(ov_mean, 1, prod))

  d$par_values <- particles
  colnames(d$par_values) <- par_names

  d$ov_values <- ov_mean

  # get wave
  ww <- dimnames(overlap_array)[[1]]
  w <- sapply(ww, function(x) strsplit(substr(x, 5, 10), "_")[[1]][1]) %>% unname %>% as.numeric
  d$wave <- w

  return(d)
}


# Function to reproduce particles for the next wave. As the simulation is
# expensive, it is done before simulating fires, so particles can be checked.
# This is similar to explore_likelihood(), but simpler, and does not simulate
# fires.
# The best particles, defined by n_best, p_best or accept_thres,
# are reproduced with or without weights. If weights are used, the overlap_prod
# is the initial weight, but if less than n_distinct unique particles are
# resampled, weights are iteratively flattened to get at least n_distinct
# unique particles. In this case, n particles are resampled from the n_best
# (or equivalent), and from that sample the MVN kernel is estimated.
# If use_weights = FALSE, the n_best particles (or equivalent) are used to
# create estimate a MVN kernel.
reproduce_particles <- function(data, n = 600,
                                n_best = "all", p_best = NULL, accept_thres = NULL,
                                var_factor = 1, centre = "global",
                                sobol = TRUE, sobol_init = F,
                                use_weights = TRUE, n_distinct = 30,
                                support) {

  ### Testo
  # data = wave1; n = 500; n_best = 100
  # var_factor = 1; centre = "global"; sobol = TRUE
  # use_weights = TRUE; n_distinct = 30
  # support = sup
  # accept_thres = NULL; p_best = NULL
  ### endo testo

  # if threshold is used, ignore p_best
  if(!is.null(accept_thres)) {
    p_best <- NULL
    n_best <- sum(data[, response] >= accept_thres)
  }

  # which particles will be resampled?
  if(!is.null(p_best)) {
    n_best <- round(nrow(data) * p_best)
    if(n_best < 100) {
      n_best <- 500
    }
  }

  if(n_best == "all") {
    n_best <- nrow(data)
  } else {
    data <- data[order(data$overlap_log_sum, decreasing = TRUE), ]
  }

  # resample the best particles
  if(use_weights) {
    # sometimes the best particles are too few, and the computation of the vcov
    # fails. In these cases we need to smooth a bit the weights, so n_distinct
    # different particles are obtained.
    t <- 0
    l <- 0
    while(l < n_distinct) {
      ww <- data$overlap_prod[1:n_best] ^ (1 / (exp(t)))
      ids_rep <- sample(1:n_best, size = n, replace = T, prob = ww)
      l <- length(unique(ids_rep))
      t <- t + 1
    }

    dexp <- data[ids_rep, ]
  } else {
    dexp <- data[1:n_best, ] # using flat weights
  }

  # transform to unconstrained scale
  dexp$par_values_raw <- logit_scaled(dexp$par_values, support)

  # MVN kernel
  # first, remove infinite values and redefine n
  finite <- apply(dexp$par_values_raw, 1, function(x) all(is.finite(x)))
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

  candidates_raw <- rmvn_sobol(n, mus, Vlarge, sobol, sobol_init = sobol_init) %>% t
  # transform to original scale
  candidates <- invlogit_scaled(candidates_raw, support)
  colnames(candidates) <- par_names

  return(candidates)
}


# These functions were used to explore the "likelihood" by fire. As simulations
# are more expensive now, we loop.

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

# Explore joint overlap function -------------------------------------------

# Phase 1: sobol ----------------------------------------------------------

n_sobol <- 3000
n_waves <- 6
n_pw <- n_sobol / n_waves

sobol_list <- vector("list", n_waves)

ss <- sobol(n = n_sobol, dim = n_coef_est)
particles <- scale_params(ss, sup)
waves_sobol <- rep(1:n_waves, each = n_pw)

# for(i in 1:n_waves) {
#   sobol_list[[i]] <- similarity_simulate_parallel_af(
#     particles_mat = particles[waves_sobol == i, ],
#     wave = i
#   )
# }

ff <- list.files(target_dir)
sobol_list <- lapply(ff, function(x) {
  readRDS(file.path(target_dir, x))
})

sobol_arr <- abind(sobol_list, along = 1)
wave1 <- summarise_overlap(sobol_arr, particles)

# wave_plot(wave1[wave1$overlap_log_sum > -100, ])
# wave_plot(wave1[wave1$overlap_log_sum > -100, ], response = "overlap_prod")


# Phase 2: search maximum with flat weights --------------------------------

np <- 3000     # total particles
npw <- 100     # particles per wave
nw <- np / npw # number of waves

var_factors <- rep(c(1, 2, 1.5), ceiling(nw / 3))

# list to save simulation arrays (not cumulative)
search_max_flat_list <- vector("list", nw)

# list saving summarised and cumulative datasets
data_list <- vector("list", nw + 1)
data_list[[1]] <- wave1

previous_waves <- max(wave1$wave)

# Loop by waves, filling the data_list.

for(w in 1:nw) {
  # # test
  # w = 1

  # index in the data_list before this wave
  d = w # because only one phase was ran before (in 6 waves)

  # get previous data
  dprev <- do.call("rbind", data_list[1:d])

  # get new particles from cumulative previous data
  particles <- reproduce_particles(
    data = dprev, n = npw, n_best = 100, var_factor = var_factors[w],
    centre = "global", sobol = TRUE, use_weights = F, support = sup
  )

  # simulate fire in the new particles.
  search_max_flat_list[[w]] <- similarity_simulate_parallel_af(
    particles_mat = particles,
    wave = w + previous_waves
  )

  # summarise simulations and save in data_list
  data_list[[d+1]] <- summarise_overlap(search_max_flat_list[[w]],
                                        particles)

  # print max overlap_log_sum
  mmmax <- max(c(dprev$overlap_log_sum, data_list[[d+1]]$overlap_log_sum))
  mm <- paste("Wave ", w + previous_waves, ". Highest overlap_log_sum = ",
              round(mmmax, 3), sep = "")
  message(mm)

  # write cumulative data_list
  nn <- paste("data_list_cumulative_wave_", w + previous_waves, ".rds", sep = "")
  saveRDS(data_list[1:(d+1)], file.path(target_dir, nn))
}
# endline



# Exploring batch 2 --------------------------------------------------------

batch2 <- readRDS(file.path(target_dir, "data_list_cumulative_wave_36.rds"))
batch2 <- do.call("rbind", batch2)
dim(batch2)

wave_plot(batch2[batch2$overlap_log_sum > -100, ], rc = c(2, 4))

x11()
pairs(batch2$par_values[batch2$overlap_log_sum > -73, ])
dev.off()

# Phase 3: search maximum using weights --------------------------------

np <- 3000     # total particles
npw <- 100     # particles per wave
nw <- np / npw # number of waves

var_factors <- rep(c(1, 2, 1.5), ceiling(nw / 3))

# list to save simulation arrays (not cumulative)
search_max_flat_list <- vector("list", nw)

# list saving summarised and cumulative datasets
data_list <- vector("list", nw + 1)

batch2 <- readRDS(file.path(target_dir, "data_list_cumulative_wave_36.rds"))
data_list[[1]] <- batch2

previous_waves <- max(batch2$wave)

# Loop by waves, filling the data_list.

for(w in 1:nw) {
  # # test
  # w = 1

  # index in the data_list before this wave
  d = w # because all previous waves have been flattened and saved into the
        # first list element.

  # get previous data
  dprev <- do.call("rbind", data_list[1:d])

  # get new particles from cumulative previous data
  particles <- reproduce_particles(
    data = dprev, n = npw, n_best = 100, var_factor = var_factors[w],
    centre = "global", use_weights = T, support = sup,
    n_distinct = 30,
    sobol = F
  )

  # simulate fire in the new particles.
  search_max_flat_list[[w]] <- similarity_simulate_parallel_af(
    particles_mat = particles,
    wave = w + previous_waves
  )

  # summarise simulations and save in data_list
  data_list[[d+1]] <- summarise_overlap(search_max_flat_list[[w]],
                                        particles)

  # print max overlap_log_sum
  mmmax <- max(c(dprev$overlap_log_sum, data_list[[d+1]]$overlap_log_sum))
  mm <- paste("Wave ", w + previous_waves, ". Highest overlap_log_sum = ",
              round(mmmax, 3), sep = "")
  message(mm)

  # write cumulative data_list
  nn <- paste("data_list_cumulative_wave_", w + previous_waves, ".rds", sep = "")
  saveRDS(data_list[1:(d+1)], file.path(target_dir, nn))
}
# endline


# Exploring batch 3 ----------------------------------------------

batch3 <- readRDS(file.path(target_dir, "data_list_cumulative_wave_66.rds"))
batch3 <- do.call("rbind", batch3)
dim(batch3)

wave_plot(batch3[batch3$overlap_log_sum > -80, ], rc = c(2, 4))
wave_plot(batch3,
          response = "overlap_prod", alpha = 0.1,
          rc = c(2, 4))
wave_plot(batch3,
          alpha = 0.1,
          rc = c(2, 4))

max(batch3$overlap_log_sum)
thres <- -70

x11()
pairs(batch3$par_values[batch3$overlap_log_sum > thres, ])
# batante corr

# improvement?
par(mfrow = c(2, 4))
for(vv in par_names) {
# vv = "ndvi"
filt1 <- batch3$overlap_log_sum > -100 &
         batch3$wave <= 36
plot(batch3$overlap_log_sum[filt2] ~ batch3$par_values[filt2, vv],
     pch = 19, col = rgb(0, 0, 0, 0.1),
     ylab = "overlap log sum", xlab = vv)

filt2 <- batch3$overlap_log_sum > -100 &
         batch3$wave > 36
points(batch3$overlap_log_sum[filt1] ~ batch3$par_values[filt1, vv],
     pch = 19, col = rgb(1, 0, 0, 0.1))

}
par(mfrow = c(1, 1))
# nice

# overlap histogram in the best particle
max_row <- which.max(batch3$overlap_log_sum)
hist(batch3$ov_values[max_row, ], breaks = 10)

batch3$overlap_mean <- rowMeans(batch3$ov_values)
batch3$overlap_median <- apply(batch3$ov_values, 1, median)
plot(overlap_log_sum ~ overlap_mean, batch3[batch3$overlap_log_sum > -70, ])
plot(overlap_log_sum ~ overlap_median, batch3[batch3$overlap_log_sum > -80, ])

# another round with higher var_factors?
sqrt(c(1, 1.5, 2, 3))

curve(dnorm(x), from = -25, to = 15, n = 300)
facs <- c(3, 6)
for(f in facs) {
  curve(dnorm(x, sd = f), add = T, col = f)
}

var_factors <- c(1, 3, 6) ^ 2


# Size difference for the best particles ----------------------------------

nbest <- 50
nrep <- 10
nmet <- 2
n_fires <- length(filenames)
fire_names <- character(n_fires)

size_sim <- array(
  NA, dim = c(nbest, nmet, nrep, n_fires),
  dimnames = list(
   particle = 1:nbest,
   met = c("overlap", "size"),
   replicate = 1:nrep,
   fire = 1:n_fires
  )
)
fire_size <- numeric(n_fires)

good <- order(batch3$overlap_log_sum, decreasing = T)[1:nbest]
ppbest <- batch3$par_values[good, ]

for(f in 1:n_fires) {
  full_data <- readRDS(file.path(data_dir, filenames[f]))
  spread_data <- full_data[c("landscape", "ig_rowcol",
                             "burned_layer", "burned_ids")]

  fire_name <- full_data$fire_id_spread
  print(fire_name)
  fire_names[f] <- fire_name
  fire_size[f] <- ncol(full_data$burned_ids)

  # add steps
  steps_focal <- steps_optim$steps[steps_optim$fire_id == fire_name]
  particles_copy <- cbind(ppbest, steps_focal)

  for(p in 1:nbest) {
    for(r in 1:nrep) {
      size_sim[p, , r, f] <- similarity_simulate_particle_size(
        particle = particles_copy[p, ],
        fire_data = spread_data
      )
    }
  }
}

saveRDS(size_sim,
        file.path(target_dir, "size_sim_array_after_batch3.rds"))


# turn size into size quotients relative to the observed

size_size_sim <- size_sim[, "size", , ]

sizeq <- size_size_sim
for(f in 1:n_fires) {
  sizeq[, , f] <- sizeq[, , f] / fire_size[f]
}

avg1 <- apply(sizeq, c(1, 3), mean)
sizeqmean <- colMeans(avg1)
hist(sizeqmean)
abline(v = median(sizeqmean), lty = 2)
abline(v = mean(sizeqmean), lty = 2)

hist(as.vector(avg1), breaks = 10)
hist(as.vector(sizeq))

ov_sim <- size_sim[, "overlap", , ]
hist(as.vector(ov_sim))
mean(ov_sim)

avg1ov <- apply(ov_sim, c(1, 3), mean)
avg2ov <- apply(avg1ov, 2, mean)
hist(avg2ov)
sizeqmean <- colMeans(avg1)
summary(ov_sim)

# Phase 4: search maximum, wider var ----------------------------------

