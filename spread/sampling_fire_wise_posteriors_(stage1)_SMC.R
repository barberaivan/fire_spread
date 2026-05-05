# This code inherits from
# <sampling_fire_wise_posteriors_(stage1).R>.

# Here I apply the SMC instead of the simple rejection sampling used in the
# dissertation, based on Del Moral et al. 2011: An adaptive sequential Monte 
# Carlo method for approximate Bayesian computation
# [Pierre Del Moral · Arnaud Doucet · Ajay Jasra]

# _02 uses a hard abc kernel

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)

library(Rcpp)
library(terra)
library(abind)         # manage arrays

library(randtoolbox)   # sobol sequences

library(posterior)     # manage posterior samples
library(tidybayes)     # not sure if this was used
library(bayesplot)     # visualize posteriors

library(foreach)       # parallelization
library(doMC)          # doRNG had problems with cpp functions

library(FireSpread)    # spread and similarity functions

library(microbenchmark)

source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for rast_from_mat and a few constants
 
sourceCpp(file.path("spread", "sample_triplets_weighted.cpp"))

# source("estimation_functions.R") # prior_dist and other stuff

# Constants --------------------------------------------------------------

data_dir <- file.path("data", "focal fires data", "landscapes")
filenames <- list.files(data_dir)
n_fires <- length(filenames)

# dir to save output
target_dir <- file.path("files", "posterior_samples_stage1_smc")

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

support_base <- rbind(params_lower, params_upper)

# ABC-SMC constants

# Number of particles
N <- 10000

# Proportion at which neff is reduced at each iteration
alpha <- 0.5

# Scale factor for the DE-MCMC differences
g <- 1

# Scale for DE-MCMC jitter
sigma_jitter <- 0.001

# SMC termination parameters
acc_thres <- 0.015  # minimum acceptance rate for MCMC
ndist_thres <- 300  # minimum number of distinct particles after MCMC move
maxit <- 100        # maximum number of iterations
edelta <- 0.0001    # decrease in e is considered small if < edelta
nstall <- 10        # max number of consecutive iters with small e allowed

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

## CAREFUL: below another variable was added to size_data (steps_upper)

# Functions ---------------------------------------------------------------

normalize <- function(x) x / sum(x)

# Simulate fires and compare them with the observed one using the overlap_spatial
# function. The landscape argument includes all data to simulate fire, and also
# burned and burned_ids layers.
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

# Simulate fires and compare then with the observed one using the overlap_spatial
# function. The landscape argument includes all data to simulate fire, and also
# burned and burned_ids layers.
# Returns a matrix with
#   overlap_spatial,
#   size_diff: size differences, as simulated - observed,
#   edge: number of pixels burned at the edge of the landscape.
# The last two are used to rule out bounded fires.
# Currently used only for steps used.
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

# Compute the spatial overlap for all the particles in a matrix, in pararllel.
# particles_mat: matrix with simulation coefficients (in columns) for all 
#   particles (rows), in the constrained (simulation) scale.
# fire_data: data to simulate fire.
similarity_simulate_parallel <- function(
    particles_mat = NULL,
    fire_data = NULL
) {
  
  stopifnot(is.matrix(particles_mat))
  
  # turn particle matrix into list
  particles_list <- lapply(
    seq_len(nrow(particles_mat)),
    function(i) particles_mat[i, ]
  )
  
  # run in parallel
  result <- foreach(pp = particles_list) %dopar% {
    similarity_simulate_particle(pp, fire_data)
  }
  
  ov <- unlist(result)
  
  res <- data.frame(
    overlap = ov
  )
  
  res$par_values <- particles_mat
  
  return(res)
}

# get_bounds: computes the metrics for the largest and smallest fires possible
get_bounds <- function(fire_data) {

  coef_burn_all <- c(1e6, rep(0, n_coef - 1)) # steps = 0 means infinite
  coef_burn_none <- c(-1e6, rep(0, n_coef - 2), 1)

  small_fire <- similarity_simulate_particle_metrics(
    coef_burn_none, fire_data = fire_data
  )

  large_fire <- similarity_simulate_particle_metrics(
    coef_burn_all, fire_data = fire_data
  )

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

# Transform from unconstrained [logit, log] scale to constrained one, with
# compact support. Uses invlogit scaled for all parameters except for steps, 
# for which it uses exp().
constrain <- function(xun, support) {
  xc <- xun
  
  names_logit <- colnames(xun)[colnames(xun) != "steps"]
  for(j in names_logit) {
    xc[, j] <- plogis(xun[, j]) * (support[2, j] - support[1, j]) + support[1, j]
  }
  
  xc[, "steps"] <- exp(xun[, "steps"])
  
  return(xc)
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



# abc_kernel: compute abc acceptance probability. It's zero if the overlap is 
#   below a threshols, and it's the overlap if it's >= thres.
# d: disimilarity index among simulated and observed datasets (1 - overlap here).
# e: d value below which prob = overlap, and above which the probability is
#   zero.
# log: return log prob?  
abc_kernel <- function(d, e, log = TRUE) {
  good <- as.numeric(d <= e)
  k <- (1 - d) * good # k = overlap contrained to be >= (1 - e)
  if(log) return(log(k)) else return(k)
}

# neff_abc: compute neff based on the quotient of the abc_kernel using two 
# epsilon (Del Moral et al. 2009).
# d: distance metric.
# e_new: new epsilon to try.
# e_old: old epsilon to use.
neff_abc <- function(d, e_old, e_new) {
  k_old <- abc_kernel(d, e_old, log = T)
  k_new <- abc_kernel(d, e_new, log = T)
  
  # log-weights
  logw <- k_new - k_old
  
  # log-sum-exp normalization
  logw <- logw - max(logw)
  w <- exp(logw)
  w <- w / sum(w)
  
  # neff
  neff <- 1 / sum(w^2)
  return(neff)
}

# bisect_e: Finds the new e (epsilon) to update the abc_kernel targettin 
# neff ~= alpha * N.
# d: vector of distances.
# alpha: desired proportion at which NEFF is reduced.
# N: number of particles.
# e_min: minimum value of abc_tolerance to choose.
# e_max: maximum tolerance. At first iteration, it's 1, then, it decreases.
# tol: tolerance in difference |neff - N * alpha|.
# maxit: maximum number of iterations.
bisect_e <- function(
    d, alpha, N, e_max = 1, tol = 1e-4, maxit = 100
) {
  
  lo <- min(d) + 0.0001
  hi <- e_max
  neff_target <- alpha * N
  
  for (i in seq_len(maxit)) {
    mid <- 0.5 * (lo + hi)
    
    # to compute the neff, assume the old value is hi
    neff_try <- neff_abc(d, e_old = e_max, e_new = mid)
    
    if (abs(neff_try - neff_target) < tol)
      return(list(e_new = mid, neff = neff_try))
    
    if (neff_try < neff_target) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  
  list(e_new = mid, neff = neff_try)
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

# Sample fire-wise posteriors ---------------------------------------------

for(f in 1:n_fires) {
  write_file <- ifelse(f > 30, TRUE, FALSE)

  fire_file <- size_data$file[f]
  fire_name <- size_data$fire_id[f]
  
  mm <- paste(
    "Fire ", fire_name, " [", f, "/", n_fires, "]", 
    " ------------------------------------",
    sep = ""
  )
  message(mm)
  
  full_data <- readRDS(file.path(data_dir, fire_file))
  # subset data needed for spread (to be cloned across workers)
  spread_data <- full_data[c(
    "landscape", "ig_rowcol",
    "burned_layer", "burned_ids",
    "counts_veg"
  )]

  bb <- get_bounds(fire_data = spread_data)
  sup <- rbind(params_lower, params_upper)
  sup[2, "steps"] <- bb["largest", "steps_used"]

  if(fire_name %in% c("2015_50", "1999_2140469994_r")) {
    registerDoMC(10)
  } else {
    registerDoMC(15)
  }
  
  ## Sequential Monte Carlo

  # Initialization

  # Sample the prior using a sobol sequence. The prior is a unit-scale, zero-mean
  # logistic distribution, flat over the constrained space.
  ss <- sobol(n = N, dim = n_coef)
  colnames(ss) <- par_names
  
  # particles will be the matrix of particles (par values) in the unconstrained
  # scale
  particles <- qlogis(ss) # unconstrained parameters
  
  # Simulate dissimilarity over initial particles
  particles_cons <- invlogit_scaled(particles, sup)
  sim <- similarity_simulate_parallel(particles_cons, spread_data)
  d <- 1 - sim$overlap
  
  # Set initial epsilon at 1 (max dissimilarity)
  e <- 1
  
  # Initial (log)-k
  k <- rep(0, N)
  
  acc_rate <- 1
  i <- 1
  ncalm <- 0
  
  # List to save all results
  out <- vector("list", maxit)
  
  # SMC loop
  go_on <- TRUE
  while(go_on) {
    # check N
    if(length(d) < N) stop("Some particles were lost in simulation.")
    
    ## Update e, targetting neff = N * alpha
    # d_sorted <- sort(d)
    # e_new <- d_sorted[round(N * alpha) - 1]
    bibi <- bisect_e(d, alpha, N, e_max = e)
    e_new <- bibi$e_new
    
    ## Update weights
    k_new <- abc_kernel(d, e_new, log = T)
    logw <- k_new - k
    
    # log-sum-exp normalization
    logw <- logw - max(logw)
    w <- exp(logw) |> normalize()
    
    ## DE-MCMC move
    
    # 1) Make threeplets of particle ids, as focal + 2 neighbours
    ids_move <- sample_triplets_weighted(w)
    
    # 2) Create source and proposals
    from <- particles[ids_move[, 1], ]
    p1 <- particles[ids_move[, 2], ]
    p2 <- particles[ids_move[, 3], ]
    jitter <- matrix(rnorm(N * n_coef, 0, sigma_jitter), N, n_coef)
    to <- from + g * (p1 - p2) + jitter
    
    # As from is created from a resampling based on w, with replacement, 
    # new weights are all 1, so they will not be explicitly included in the 
    # update of weights at next iteration.
    
    # Simulate fires in proposal
    ov_to <- similarity_simulate_parallel(
      particles_mat = invlogit_scaled(to, sup),
      fire_data = spread_data
    )$overlap
    
    # Compute metrics to decide acceptance
    d_to <- 1 - ov_to
    k_to <- abc_kernel(d_to, e_new, log = T)
    prior_to <- rowSums(dlogis(to, log = T))
    
    k_from <- k_new[ids_move[, 1]]
    prior_from <- rowSums(dlogis(from, log = T))
    
    # log posterior
    lp_to <- k_to + prior_to
    lp_from <- k_from + prior_from
    
    # Metropolis decision using log-sum-exp trick
    # Accept if log(U) < (to - from)
    log_u <- log(runif(N))
    
    accept <- log_u < (lp_to - lp_from)
    acc_rate <- sum(accept) / N
    
    # Update particles, dissimilarity, k, and epsilon
    particles <- from
    particles[accept, ] <- to[accept, ]
    
    d <- d[ids_move[, 1]] # d at source particles
    d[accept] <- d_to[accept] # replace d at accepted particles
    
    e_diff <- e - e_new
    e <- e_new
    
    k <- k_from
    k[accept] <- k_to[accept]
    
    # Number of unique particles
    ndist <- sum(!duplicated(particles))
    
    ## Inform acceptance rate and epsilon
    mm <- paste(
      "Iteration ", i, "\n",
      "  acceptance: ", round(acc_rate, 4), "\n", 
      "  distinct: ", ndist, "\n",
      "  epsilon: ", round(e, 4), 
      sep = ""
    )
    message(mm)
    
    # Save iteration result
    res <- list(
      iteration = i,
      acc_rate = acc_rate,
      distinct = ndist,
      epsilon = e,
      particles = particles, # unconstrained, all-logit scale!
      overlap = 1 - d
    )
    out[[i]] <- res
    
    # Update iteration
    i <- i + 1
    
    # continue?
    ncalm <- ifelse(e_diff < edelta, ncalm + 1, 0)
    
    go_on <- (
      acc_rate >= acc_thres & 
      ndist >= ndist_thres & 
      i <= maxit & 
      ncalm <= nstall
    )
    
    # Save partial or full result
    if (go_on) {
      if (write_file) { # save partial result
        fname <- paste(
          "partial_samples_", fire_name, "__", i-1, ".rds", sep = ""
        )
        fdir <- file.path(
          "files", "posterior_samples_stage1_smc", fname
        )
        saveRDS(res, fdir)
      }
    } else {
      # save out and remove all partial files
      fname <- paste(
        "full_samples_history_", fire_name, ".rds", sep = ""
      )
      fdir <- file.path(
        "files", "posterior_samples_stage1_smc", fname
      )
      saveRDS(out[1:(i-1)], fdir)
      
      # remove partials
      dirr <- file.path(
        "files", "posterior_samples_stage1_smc"
      )
      fn <- list.files(dirr)
      fn <- fn[grep("partial_", fn)]
      unlink(file.path(dirr, fn))
    }
  }
  
  # clean
  remove(full_data, spread_data, out, res)
  gc()
}

# start at 26/01/2026, 11:40 h
# ended at 
# ~2.4 days
 
# Particles are stored at logit scale, including steps (not log)

# Merge samples in a single array -----------------------------------------

target_dir <- file.path("files", "posterior_samples_stage1_smc")
samples_files <- list.files(target_dir, pattern = "full_samples_history_")
J1 <- length(samples_files)

fire_names_0 <- sapply(samples_files, function(i) {
  strsplit(i, split = "full_samples_history_")[[1]][2]
})
fire_ids <- sapply(fire_names_0, function(i) {
  strsplit(i, split = ".rds")[[1]][1]
}) |> unname()

# Make a list with matrices of samples
samples_list <- vector("list", J1)
names(samples_list) <- fire_ids

for (i in 1:J1) {
  print(i)
  
  rrr <- readRDS(file.path(target_dir, samples_files[i]))
  rr <- rrr[[length(rrr)]] # get latest iteration
  
  out <- rr$particles
  finite <- sapply(1:nrow(out), function(j) all(is.finite(out[j, ])))
  out <- out[finite, ]
  
  if (nrow(out) != N) {
    warning("Missing samples")
  }
  
  # get support to scale parameters
  sup_local <- support_base
  sup_local[2, "steps"] <- size_data$steps_upper[size_data$fire_id == fire_ids[i]]
  
  # transform from all-logit to constrained space
  out_cons <- invlogit_scaled(out, sup_local) 
  # transform to unconstrained (logit-log) space
  out_unc <- unconstrain(out_cons, sup_local)
  
  samples_list[[i]] <- out_unc
}

# Rows removed?
all(sapply(samples_list, function(x) dim(x))[1, ] == N)

# Turn into 3D array
samples_arr <- abind::abind(samples_list, along = 3)
str(samples_arr)
dimnames(samples_arr) <- list(
  iter = 1:N,
  param = par_names,
  fire_id = fire_ids
)

saveRDS(samples_arr, file.path(target_dir, "samples_all_fires.rds"))


# Check output structure with previous file --------------------------------

smc_out <- readRDS(file.path(target_dir, "samples_all_fires.rds"))
imp_out <- readRDS(file.path("files", "posterior_samples_stage1", "samples_boot_all_fires.rds"))
str(smc_out)
str(imp_out)
all(dimnames(smc_out)[[3]] == dimnames(imp_out)[[3]])
# OK

# Compare overlap among SMC and importance sampling -----------------------

ov_mean_smc <- sapply(1:J1, function(i) {
  # i = 1
  rrr <- readRDS(file.path(target_dir, samples_files[i]))
  rr <- rrr[[length(rrr)]]
  return(mean(rr$overlap))
})

target_dir_imp <- file.path("files", "posterior_samples_stage1")
samples_files_imp <- list.files(target_dir_imp, pattern = "-samples.rds")

ov_mean_imp <- sapply(1:J1, function(i) {
  rrr <- readRDS(file.path(target_dir_imp, samples_files_imp[i]))
  return(rrr$overlap_mean)
})

plot(ov_mean_smc ~ ov_mean_imp)
abline(0, 1)
# OK, slightly higher for SMC
mean(ov_mean_smc - ov_mean_imp)
range(ov_mean_smc - ov_mean_imp)