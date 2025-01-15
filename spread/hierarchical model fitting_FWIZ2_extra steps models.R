# Inherits from <hierarchical model fitting.R>. Now, simply use standardized
# FWI (cumulative, as before), and sample longer, fewer chains

# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(deeptime)
theme_set(theme_bw())

library(bayesplot)
library(bayestestR)    # hdi
library(posterior)
library(abind)

library(terra)

library(LaplacesDemon) # rinvwishart
library(invgamma)
library(MASS)          # mvrnorm
library(truncnorm)
library(truncreg)      # tnorm regression, for inits
library(truncnorm)

library(mgcv)          # models to tune proposals
library(DHARMa)        # model assessment

library(FireSpread)    # simulate fires to test model
library(foreach)       # parallelization
library(doMC)          # parallelization

library(rstan)         # exploratory logitnorm steps model
library(logitnorm)     # logtinorm model for steps
source("spread/mcmc_functions.R")

# Functions ---------------------------------------------------------------

summarise <- function(x) {
  q <- quantile(x, probs = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
                method = 8)
  names(q) <- c("eti_lower_95", "eti_lower_90", "eti_lower_80", "median",
                "eti_upper_80", "eti_upper_90", "eti_upper_95")

  hdi_95 <- hdi(x, ci = 0.95)
  hdi_90 <- hdi(x, ci = 0.90)
  hdi_80 <- hdi(x, ci = 0.80)

  hdis <- c(
    "hdi_lower_95" = hdi_95$CI_low, "hdi_lower_90" = hdi_90$CI_low,
    "hdi_lower_80" = hdi_80$CI_low,
    "hdi_upper_80" = hdi_80$CI_high, "hdi_upper_90" = hdi_90$CI_high,
    "hdi_upper_95" = hdi_95$CI_high
  )

  res <- c("mean" = mean(x), hdis, q)

  return(res)
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

constrain_vec <- function(x, v, support) {
  if(v != "steps") {
    return(plogis(x) * (support[2, v] - support[1, v]) + support[1, v])
  } else {
    return(exp(x))
  }
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

# The same but vectorized (v), use for steps parameters
logit_scaledv <- function(x, L, U) {
  qlogis((x - L) / (U - L))
}

# and its inverse
invlogit_scaledv <- function(x, L, U) {
  plogis(x) * (U - L) + L
}

# MCMC: sampler for the joint posterior

# nsim: number of iterations to run.
# sd_jump: list containing the proposal sd for each parameter that is updated
#   through m-h, using either normal or truncated normal in the case of sigma.
#   These parameters are coefficients for area ~ steps regression and the steps
#   of fires for which spread was not simulated (unknown ignition point).
# start: list with initial values.
# samples: list with the output of a previous mcmc run, used to get starting
#   values.
# jump_factor: in the case of missing(sd_jump), the proposal sigmas are computed
#   from the MLE distribution, but multiplying them by the jump_factor.
# sd_jump_out: return the used sd_jump list?

# returns a list with 3 arrays of samples: fixef, ranef and steps.
mcmc <- function(nsim = 50, thin = 1, sd_jump, start, samples, jump_factor = 4,
                 sd_jump_out = F, progress = F) {

  # #_____ TEST
  # nsim = 10; thin = 1; jump_factor = 4; sd_jump_out = F
  # #_____ End test

  #### Allocate memory to save samples,
  nsave <- floor(nsim / thin)
  # replace nsim to avoid computing iterations that won't be save
  nsim <- nsave * thin

  # fixef:       [parname, partype (a, b, s2), nsave]
  # rho:         [parname, parname,            nsave] # correlation matrix for ranef
  # ranef:       [parname, fire_id,            nsave]
  # steps:       [1,       fire_id,            nsave]

  # fixef include the area ~ steps regression parameters.

  fixef_save <- array(NA, dim = c(n_coef + 1, 3, nsave),
                      dimnames = list(
                        par_names = par_names_all,
                        par_class = c("a", "b", "s2"),
                        iter = 1:nsave
                      ))

  rho_save <- array(NA, dim = c(n_coef, n_coef, nsave),
                    dimnames = list(
                      par_names = par_names,
                      par_names = par_names,
                      iter = 1:nsave
                    ))

  ranef_save <- array(NA, dim = c(n_coef, nfires_spread, nsave),
                      dimnames = list(
                        par_names = par_names,
                        fire_id = fire_ids,
                        iter = 1:nsave
                      ))

  steps_save <- array(NA, dim = c(nfires_nonspread, nsave),
                            dimnames = list(
                              fire_id = fire_ids_nonspread,
                              iter = 1:nsave
                            ))

  #### Define proposals sd (if missing)
  if(missing(sd_jump)) {

    tmp <- par_start$fixef["area", , ]
    tmp[3, ] <- sqrt(tmp[3, ]) # sd_jump for sigma, not sigma2
    area_jump <- apply(tmp, 1, sd) * jump_factor
    steps_jump <- apply(par_start$steps, 1, sd) * jump_factor

    if(sd_jump_out) {
      sd_jump <- list(area = area_jump,
                      steps = steps_jump)
    }

  } else {
    area_jump <- sd_jump$area
    steps_jump <- sd_jump$steps
  }

  #### Define initial values (if missing)

  if(!missing(start)) {
    fixef <- start$fixef
    ranef <- start$ranef
    steps <- start$steps
  }

  # If there is no start, but there is a previous sample, set as start
  # the last value of the previous sample
  if(missing(start) & !missing(samples)) {
    nn <- dim(samples$fixef)[3]
    fixef <- samples$fixef[, , nn]
    ranef <- samples$ranef[, , nn]
    steps <- samples$steps[, nn]
  }

  # if no start nor samples are provided, make random start from the MLEs
  if(missing(start) & missing(samples)) {
    is <- sample(1:5000, size = 1)
    fixef <- par_start$fixef[, , is]
    ranef <- par_start$ranef[, , is]
    steps <- par_start$steps[, is]
  }

  #### Transient matrices to store expected value of random effects
  mu <- ranef
  mu[,] <- NA
  mu_steps <- steps
  mu_steps[] <- NA

  # MCMC loop -------------------------------------------------------------
  for(k in 1:nsim) {

    # Spread fixed effects ------------------------------------------------
    # (i.e., not area parameters)
    for(v in 1:n_coef) {

      if(v < n_coef) {
        y <- ranef[v, ]
        X <- X_
        tXX <- tXX_
      } else {
        # for steps, merge data for spread and non-spread fires
        y <- c(ranef[v, ], steps)
        X <- Xlong
        tXX <- tXXlong
      }

      fixef[v, ] <- update_lm(
        y = y, X = X, tXX = tXX, s2 = fixef[v, "s2"],
        b0 = b0[, v], S0_inv = S0_inv, t0 = t0, d0 = d0
      )

      # update mu
      mu_temp <- X %*% fixef[v, 1:2]
      if(v < n_coef) {
        mu[v, ] <- mu_temp
      } else {
        mu[v, ] <- mu_temp[ids1]
        mu_steps <- mu_temp[ids2]
      }
    }

    # Area ~ steps parameters -------------------------------------------------
    steps_all <- c(ranef[n_coef, ], steps)
    aa <- n_coef+1
    fixef[aa, ] <- update_truncnorm(
      y = area_all, x = steps_all, coef = fixef[aa, ], L = areaL,
      b0 = b0[, aa], S0 = S0, t0 = t0, d0 = d0,
      sd_jump = area_jump
    )

    # Correlation matrix -------------------------------------------------
    error_mat <- t(ranef - mu)
    rho <- update_corr(error_mat)
    sdmat <- diag(sqrt(fixef[1:n_coef, 3]))
    Sigma <- sdmat %*% rho %*% sdmat

    # Random effects (Lunn method) ---------------------------------------
    ranef <- update_ranef(
      Y = ranef, mu = mu, Sigma = Sigma, mu0 = mu0, Sigma0 = Sigma0,
      Ytry = Ytry, steps_bounds = steps_bounds,
      area = area1, areaL = areaL, area_coef = fixef[n_coef+1, ],
      s = n_coef
    )

    # Random effects (steps) ---------------------------------------------
    steps <- update_steps(
      steps, mu_steps, s = sqrt(fixef[n_coef, 3]),
      area = area2, areaL, area_coef = fixef[n_coef+1, ], sd_jump = steps_jump
    )

    #### Save samples if thin iterations have passed
    if(k %% thin == 0) {
      if(progress) print(k)
      s <- k / thin
      fixef_save[, , s] <- fixef
      rho_save[, , s] <- rho
      ranef_save[, , s] <- ranef
      steps_save[, s] <- steps
    }
  }

  #### Merge samples
  out <- list(
    fixef = fixef_save,
    rho = rho_save,
    ranef = ranef_save,
    steps = steps_save
  )

  if(sd_jump_out) out$sd_jump <- sd_jump

  return(out)
}

# function to run mcmc in parallel. n_cores is the number of cores and chains.
# the result is the same as for mcmc, but the arrays have a fourth dimension,
# the chain.
# Chains are initialized from pre-selected samples of a long mcmc run.
mcmc_parallel <- function(nsim = 50, thin = 1, n_cores = 8, sd_jump,
                          start_samples) {

  # ### TESTO
  # n_cores <- 8
  # iii <- sample(1:dim(run0_thin$fixef)[3], size = n_cores, replace = F)
  # nsim = 50; thin = 1; n_cores = 8; sd_jump = sd_jump_tune
  # start_samples <- list(
  #   fixef = run0_thin$fixef[, , iii],
  #   ranef = run0_thin$ranef[, , iii],
  #   steps = run0_thin$steps[, iii]
  # )
  # ###
  registerDoMC(n_cores)

  # turn starting values into list
  start_list <- vector("list", n_cores)
  for(cc in 1:n_cores) {
    ll <- list(
      fixef = start_samples$fixef[, , cc],
      ranef = start_samples$ranef[, , cc],
      steps = start_samples$steps[, cc]
    )
    start_list[[cc]] <- ll
  }

  runs <- foreach(ss = start_list) %dopar% {
    mcmc(nsim = nsim, thin = thin, sd_jump = sd_jump, start = ss)
  }

  # extract lists
  fixef_l <- vector("list", n_cores)
  rho_l <- vector("list", n_cores)
  ranef_l <- vector("list", n_cores)
  steps_l <- vector("list", n_cores)

  for(cc in 1:n_cores) {
    fixef_l[[cc]] <- runs[[cc]]$fixef
    rho_l[[cc]] <- runs[[cc]]$rho
    ranef_l[[cc]] <- runs[[cc]]$ranef
    steps_l[[cc]] <- runs[[cc]]$steps
  }

  # tidy runs
  fixef <- abind::abind(fixef_l, along = 4)
  rho <- abind::abind(rho_l, along = 4)
  ranef <- abind::abind(ranef_l, along = 4)
  steps <- abind::abind(steps_l, along = 3)

  dimnames(fixef) <- c(dimnames(fixef_l[[1]]), list("chain" = as.character(1:n_cores)))
  dimnames(rho) <- c(dimnames(rho_l[[1]]), list("chain" = as.character(1:n_cores)))
  dimnames(ranef) <- c(dimnames(ranef_l[[1]]), list("chain" = as.character(1:n_cores)))
  dimnames(steps) <- c(dimnames(steps_l[[1]]), list("chain" = as.character(1:n_cores)))

  out <- list(
    fixef = fixef,
    rho = rho,
    ranef = ranef,
    steps = steps
  )

  return(out)
}

# function to count succesive changes in a vector
count_changes <- function(x) sum(abs(diff(x)) > 0)

# Compute the acceptance rate for the m-h updated parameters.
acceptance <- function(samples) {
  nn <- dim(samples$fixef)[3]
  nt <- nn - 1 # transitions
  parea <- apply(samples$fixef[n_coef+1, , , drop = F], 1:2, count_changes) / nt
  psteps <- unname(apply(samples$steps, 1, count_changes)) / nt
  out <- list(area = parea, steps = psteps)
  return(out)
}

thin <- function(samples, rate = 100) {
  nn <- dim(samples$fixef)[3]
  if(nn < (rate * 10)) return(samples)

  ids <- rev(seq(nn, 1, by = -rate))

  out <- list(
    fixef = samples$fixef[, , ids],
    rho = samples$rho[, , ids],
    ranef = samples$ranef[, , ids],
    steps = samples$steps[, ids]
  )
  return(out)
}

# Functions to simulate fires (assess model)

# function to simulate a fire and compute the overlap and size, including size
# by veg type.
simulate_metrics <- function(particle, fire_data = NULL) {

  metrics <- numeric(nmet)
  names(metrics) <- met_names

  coef_intercepts <- rep(particle[1], n_veg)
  coef_nd <- particle[2:3]
  coef_terrain <- particle[4:(n_coef-1)]
  steps <- particle[n_coef]

  fire_sim <- simulate_fire_compare_veg(
    layer_vegetation = fire_data$landscape[, , "veg"],
    layer_nd = fire_data$landscape[, , nd_variables],
    layer_terrain = fire_data$landscape[, , terrain_variables],
    coef_intercepts = coef_intercepts,
    coef_nd = coef_nd,
    coef_terrain = coef_terrain,
    ignition_cells = fire_data$ig_rowcol,
    upper_limit = upper_limit,
    steps = steps,
    n_veg = n_veg
  )

  metrics["overlap"] <- overlap_spatial(
    fire_sim, fire_data[c("burned_layer", "burned_ids")]
  )
  metrics["size"] <- sum(fire_sim$counts_veg)
  metrics[3:nmet] <- fire_sim$counts_veg

  return(metrics)
}

# The same, in parallel
simulate_metrics_parallel <- function(particles, # matrix with parameter vectors (rows)
                                      fire_data = NULL) {

  particles_list <- lapply(1:nrow(particles), function(x) particles[x, ])

  # simulate in parallel
  result <- foreach(pp = particles_list) %dopar% {
    simulate_metrics(pp, fire_data = fire_data)
  }

  # inherit previous dimnames
  names_single <- names(result[[1]])

  # rbind list result
  res <- do.call("rbind", result)
  colnames(res) <- names_single

  return(res)
}

nice_theme <- function() {
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),

    axis.line = element_line(linewidth = 0.3),

    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),

    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "white", color = "white")
  )
}

# Climatic data ------------------------------------------------------------

# fwi_data <- read.csv(file.path("data", "climatic_data_by_fire_fwi-fortnight-cumulative_FWIZ.csv"))
fwi_data <- read.csv(file.path("data", "climatic_data_by_fire_fwi-fortnight-cumulative_FWIZ2.csv"))

# Fire polygons (to get area) ---------------------------------------------

ff <- vect("data/patagonian_fires_spread.shp")
ff$area_ha <- expanse(ff) / 1e4 # turn m2 to ha

# Constants ---------------------------------------------------------------

# constants for fire spread simulation
upper_limit <- 1
n_veg <- 5
veg_names <- c("wet", "subalpine", "dry", "shrubland", "grassland")
veg_levels <- c("Wet forest", "Subalpine forest", "Dry forest", "Shrubland", "Grassland")

n_terrain <- 2
terrain_names <- c("slope", "wind")
terrain_variables <- c("elevation", "wdir", "wspeed")
n_nd <- n_fi <- 2        # flammability indices
nd_variables <- c("vfi", "tfi")

par_names <- c("intercept", nd_variables, terrain_names, "steps")
n_coef <- length(par_names)

par_names_all <- c(par_names, "area")

n_par <- length(par_names_all)
n_pt <- 3 # b0, b1, s2

# load file with constants to standardize
ndvi_params <- readRDS(file.path("data", "flammability indices",
                                 "ndvi_optim_and_proportion.rds"))
slope_sd <- ndvi_params$slope_term_sd

# support for parameters
ext_alpha <- 50
ext_beta <- 30

params_lower <- c(-ext_alpha, rep(0, n_coef-2), 5)
params_upper <- c(ext_alpha, rep(ext_beta, n_coef-2), NA)
names(params_lower) <- names(params_upper) <- par_names
params_upper["slope"] <- ext_beta / slope_sd

support <- rbind(params_lower, params_upper)
colnames(support) <- names(params_lower) <- names(params_upper) <- par_names
support_width <- apply(support, 2, diff)

# flammability indices parameters
fi_params <- readRDS(file.path("data", "flammability indices",
                               "flammability_indices.rds"))

# summary of predictors that make the flammability indices
data_summ <- readRDS(file.path("data", "flammability indices",
                               "ndvi_elevation_summary.rds"))

# to simulate fires and check model
lands_dir <- file.path("data", "focal fires data", "landscapes")
nmet <- n_veg + 2 # size by veg, size total, overlap
met_names <- c("overlap", "size",
               "wet", "subalpine", "dry", "shrubland", "grassland")

# Load and prepare data ---------------------------------------------------

# data with steps bounds
size_data <- readRDS(file.path("data", "focal fires data", "fire_size_data.rds"))
rownames(size_data) <- size_data$fire_id

# dir to load files
target_dir <- file.path("files", "posterior_samples_stage1")
Ytry <- readRDS(file.path(target_dir, "samples_boot_all_fires.rds"))
Ytry <- aperm(Ytry, c(2, 3, 1))
N1 <- dim(Ytry)[3] # number of samples from stage 1
fire_ids <- dimnames(Ytry)[[2]]

# merge with FWI data
# two fires were split, but they have the same FWI.
fires_data_0 <- fwi_data[!(fwi_data$fire_id %in% c("2011_19", "2015_47")), ]
fires_data_1 <- fwi_data[fwi_data$fire_id %in% c("2011_19", "2015_47"), ]
fires_data_2 <- fires_data_1[c(1, 1, 2, 2), ]

fires_data_2$fire_id[fires_data_2$fire_id == "2011_19"] <-
  fire_ids[grep("2011_19", fire_ids)]

fires_data_2$fire_id[fires_data_2$fire_id == "2015_47"] <-
  fire_ids[grep("2015_47", fire_ids)]

# bring area of separated fires
for(i in 1:nrow(fires_data_2)) {
  fires_data_2$area_ha[i] <- ff$area_ha[ff$fire_id == fires_data_2$fire_id[i]]
}

# put together all fires
fires_data <- rbind(fires_data_0, fires_data_2)
rownames(fires_data) <- NULL

# add log area and scaled fwi
fires_data$area_ha_log <- log(fires_data$area_ha)

fwi_mean <- mean(fires_data$fwi_fort_expquad) # more than zero
fwi_sd <- sd(fires_data$fwi_fort_expquad)     # less than 1

### export fwi scale
# fwi_spread_mean_sd <- list(fwi_mean = fwi_mean,
#                            fwi_sd = fwi_sd)
# saveRDS(fwi_spread_mean_sd,
#         file.path("files", "hierarchical_model_FWIZ", "fwi_mean_sd_spread.rds"))
###

fires_data$fwi <- (fires_data$fwi_fort_expquad - fwi_mean) / fwi_sd
rownames(fires_data) <- fires_data$fire_id

fires_data_spread <- fires_data[fires_data$fire_id %in% fire_ids, ]
fires_data_spread <- fires_data_spread[fire_ids, ]

fires_data_nonspread <- fires_data[!(fires_data$fire_id %in% fire_ids), ]

fire_ids_nonspread <- fires_data_nonspread$fire_id
nfires_nonspread <- length(fire_ids_nonspread)
nfires_spread <- length(fire_ids)

# lower bound for fire area (10 ha)
areaL <- log(10-1e-6)

# N of random effects
J1 <- nfires_spread
J2 <- nfires_nonspread
J <- J1 + J2
ids1 <- 1:J1
ids2 <- (J1+1):J

# variables for both subsets of fires
fwi1 <- fires_data_spread$fwi
fwi2 <- fires_data_nonspread$fwi
fwi_all <- c(fwi1, fwi2)

area1 <- fires_data_spread$area_ha_log
area2 <- fires_data_nonspread$area_ha_log
area_all <- c(area1, area2)

# design matrices
X_ <- cbind(rep(1, J1), fwi1)
tXX_ <- t(X_) %*% X_

Xsteps <- cbind(rep(1, J2), fwi2)
Xlong <- rbind(X_, Xsteps)
tXXlong <- t(Xlong) %*% Xlong

# order size data based on fire_ids, to get steps_bounds
size_data <- size_data[fire_ids, ]
steps_bounds <- cbind(rep(5, J1), size_data$steps_upper)
rownames(steps_bounds) <- fire_ids

# GET FWIZ range in the study period --------------------------------------

pnnh <- vect(file.path("data", "protected_areas", "apn_limites.shp"))
pnnh <- pnnh[pnnh$nombre == "Nahuel Huapi", ]
pnnh <- project(pnnh, "EPSG:5343")

# FWI
fwi_fort <- rast(file.path("data", "fwi_daily_1998-2022", "24km", # change if using 52 km
                           "fwi_fortnights_19970701_20230630_standardized.tif"))
fwi_fort <- project(fwi_fort, "EPSG:5343")
fwi_fort <- crop(fwi_fort, pnnh, snap = "out")
fwiz_range_ori <- range(values(fwi_fort))

fwiz_range <- (fwiz_range_ori - fwi_mean) / fwi_sd
fwiz_range 
# -1.786781  3.401912 // 52 km
# -1.963114  3.671563 // 24 km

# Steps model with Stan (get inits) ---------------------------------------

# Ytry contains log-steps samples
steps_fitted1 <- as.numeric(apply(Ytry["steps", , ], 1, function(x) mean(exp(x))))

sdata1 <- list(
  N1 = J1, N2 = J2,
  
  steps1 = steps_fitted1,
  fwi1 = fwi1,
  fwi2 = fwi2,
  
  areaL = areaL,
  
  area1 = area1,
  area2 = area2,
  
  L = 2, # forces to simulate spread outside the ignition
  Umax = 60000 / 30, # 60 km max
  Umin = ceiling(max(steps_fitted1)) + 0.001,
  
  prior_steps_int_mn = 0,
  prior_steps_int_sd = 100,
  prior_steps_b_sd = 10,
  prior_steps_sigma_sd = 2,
  
  prior_area_int_mn = mean(c(area1, area2)),
  prior_area_int_sd = sd(c(area1, area2)),
  prior_area_b_sd = 3,
  prior_area_sigma_sd = sd(c(area1, area2)) * 2
)
L <- sdata1$L


smodel <- stan_model("spread/steps_model_logitnorm_ls-fixed.stan")

stansteps <- sampling(
  smodel, sdata1, seed = 234669, refresh = 100,
  cores = 4, chains = 4, iter = 4000, warmup = 1000
)

mcmc_dens(stansteps, pars = c("steps_int", "steps_b", "U"),
          facet_args = list(nrow = 1))

mcmc_pairs(stansteps, pars = c("steps_int", "steps_b", "U", "steps_sigma"))

# Plot steps curves
spar <- as.matrix(stansteps, pars = c("steps_int", "steps_b", "U",
                               "steps_sigma"))

ni <- 500
ids <- sample(1:nrow(spar), ni, F)
nr <- 200
fwiseq <- seq(min(c(fwi1, fwi2)), 10, length.out = nr)
mumat_logit <- matrix(NA, nr, ni)
mumat <- matrix(NA, nr, ni)
ylims <- matrix(NA, nr, 2)

for(i in 1:ni) {
  it = ids[i]
  mumat_logit[, i] <- spar[it, "steps_int"] + spar[it, "steps_b"] * fwiseq
  for(j in 1:nr) {
    mumat[j, i] <- 
      momentsLogitnorm(mumat_logit[j, i], spar[it, "steps_sigma"])["mean"] *
      (spar[it, "U"] - L) + L
  }
}

# limits of predictive distribution
for(r in 1:nr) {
  # r = 1
  yl <- rnorm(ni * 4, mumat_logit[r, ], spar[ids, "steps_sigma"])
  yy <- plogis(yl) * (spar[ids, "U"] - L) + L
  ylims[r, ] <- quantile(yy, probs = c(0.025, 0.975), method = 8)
}

pred1 <- data.frame(
  fwi = fwiseq,
  mu = rowMeans(mumat),
  mu_lower = apply(mumat, 1, quantile, probs = 0.025, method = 8),
  mu_upper = apply(mumat, 1, quantile, probs = 0.975, method = 8),
  y_lower = ylims[, 1],
  y_upper = ylims[, 2]
)

datadf1 <- data.frame(
  steps = c(
    steps_fitted1, 
    as.matrix(stansteps, "steps2") |> colMeans()
  ),
  fwi = c(fwi1, fwi2),
  type = rep(c("ABC-area-fit", "only-area-fit"), c(J1, J2))
)

p1 <- ggplot(pred1, aes(fwi, mu, ymin = mu_lower, ymax = mu_upper)) +
  geom_ribbon(mapping = aes(fwi, mu, ymin = y_lower, ymax = y_upper),
              alpha = 0.15) +
  geom_ribbon(alpha = 0.3) + 
  geom_line() + 
  geom_point(data = datadf1, mapping = aes(x = fwi, y = steps, color = type),
             inherit.aes = F, shape = 21, size = 2, stroke = 0.65) + 
  scale_color_viridis(discrete = T, end = 0.5) +
  geom_vline(xintercept = max(fwiz_range), 
             color = viridis(1, option = "C", begin = 0.5)) +
  scale_y_continuous(limits = c(0, 2000), expand = c(0.001, 00.01)) +
  theme(panel.grid.minor = element_blank()) + 
  ylab("Steps") + 
  xlab("FWI anomaly")
p1

ggsave("spread/figures_FWIZ2/step_model_stan_logitnorm.png",
       width = 15, height = 12, units = "cm", plot = p1)


# check area ~ steps:
datadf1$steps_log <- log(datadf1$steps)
datadf1$area_log <- area_all
plot(area_log ~ steps_log, datadf1)
mm <- lm(area_log ~ steps_log, data = datadf1)
abline(coef(mm), col = 2, lwd = 2)


# Steps model with Stan (fixed upper) -------------------------------------

# Ytry contains log-steps samples
steps_fitted1 <- as.numeric(apply(Ytry["steps", , ], 1, function(x) mean(exp(x))))

sdata1a <- list(
  N1 = J1, N2 = J2,
  
  steps1 = steps_fitted1,
  fwi1 = fwi1,
  fwi2 = fwi2,
  
  area1 = area1,
  area2 = area2,
  
  L = 2, # forces to simulate spread outside the ignition
  U = 2000, # 60000 / 30, # 60 km max
  
  prior_steps_int_sd = 100,
  prior_steps_b_sd = 10,
  prior_steps_sigma_sd = 2,
  
  prior_area_int_mn = mean(c(area1, area2)),
  prior_area_int_sd = sd(c(area1, area2)),
  prior_area_b_sd = 3,
  prior_area_sigma_sd = sd(c(area1, area2)) * 2
)
L <- sdata1a$L
U <- sdata1a$U

smodel1a <- stan_model("spread/steps_model_logitnorm_fixedU.stan")

stansteps1a <- sampling(
  smodel1a, sdata1a, seed = 234669, refresh = 50,
  cores = 4, chains = 4, iter = 4000, warmup = 1000
)

mcmc_dens(stansteps1a, pars = c("steps_int", "steps_b", "steps_sigma"),
          facet_args = list(nrow = 1))

mcmc_pairs(stansteps1a, pars = c("steps_int", "steps_b", "steps_sigma"))

# Plot steps curves
spar <- as.matrix(stansteps1a, pars = c("steps_int", "steps_b", "steps_sigma"))

ni <- 500
ids <- sample(1:nrow(spar), ni, F)
nr <- 200
fwiseq <- seq(min(c(fwi1, fwi2)), 10, length.out = nr)
mumat_logit <- matrix(NA, nr, ni)
mumat <- matrix(NA, nr, ni)
ylims <- matrix(NA, nr, 2)

for(i in 1:ni) {
  it = ids[i]
  mumat_logit[, i] <- spar[it, "steps_int"] + spar[it, "steps_b"] * fwiseq
  for(j in 1:nr) {
    mumat[j, i] <- 
      momentsLogitnorm(mumat_logit[j, i], spar[it, "steps_sigma"])["mean"] *
      (U - L) + L
  }
}

# limits of predictive distribution
for(r in 1:nr) {
  # r = 1
  yl <- rnorm(ni * 4, mumat_logit[r, ], spar[ids, "steps_sigma"])
  yy <- plogis(yl) * (U - L) + L
  ylims[r, ] <- quantile(yy, probs = c(0.025, 0.975), method = 8)
}

pred1 <- data.frame(
  fwi = fwiseq,
  mu = rowMeans(mumat),
  mu_lower = apply(mumat, 1, quantile, probs = 0.025, method = 8),
  mu_upper = apply(mumat, 1, quantile, probs = 0.975, method = 8),
  y_lower = ylims[, 1],
  y_upper = ylims[, 2]
)

datadf1 <- data.frame(
  steps = c(
    steps_fitted1, 
    as.matrix(stansteps1a, "steps2") |> colMeans()
  ),
  fwi = c(fwi1, fwi2),
  type = rep(c("ABC-area-fit", "only-area-fit"), c(J1, J2))
)

p1a <- ggplot(pred1, aes(fwi, mu, ymin = mu_lower, ymax = mu_upper)) +
  geom_ribbon(mapping = aes(fwi, mu, ymin = y_lower, ymax = y_upper),
              alpha = 0.15) +
  geom_ribbon(alpha = 0.3) + 
  geom_line() + 
  geom_point(data = datadf1, mapping = aes(x = fwi, y = steps, color = type),
             inherit.aes = F, shape = 21, size = 2, stroke = 0.65) + 
  scale_color_viridis(discrete = T, end = 0.5) +
  geom_vline(xintercept = max(fwiz_range), 
             color = viridis(1, option = "C", begin = 0.5)) +
  scale_y_continuous(limits = c(0, 2000), expand = c(0.001, 00.01)) +
  theme(panel.grid.minor = element_blank()) + 
  ylab("Steps") + 
  xlab("FWI anomaly")
p1a

ggsave("spread/figures_FWIZ2/step_model_stan_logitnorm_fixedU.png",
       width = 15, height = 12, units = "cm", plot = p1)


ggarrange2(p1 + theme(legend.position = "none") + ggtitle("U unknown"), 
           p1a + ggtitle("U fixed"), 
           nrow = 1)

# Steps model 1b (exploratory) with Stan -------------------------------------

# fwiz_range

steps_fitted1 <- as.numeric(apply(Ytry["steps", , ], 1, function(x) mean(exp(x))))

sdata1 <- list(
  N1 = J1, N2 = J2,
  
  steps1 = steps_fitted1,
  fwi1 = fwi1,
  fwi2 = fwi2,
  
  area1 = area1,
  area2 = area2,
  
  steps_lower = 2,
  
  upper_max = 60000 / 30, # 60 km max
  upper_min = ceiling(max(steps_fitted1)) + 0.001,
  
  prior_steps_mid_mn = mean(fwi1),
  prior_steps_mid_sd = 20, #sd(fwi1),
  prior_steps_b_sd = 10,
  prior_steps_sigma_sd = 2,
  
  prior_area_int_mn = mean(c(area1, area2)),
  prior_area_int_sd = sd(c(area1, area2)),
  prior_area_b_sd = 3,
  prior_area_sigma_sd = sd(c(area1, area2)) * 2
)


smodel1 <- stan_model("spread/steps_model_logitnorm2.stan")

m1 <- sampling(
  smodel1, sdata1, seed = 234669, refresh = 50,
  #cores = 1, chains = 1, iter = 5
  cores = 4, chains = 4, iter = 2000
)

mcmc_dens(m1, pars = c("steps_mid", "steps_b", "steps_upper"),
          facet_args = list(nrow = 1))

mcmc_pairs(m1, pars = c("steps_mid", "steps_b", "steps_upper", "steps_sigma"))

# Plot steps curves
spar <- as.matrix(m1, pars = c("steps_mid", "steps_b", "steps_upper",
                               "steps_sigma"))

ni <- 500
ids <- sample(1:nrow(spar), ni, F)
nr <- 200
fwiseq <- seq(min(c(fwi1, fwi2)), 10, length.out = nr)
mumat_logit <- matrix(NA, nr, ni)
mumat <- matrix(NA, nr, ni)
ylims <- matrix(NA, nr, 2)

for(i in 1:ni) {
  it = ids[i]
  mumat_logit[, i] <- spar[it, "steps_b"] * (fwiseq - spar[it, "steps_mid"])
  for(j in 1:nr) {
    mumat[j, i] <- 
      momentsLogitnorm(mumat_logit[j, i], spar[it, "steps_sigma"])["mean"] *
      (spar[it, "steps_upper"] - 1) + 1
  }
}

# limits of predictive distribution
for(r in 1:nr) {
  # r = 1
  yl <- rnorm(ni * 4, mumat_logit[r, ], spar[ids, "steps_sigma"])
  yy <- plogis(yl) * (spar[ids, "steps_upper"] - 1) + 1
  ylims[r, ] <- quantile(yy, probs = c(0.025, 0.975), method = 8)
}

pred1 <- data.frame(
  fwi = fwiseq,
  mu = rowMeans(mumat),
  mu_lower = apply(mumat, 1, quantile, probs = 0.025, method = 8),
  mu_upper = apply(mumat, 1, quantile, probs = 0.975, method = 8),
  y_lower = ylims[, 1],
  y_upper = ylims[, 2]
)

datadf1 <- data.frame(
  steps = c(
    steps_fitted1, 
    as.matrix(m1, "steps2") |> colMeans()
  ),
  fwi = c(fwi1, fwi2),
  type = rep(c("ABC-area-fit", "only-area-fit"), c(J1, J2))
)

p1 <- ggplot(pred1, aes(fwi, mu, ymin = mu_lower, ymax = mu_upper)) +
  geom_ribbon(mapping = aes(fwi, mu, ymin = y_lower, ymax = y_upper),
              alpha = 0.15) +
  geom_ribbon(alpha = 0.3) + 
  geom_line() + 
  geom_point(data = datadf1, mapping = aes(x = fwi, y = steps, color = type),
             inherit.aes = F, shape = 21, size = 2, stroke = 0.65) + 
  scale_color_viridis(discrete = T, end = 0.5) +
  geom_vline(xintercept = max(fwiz_range), 
             color = viridis(1, option = "C", begin = 0.5)) +
  scale_y_continuous(limits = c(0, 2000), expand = c(0.001, 00.01)) +
  theme(panel.grid.minor = element_blank()) + 
  ylab("Steps") + 
  xlab("FWI anomaly")
p1

ggsave("spread/figures_FWIZ2/step_model_stan_logitnorm2_priorwide.png",
       width = 15, height = 12, units = "cm", plot = p1)


# Steps model 2 (exploratory) with Stan -------------------------------------

sdata2 <- list(
  N1 = J1, N2 = J2,
  
  steps1 = steps_fitted1,
  fwi1 = fwi1,
  fwi2 = fwi2,
  
  area1 = area1,
  area2 = area2,
  
  prior_steps_int_mn = mean(log(steps_fitted1)),
  prior_steps_int_sd = sd(log(steps_fitted1)),
  prior_steps_b_sd = 5,
  prior_steps_sigma_sd = sd(log(steps_fitted1)) * 2,
  
  prior_area_int_mn = mean(c(area1, area2)),
  prior_area_int_sd = sd(c(area1, area2)),
  prior_area_b_sd = 3,
  prior_area_sigma_sd = sd(c(area1, area2)) * 2
)


smodel2 <- stan_model("spread/steps_model_lognorm2.stan")

m2 <- sampling(
  smodel2, sdata2, seed = 234669, refresh = 50,
  #cores = 1, chains = 1, iter = 5
  cores = 4, chains = 4, iter = 2000
)

# Plot steps curves
spar2 <- as.matrix(m2, pars = c("steps_int", "steps_b", "steps_sigma"))

ni <- 500
ids <- sample(1:nrow(spar2), ni, F)
nr <- 200
fwiseq <- seq(min(c(fwi1, fwi2)), 10, length.out = nr)
mumat_log <- matrix(NA, nr, ni)
mumat2 <- matrix(NA, nr, ni)
ylims <- matrix(NA, nr, 2)

for(i in 1:ni) {
  it = ids[i]
  mumat_log[, i] <- spar2[it, "steps_int"] + spar2[it, "steps_b"] * fwiseq
  mumat2[, i] <- exp(mumat_log[, i] + 0.5 * spar2[it, "steps_sigma"] ^ 2)
}

# limits of predictive distribution
for(r in 1:nr) {
  # r = 1
  yy <- rlnorm(ni * 4, mumat_log[r, ], spar2[ids, "steps_sigma"])
  ylims[r, ] <- quantile(yy, probs = c(0.025, 0.975), method = 8)
}

pred2 <- data.frame(
  fwi = fwiseq,
  mu = rowMeans(mumat2),
  mu_lower = apply(mumat2, 1, quantile, probs = 0.025, method = 8),
  mu_upper = apply(mumat2, 1, quantile, probs = 0.975, method = 8),
  y_lower = ylims[, 1],
  y_upper = ylims[, 2]
)

datadf2 <- data.frame(
  steps = c(
    steps_fitted1, 
    exp(as.matrix(m2, "steps_log2")) |> colMeans()
  ),
  fwi = c(fwi1, fwi2),
  type = rep(c("ABC-area-fit", "only-area-fit"), c(J1, J2))
)

p2 <- ggplot(pred2, aes(fwi, mu, ymin = mu_lower, ymax = mu_upper)) +
  geom_ribbon(mapping = aes(fwi, mu, ymin = y_lower, ymax = y_upper),
              alpha = 0.15) +
  geom_ribbon(alpha = 0.3) + 
  geom_line() + 
  geom_point(data = datadf2, mapping = aes(x = fwi, y = steps, color = type),
             inherit.aes = F, shape = 21, size = 2, stroke = 0.65) + 
  scale_color_viridis(discrete = T, end = 0.5) +
  geom_vline(xintercept = max(fwiz_range), 
             color = viridis(1, option = "C", begin = 0.5)) +
  scale_y_continuous(expand = c(0.001, 00.01)) +
  theme(panel.grid.minor = element_blank()) + 
  ylab("Steps") + 
  xlab("FWI anomaly")

ggsave("spread/figures_FWIZ2/step_model_stan_lognorm_raw.png",
       width = 15, height = 12, units = "cm")

pboth <- deeptime::ggarrange2(
  p1 + 
    scale_y_continuous(limits = c(0, 2000), expand = c(0.001, 00.01)) +
    theme(legend.position = "none"), 
  p2 + 
    scale_y_continuous(limits = c(0, 2000), expand = c(0.001, 00.01)) +
    theme(legend.position = "bottom"), 
  nrow = 1
)

ggsave("spread/figures_FWIZ2/step_models_stan_logit-log-norm_priorwide.png",
       width = 18, height = 11, units = "cm",
       plot = pboth)


p3 <- p2 + 
  scale_y_continuous(limits = c(0, 5000), expand = c(0.001, 00.01))
ggsave("spread/figures_FWIZ2/step_model_stan_lognorm_raw2.png",
       width = 18, height = 11, units = "cm",
       plot = p3)


# MLE estimates -----------------------------------------------------------

# using samples from the first-step abc, get a point estimate of all parameters.
par_start <- vector("list", 3)
names(par_start) <- c("fixef", "ranef", "steps")
nsim <- 5000
set.seed(123)
ii <- sample(1:N1, size = nsim)

par_start[["fixef"]] <- array(NA, dim = c(n_coef + 1, 3, nsim),
                    dimnames = list(
                      par_names = c(par_names, "area"),
                      par_class = c("a", "b", "s2"),
                      iter = 1:nsim
                    ))

par_start[["ranef"]] <- Ytry[, , 1:nsim]
par_start[["ranef"]][, , ] <- NA

mm <- matrix(NA, J2, nsim)
rownames(mm) <- fire_ids_nonspread
par_start[["steps"]] <- mm

# # Compute estimates
# for(i in 1:nsim) {
#   # i = 1
#   print(i)
#
#   ranef_ <- Ytry[, , ii[i]]
#   par_start[["ranef"]][, , i] <- ranef_
#
#   # simple parameters (before steps)
#   for(v in 1:(n_coef-1)) {
#     parvals <- ranef_[v, ]
#     mod <- lm(parvals ~ fwi1)
#     par_start[["fixef"]][v, , i] <- c(coef(mod), sigma(mod) ^ 2)
#   }
#
#   # steps ~ area model
#   # used to provide initial values for the non-spread steps parameters
#   ss <- ranef_["steps", ]
#   dd <- data.frame(ss = ss, aa = area1)
#   mod <- lm(ss ~ aa, data = dd)
#   # simulate steps in non-spread fires (careful: do not predict, simulate)
#   mu <- predict(mod, data.frame(aa = area2))
#   sigma <- sigma(mod)
#   steps_sim <- rnorm(J2, mu, sigma)
#   par_start[["steps"]][, i] <- steps_sim
#
#   # estimate steps ~ fwi and area ~ steps using all fires
#   ss_full <- c(ss, steps_sim)
#   # steps ~ fwi
#   mod <- lm(ss_full ~ fwi_all)
#   par_start[["fixef"]]["steps", , i] <- c(coef(mod), sigma(mod) ^ 2)
#   # area ~ steps (truncated normal regression)
#   mod <- truncreg(area_all ~ ss_full, point = areaL, direction = "left")
#   par_start[["fixef"]]["area", , i] <- c(coef(mod)[1:2], coef(mod)[3] ^ 2)
# }
# saveRDS(par_start, file.path("files", "hierarchical_model_FWIZ", "par_start.rds"))
par_start <- readRDS(file.path("files", "hierarchical_model_FWIZ", "par_start.rds"))

# Priors for hyperparameters ----------------------------------------------

# priors for intercepts and slopes. b0 has the means, S0 has the vcov (diagonal).
b0 <- matrix(0, 2, n_coef + 1)
colnames(b0) <- par_names_all
rownames(b0) <- c("a", "b")
b0["a", "steps"] <- mean(par_start$fixef["steps", "a", ])
b0["a", "area"] <- mean(par_start$fixef["area", "a", ])

# The prior sd for regression coefficients is the same for all parameters
# (sd = 10)
S0 <- diag(rep(10 ^ 2, 2)) # used for truncnorm regression (area), updated with MH
S0_inv <- solve(S0)        # used for gibbs updates

# s2 for all parameters has inv-gamma prior, for conjugacy.
# invgamma::dinvgamma(x, t0, d0)
t0 <- 1; d0 <- 1 / 1000

# priors used in stage1, to be subtracted. Apply at the logit scale, MVN.
priors1 <- readRDS(file.path(target_dir, "priors_stage1.rds"))
mu0 <- priors1$mu0[, fire_ids]
Sigma0 <- priors1$S0[, , fire_ids]

# MCMC adaptation ---------------------------------------------------------

# run a long chain to get good starting points
nsim <- 10000

system.time(
  run0 <- mcmc(nsim = nsim, sd_jump_out = T)
)
# 50.113 s / 1000 iter
# 50.113 * 1000 / 3600 # 13 h para correr 1e6 muestras. OK.
xiter <- 1:nsim

# 489.122 / 10000 # 0.0489122 / iter ## cuando fueron 10000, sin print
# 1e6 * 0.0489122 / 3600 # sí, 13.59 h

par(mfrow = c(2, 4))
for(v in 1:(n_coef+1)) {
  # v = 1
  yy <- range(c(run0$fixef[v, 1:2, ], sqrt(run0$fixef[v, "s2", ])))
  plot(run0$fixef[v, "a", ] ~ xiter, type = "l", ylim = yy,
       main = par_names_all[v])
  lines(run0$fixef[v, "b", ] ~ xiter, col = "red")
  lines(sqrt(run0$fixef[v, "s2", ]) ~ xiter, col = "green")
}
par(mfrow = c(1, 1))

saveRDS(run0, file.path("files", "hierarchical_model_FWIZ", "run0.rds"))
run0 <- readRDS(file.path("files", "hierarchical_model_FWIZ", "run0.rds"))

# Initial proposal sigma:
tmp <- run0$fixef["area", , ]
tmp[3, ] <- sqrt(tmp[3, ]) # sd_jump for sigma, not sigma2
area_jump <- apply(tmp, 1, sd)

sd_jump1 <- list(
  area = sqrt(area_jump ^ 2 * 2),
  steps = sqrt(apply(run0$steps, 1, sd) ^ 2 * 2)
)
# posterior sd after a long run, just for curiosity
sd_post <- list(
  area = apply(run0$fixef[7, , 1000:nsim, drop = F], 1:2, sd),
  steps = apply(run0$steps[, 1000:nsim], 1, sd)
)

ns <- 1000 # iter by step

# first run, to get acceptance
run1 <- mcmc(nsim = ns, sd_jump = sd_jump1, samples = run0, sd_jump_out = T)
a1 <- acceptance(run1)

# acceptance trial steps
K <- 20

# make space for acceptance and sigma
accept_track <- list(
  area = run1$fixef[n_coef+1, , 1:K],
  steps = run1$steps[, 1:K]
)
sigma_track <- list(
  area = run1$fixef[n_coef+1, , 1:K],
  steps = run1$steps[, 1:K]
)
accept_track$steps[, ] <- NA
sigma_track$steps[, ] <- NA
# just make space

# fill sigma and accept in the first run
accept_track$area[, 1] <- a1$area
accept_track$steps[, 1] <- a1$steps
sigma_track$area[, 1] <- sd_jump1$area
sigma_track$steps[, 1] <- sd_jump1$steps

# Fill 4 sigma values below and above the first, to fit the first regression
# using 5 data points
factors <- seq(0.1, 2, by = 0.2)
lf <- length(factors)
for(k in 2:(lf+1)) {
  print(k)
  sigma_track$area[, k] <- sd_jump1$area * factors[k-1]
  sigma_track$steps[, k] <- sd_jump1$steps * factors[k-1]

  # run MCMC
  sss <- list(area = sigma_track$area[, k],
              steps = sigma_track$steps[, k])
  run_k <- mcmc(ns, sd_jump = sss, samples = run0)

  # compute and store acceptance
  a_k <- acceptance(run_k)
  accept_track$area[, k] <- a_k$area
  accept_track$steps[, k] <- a_k$steps
}

# iterative fitting
for(k in (lf+2):K) {
  # k = 12
  print(k)
  # use previous runs to fit a regression of sigma ~ accept, and choose the
  # predicted sigma for accept = 0.44

  # area parameters
  for(j in 1:3) { # loop over (a, b, s2)
    # j = 1
    dd <- data.frame(
      aa = accept_track$area[j, 1:(k-1)],
      ss = sigma_track$area[j, 1:(k-1)]
    )

    # remove sigma too close to zero
    dd <- dd[dd$ss >= 1e-4, ]
    dd$lss = log(dd$ss)

    # fit regression
    if(k < 10) {
      mm <- lm(lss ~ aa + I(aa ^ 2), data = dd)
    } else {
      mm <- gam(lss ~ s(aa, k = 6), data = dd, method = "REML")
    }

    ss_pred <- predict(mm, newdata = data.frame(aa = 0.44),
                       se.fit = F) |> exp()
    sigma_track$area[j, k] <- ifelse(ss_pred < 1e-4, 1e-4, ss_pred)
  }

  # steps
  for(j in 1:J2) {
    aa <- accept_track$steps[j, 1:(k-1)]
    lss <- log(sigma_track$steps[j, 1:(k-1)])

    # get previous data
    dd <- data.frame(
      aa = accept_track$steps[j, 1:(k-1)],
      ss = sigma_track$steps[j, 1:(k-1)]
    )

    # remove sigma too close to zero
    dd <- dd[dd$ss >= 1e-4, ]
    dd$lss = log(dd$ss)

    # fit regression
    if(k < 10) {
      mm <- lm(lss ~ aa + I(aa ^ 2), data = dd)
    } else {
      mm <- gam(lss ~ s(aa, k = 6), data = dd, method = "REML")
    }

    ss_pred <- predict(mm, newdata = data.frame(aa = 0.44),
                       se.fit = F) |> exp()
    sigma_track$steps[j, k] <- ifelse(ss_pred < 1e-4, 1e-4, ss_pred)
  }

  # run MCMC
  sss <- list(area = sigma_track$area[, k],
              steps = sigma_track$steps[, k])
  run_k <- mcmc(ns, sd_jump = sss, samples = run0)

  # compute and store acceptance
  a_k <- acceptance(run_k)
  accept_track$area[, k] <- a_k$area
  accept_track$steps[, k] <- a_k$steps
}

# Fit reverse model (acceptance ~ sigma) and choose that value.
sd_jump_tune <- list(
  area = sigma_track$area[, K],  # place-holder
  steps = sigma_track$steps[, K]
)

# area parameters
for(j in 1:3) {
  # get data
  dd <- data.frame(
    aa = accept_track$area[j, 1:K],
    ss = sigma_track$area[j, 1:K]
  )

  # fit gam
  mm <- gam(aa ~ s(ss, k = 6), family = betar(),
            data = dd, method = "REML")

  fn <- function(ss) {
    apred <- predict(mm, newdata = data.frame(ss = ss), type = "response")
    return((apred - 0.44) ^ 2)
  }

  opt <- optim(sigma_track$area[j, K], fn, method = "Brent",
               lower = min(dd$ss), upper = max(dd$ss))

  sd_jump_tune$area[j] <- opt$par

  # Visualize
  tit <- paste("area", colnames(sd_jump1$area)[j], sep = "; ")
  ppp <- data.frame(ss = seq(min(dd$ss), max(dd$ss), length.out = 100))
  ppp$y <- predict(mm, ppp, type = "response")
  plot(aa ~ ss, data = dd, main = tit, xlab = "Sigma", ylab = "Acceptance")
  lines(y ~ ss, data = ppp)
  abline(v = opt$par, col = 2, lty = 2)
  abline(h = 0.44, col = 4, lty = 2)
}

# steps
for(j in 1:J2) {
  dd <- data.frame(
    aa = accept_track$steps[j, 1:K],
    ss = sigma_track$steps[j, 1:K]
  )

  # fit gam
  mm <- gam(aa ~ s(ss, k = 6), family = betar(),
            data = dd, method = "REML")

  fn <- function(ss) {
    apred <- predict(mm, newdata = data.frame(ss = ss), type = "response")
    return((apred - 0.44) ^ 2)
  }

  opt <- optim(sigma_track$steps[j, K], fn, method = "Brent",
               lower = min(dd$ss), upper = max(dd$ss))

  sd_jump_tune$steps[j] <- opt$par
}

# compare posterior sd with tune sd
plot(sd_jump_tune$steps ~ sd_post$steps)
mm <- lm(sd_jump_tune$steps ~ sd_post$steps - 1)
abline(c(0, coef(mm))) # coef(mm) = 2.3

fftune <- as.vector(sd_jump_tune$area)
ffpost <- as.vector(sd_post$area)
plot(fftune ~ ffpost)
mm <- lm(fftune ~ ffpost - 1)
abline(c(0, coef(mm))) # coef(mm) = 0.28

saveRDS(sd_jump_tune, file.path("files", "hierarchical_model_FWIZ", "sd_jump_tune.rds"))

# Run MCMC in parallel ----------------------------------------------------

run0 <- readRDS(file.path("files", "hierarchical_model_FWIZ", "run0.rds"))
sd_jump_tune <- readRDS(file.path("files", "hierarchical_model_FWIZ", "sd_jump_tune.rds"))

n_cores <- 8 # takes the same time as with more cores
nsave <- 150
thin <- 1000
nb <- 10     # it takes long, run in batches

# Prepare init samples
eend <- dim(run0$fixef)[3]
bbegin <- floor(eend * 0.7)
iii <- sample(bbegin:eend, size = n_cores, replace = F)
start_samples8 <- list(
  fixef = run0$fixef[, , iii],
  ranef = run0$ranef[, , iii],
  steps = run0$steps[, iii]
)

Sys.time()
for(batch in 1:nb) {
  print(batch)

  if(batch == 1) {
    init_samples <- start_samples8
  } else {
    # load previous batch
    nump <- stringr::str_pad(batch-1, 2, side = "left", pad = "0")
    fname_prev <- paste("draws_batch_", nump, ".rds", sep = "")
    draws_prev <- readRDS(file.path("files", "hierarchical_model_FWIZ", fname_prev))
    # get last iters
    last <- dim(draws_prev$fixef)[3]
    init_samples <- list(
      fixef = draws_prev$fixef[, , last, ],
      ranef = draws_prev$ranef[, , last, ],
      steps = draws_prev$steps[, last, ]
    )
  }

  draws <- mcmc_parallel(nsim = nsave * thin, thin = thin, n_cores = n_cores,
                         start_samples = init_samples, sd_jump = sd_jump_tune)
  num <- stringr::str_pad(batch, 2, side = "left", pad = "0")
  fname <- paste("draws_batch_", num, ".rds", sep = "")
  saveRDS(draws, file.path("files", "hierarchical_model_FWIZ", fname))
}
Sys.time()
# start: 2024-11-07 04:44:55 -03
# end:   2024-11-08 05:31:00 -03

# Tidy samples ------------------------------------------------------------

ff <- list.files(file.path("files", "hierarchical_model_FWIZ"), pattern = "draws_batch_")
nsave <- 150

# assign growing iter_id to samples from consecutive batches
dlist <- lapply(1:length(ff), function(i) {
  x <- readRDS(file.path("files", "hierarchical_model_FWIZ", ff[i]))
  # assign NOT growing chain_id
  iter_names <- as.character((1:nsave) + (i-1) * nsave)
  print(i)
  nitems <- length(x)
  for(j in 1:nitems) {
    dimnames(x[[j]])[["iter"]] <- iter_names
  }
  return(x)
})


# Merge batches in a list with arrays (draws), which will be exported.
# In addition, create draws_arrays to summarize with posterior package.
partypes <- c("fixef", "rho", "ranef", "steps")
nitems <- length(partypes)
draws <- vector("list", nitems)
names(draws) <- partypes
dimpar <- c(2, 2, 2, 1)
# dimension of a single posterior sample for each parameter type
draws_arrays <- vector("list", nitems)
names(draws_arrays) <- partypes
nc <- 8
ni <- 1500
npost <- nc * ni

for(p in 1:nitems) {
  draws_temp <- abind(lapply(dlist, function(x) x[[p]]),
                      along = dimpar[p] + 1) # along iter dimension

  if(dimpar[p] > 1) {
    nvar <- prod(dim(draws_temp)[1:dimpar[p]]) # number of (collapsed) parameters
    drarr <- array(draws_temp, c(nvar, ni, nc))
    # this collapses the first two dimensions, which in most cases have parameters
    # in matrix form

    # # check ordering:
    # draws_temp[, , 2, 1, drop = F] # third col, first row
    # drarr[1:nvar, 2, 1, drop = F] # element 6 * 2 + 1

    # set name to collapsed variable
    vardf <- expand.grid(
      v1 = dimnames(draws_temp)[[1]],
      v2 = dimnames(draws_temp)[[2]]
    )

    dimnames(drarr) <- list(
      variable = paste(vardf$v1, vardf$v2, sep = "__"),
      iteration = dimnames(draws_temp)[["iteration"]],
      chain = dimnames(draws_temp)[["chain"]]
    )
  } else {
    drarr <- draws_temp
    names(dimnames(drarr)) <- c("variable", "iteration", "chain")
  }

  # remove redundant parameters
  if(partypes[p] == "rho") {
    rows_keep <- which(lower.tri(matrix(1:n_coef^2, n_coef, n_coef)))
    drarr <- drarr[rows_keep, , ]
  }

  draws_arrays[[p]] <- aperm(drarr, c(2, 3, 1)) # get iteration, chain, variable.

  # collapse iterations and chains, in a single "iterations" dimension
  par_dims <- 1:dimpar[p]
  iter_dims <- (1:2) + max(par_dims)

  # permute to collapse (dims to collapse must be the first ones)
  perm <- c(iter_dims, par_dims)
  draws2 <- aperm(draws_temp, perm)
  draws3 <- array(draws2, c(ni*nc, dim(draws_temp)[par_dims]))
  new_par_dims <- par_dims + 1
  draws4 <- aperm(draws3, c(new_par_dims, 1))

  par_dn <- dimnames(draws_temp)[par_dims]
  dn <- c(par_dn, list(1:(ni*nc)))
  names(dn) <- c(rep("variable", dimpar[p]), "iteration")
  dimnames(draws4) <- dn

  draws[[p]] <- draws4
}

saveRDS(draws, file.path("files", "hierarchical_model_FWIZ",
                         "spread_model_samples.rds"))

# summaries
summaries <- lapply(draws_arrays, function(x) summarise_draws(x))
for(i in 1:nitems) {
  print(partypes[i])
  mm <- apply(summaries[[i]][, c("ess_tail", "ess_bulk", "rhat")],
              2, range, na.rm = T)
  print(mm)
}
# [1] "fixef"
#      ess_tail  ess_bulk      rhat
# [1,]  9283.017  5394.024 0.9997468
# [2,] 12113.424 12154.234 1.0009871

# [1] "rho"
#      ess_tail  ess_bulk      rhat
# [1,] 10540.83  9121.035 0.9998474
# [2,] 12037.66 12269.805 1.0013950

# [1] "ranef"
#      ess_tail  ess_bulk      rhat
# [1,]  1084.75  1329.899 0.9996582
# [2,] 12348.54 12441.650 1.0027232

# [1] "steps"
#      ess_tail ess_bulk      rhat
# [1,] 10501.14 10940.43 0.9995957
# [2,] 12240.31 12496.06 1.0008176


# Load tidy samples -------------------------------------------------------

draws <- readRDS(file.path("files", "hierarchical_model_FWIZ",
                           "spread_model_samples.rds"))
npost <- ncol(draws$steps)

# Correlation plots -------------------------------------------------------

# Make density of correlation coefficient between parameters, conditional and
# marginal to FWI. (Conditional fixes it at the all-fires mean = 0).

combs <- as.data.frame(combn(par_names, 2) |> t())
combs$name <- paste(combs$V1, combs$V2, sep = "_")
ncomb <- nrow(combs)

dcorr_list <- vector("list", ncomb)

for(i in 1:ncomb) {
  # i = 1
  print(i)

  v1 <- combs$V1[i]
  v2 <- combs$V2[i]

  # Make Vcov
  rho <- draws$rho[v1, v2, ]
  s2_1 <- draws$fixef[v1, "s2", , drop = T]
  s2_2 <- draws$fixef[v2, "s2", , drop = T]
  covv <- sqrt(s2_1) * sqrt(s2_2) * rho

  V <- array(NA, dim = c(2, 2, npost))
  V[1, 1, ] <- s2_1
  V[2, 2, ] <- s2_2
  V[1, 2, ] <- V[2, 1, ] <- covv

  # Subset fixed effects
  coef1 <- draws$fixef[v1, -3, ] # get a,b, remove s2
  coef2 <- draws$fixef[v2, -3, ]

  # get mu at unconstrained scale
  mu1 <- Xlong %*% coef1
  mu2 <- Xlong %*% coef2
  mus <- abind(mu1, mu2, along = 3)
  mus <- aperm(mus, c(1, 3, 2))

  # Simulate parameters at unconstrained scale
  sims_marg <- sapply(1:npost,
    function(j) {
      rmvn(J, mus[, , j], V[, , j])
    },
    simplify = "array"
  )

  zeroes <- matrix(0, J, 2)
  sims_cond <- sapply(1:npost,
    function(j) {
      rmvn(J, zeroes, V[, , j])
    },
    simplify = "array"
  )

  # Scale to unconstrained space
  sims_marg2 <- sims_marg
  sims_cond2 <- sims_cond

  vv <- c(v1, v2)
  for(p in 1:2) {
    sims_marg2[, p, ] <- constrain_vec(sims_marg[, p, ], vv[p], support)
    sims_cond2[, p, ] <- constrain_vec(sims_cond[, p, ], vv[p], support)
  }

  # Compute correlation coefficients
  corr_marg <- apply(sims_marg2, 3, function(x) cor(x)[1, 2])
  corr_cond <- apply(sims_cond2, 3, function(x) cor(x)[1, 2])

  # Compute density
  dmarg <- density(corr_marg, from = -1, to = 1, n = 2 ^ 10, adjust = 1.5)
  dcond <- density(corr_cond, from = -1, to = 1, n = 2 ^ 10, adjust = 1.5)

  dcor <- data.frame(
    dens = c(dmarg$y, dcond$y),
    x = c(dmarg$x, dcond$x),
    type = factor(rep(c("Marginal to FWI", "Conditional to FWI"),
                      each = length(dmarg$y)),
                  levels = c("Marginal to FWI", "Conditional to FWI")),
    V1 = v1,
    V2 = v2,
    name = combs$name[i]
  )
  dcorr_list[[i]] <- dcor
}

# Make a list of plots where the missing combinations are just blank

# nested list to include blank plots
plist <- vector("list", 5)
for(p in 1:5) plist[[p]] <- vector("list", 5)

# layout matrix
lay <- matrix(
  c(
    1, 2, 3, 4, 5,
    NA, 6, 7, 8, 9,
    NA, NA, 10, 11, 12,
    NA, NA, NA, 13, 14,
    NA, NA, NA, NA, 15
  ),
  ncol = 5
)

col_labs <- par_names[-n_coef]
row_labs <- par_names[-1]

lwd_y <- 0.4
lwd_x <- 0.1
bgcol <- "#1a1a1aff"

for(col in 1:5) {
  for(row in 1:5) {
    # col = 1; row = 2
    ii <- lay[row, col]

    if(is.na(ii)) {
      data <- data.frame(V2 = row_labs[row],
                         V1 = col_labs[col])

      # blank plot
      plotcito <- ggplot(data) + geom_blank() +
        theme(panel.border = element_blank())

      if(row == 1) {
        plotcito <- plotcito +
          facet_grid(V2 ~ V1, switch = "y") +
          theme(strip.background = element_rect(color = bgcol, fill = bgcol),
                strip.text = element_text(color = "white", size = 11),
                strip.background.y = element_blank(),
                strip.text.y = element_blank())
      }

    } else {
      data <- dcorr_list[[ii]]

      plotcito <- ggplot(data) +
        geom_vline(xintercept = 0, linetype = "dotted",
                   linewidth = 0.55) +
        geom_ribbon(mapping = aes(x = x, ymin = 0, ymax = dens, fill = type),
                    color = NA, alpha = 0.4) +
        geom_line(mapping = aes(x = x, y = dens, color = type), linewidth = 0.5) +
        scale_color_viridis(option = "C", discrete = T, end = 0.5) +
        scale_fill_viridis(option = "C", discrete = T, end = 0.5) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0), limits = c(-1, 1)) +
        expand_limits(y = max(data$dens) * 1.05) +
        facet_grid(V2 ~ V1, switch = "y") +
        theme(legend.title = element_blank(),
              legend.position = "none",
              panel.grid = element_blank(),
              panel.spacing.x = unit(0,"line"),
              panel.border = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title = element_blank(),
              axis.line.x = element_line(linewidth = lwd_x),
              axis.line.y = element_line(linewidth = lwd_y),
              strip.background = element_rect(color = bgcol, fill = bgcol),
              strip.text = element_text(color = "white", size = 11),
              strip.placement = "inside")
      # plotcito

      # remove details from inner plots
      if(row > 1) {
        plotcito <- plotcito +
          theme(strip.background.x = element_blank(),
                strip.text.x = element_blank())
      }

      if(col > 1) {
        plotcito <- plotcito +
          theme(strip.background.y = element_blank(),
                strip.text.y = element_blank(),
                strip.text.x = element_text(color = "white", size = 11))
                # because ggarrange2 separates the line from the strip
      }

      if(row < 5) {
        plotcito <- plotcito +
          theme(axis.ticks.x = element_blank(),
                axis.text.x = element_blank())
      }

      if(ii == 12) {
        plotcito <- plotcito +
          theme(legend.position = "bottom",
                axis.title.x = element_text()) +
          xlab("Correlation coefficient")
      }
    }

    plist[[col]][[row]] <- plotcito
  }
}

# longanize
plist_long <- do.call("c", plist)

# plot_cor <- ggarrange2(plots = plist_long, byrow = F)
plot_cor <- egg::ggarrange(plots = plist_long, byrow = F)
ggsave("spread/figures_FWIZ/pairs_plot_correlations_posteriors.png",
       plot = plot_cor,
       width = 22, height = 20, units = "cm")



# Spread parameters as function of FWI ------------------------------------

# a panel for each parameter, at the constrained scale, as a function of FWI
# at its original scale. Use the FWI within the simulated fires?

nr <- 200 # number of random effects to simulate in order to compute mean
npred <- 150

fwi_seq <- seq(min(fwi_all), max(fwi_all), length.out = npred)
X <- cbind(rep(1, npred), fwi_seq)
fwi_seq_ori <- (fwi_seq * fwi_sd) + fwi_mean

ranef_raw <- matrix(rnorm(nr * n_coef), ncol = n_coef)

# random effects array, placeholders
ranef_tmp <- array(NA, dim = c(npred, n_coef, nr))
ranef_cons <- array(NA, dim = c(npred, n_coef, nr))

# array with posterior samples for the predicted average across raneffs
mu_samples <- array(NA, dim = c(npred, n_coef, npost))

# Loop to compute means
for(i in 1:npost) {
  if(i %% 100 == 0) print(i)

  # mu at unconstrained scale
  mumat <- X %*% t(draws$fixef[1:n_coef, c("a", "b"), i])

  # Compute choleski factor of vcov matrix for random effects
  sds <- draws$fixef[1:n_coef, "s2", i] |> sqrt()
  rho <- draws$rho[, , i]
  V <- diag(sds) %*% rho %*% diag(sds)
  Vchol_U <- chol(V)

  # unconstrained centred random effects
  ranef_centred <- ranef_raw %*% Vchol_U

  # unconstrained random effects
  for(j in 1:nr) {
    ranef_tmp[, , j] <- t(t(mumat) + ranef_centred[j, ])
  }

  # constrain (only logit-link ones, leave steps at log)
  for(v in 1:(n_coef-1)) {
    # v = 1
    ranef_cons[, v, ] <- constrain_vec(ranef_tmp[, v, ], v = par_names[v],
                                         support = support)
  }
  ranef_cons[, n_coef, ] <- ranef_tmp[, n_coef, ] # steps not modified

  # average raneffs
  mu_samples[, , i] <- apply(ranef_cons, 1:2, mean)
}
saveRDS(mu_samples, file.path("files", "hierarchical_model_FWIZ", "mu_samples_prediction.rds"))
mu_samples <- readRDS(file.path("files", "hierarchical_model_FWIZ", "mu_samples_prediction.rds"))

# summarize posterior and longanize.
str(mu_samples)
mu_summ <- apply(mu_samples, 1:2, summarise)
dimnames(mu_summ) <- list(
  metric = dimnames(mu_summ)[[1]],
  row = 1:npred,
  par_name = par_names
)
mu_summ_sub <- mu_summ[c("mean", "eti_lower_95", "eti_upper_95"), , ]
mu_df1 <- as.data.frame.table(mu_summ_sub, responseName = "par_value")
mu_df <- pivot_wider(mu_df1, names_from = "metric", values_from = "par_value")
mu_df$row <- as.numeric(as.character(mu_df$row))

# merge with fwi data
df_fwi <- data.frame(row = 1:npred, fwi_z = fwi_seq,
                     fwi = fwi_seq_ori)
mu_df <- left_join(mu_df, df_fwi, by = "row")

# Prepare random effects (points)

# Constrain
ranef_cons <- draws$ranef
for(i in 1:npost) {
  ranef_cons[, , i] <- t(constrain(t(draws$ranef[, , i]), support))
  ranef_cons[n_coef, , i] <- draws$ranef[n_coef, , i] # keep log steps
}

# summarise
ranef_summ <- apply(ranef_cons, 1:2, summarise)
names(dimnames(ranef_summ)) <- c("type", "par_name", "fire_id")
ranef_df1 <- as.data.frame.table(ranef_summ, responseName = "par_value")
ranef_df <- pivot_wider(ranef_df1, names_from = "type", values_from = "par_value")

# merge with fwi data
df_fwi_ranef <- fires_data_spread[, c("fwi", "fwi_fort_expquad")]
colnames(df_fwi_ranef) <- c("fwi_z", "fwi")
df_fwi_ranef$fire_id <- rownames(df_fwi_ranef)

ranef_df <- left_join(ranef_df, df_fwi_ranef, by = "fire_id")
ranef_df$par_name <- factor(ranef_df$par_name, levels = par_names)

# compute probabilities of FWI having a positive slope at the unconstrained scale.
b_probs <- apply(draws$fixef[1:n_coef, "b", ], 1, function(x) {
  sum(x > 0) / length(x)
})
b_probs_text <- paste(round(b_probs * 100, 2), "%")

# get higher values for lineranges to define text position
high <- aggregate(eti_upper_95 ~ par_name, ranef_df, max)[, "eti_upper_95"]
low <- aggregate(eti_lower_95 ~ par_name, ranef_df, min)[, "eti_lower_95"]
width <- high - low
ypos <- high - width * 0.15

probs_data <- data.frame(par_name = factor(par_names, levels = par_names),
                         prob = b_probs_text,
                         x = 2,
                         y = ypos)

# Plot
lwd_y <- 0.4
lwd_x <- 0.1
bgcol <- "#1a1a1aff"

vir1 <- viridis(1, option = "C")
vir2 <- viridis(1, begin = 0.5, option = "C")

ggplot(mu_df, aes(fwi, mean, ymin = eti_lower_95, ymax = eti_upper_95)) +
  geom_ribbon(color = NA, alpha = 0.4, fill = vir1) +
  geom_line(color = vir1) +

  geom_linerange(data = ranef_df, alpha = 0.5,
                 color = vir2) +
  geom_point(data = ranef_df, alpha = 0.7, shape = 21, stroke = 0.35,
             color = vir1, fill = vir2) +

  geom_text(aes(x, y, label = prob), data = probs_data, inherit.aes = F,
            size = 7/.pt) +

  facet_wrap(vars(par_name), scales = "free_y", strip.position = "left",
             axes = "all") +

  xlab("Fire Weather Index anomaly") +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  theme(strip.placement = "outside",
        strip.background = element_rect(color = "white", fill = "white"),
        strip.text = element_text(margin = margin(r = 0, l = 4, unit = "mm"),
                                  size = 11),
        panel.spacing.x = unit(0, "mm"),
        panel.spacing.y = unit(6, "mm"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(linewidth = 0.3))

ggsave("spread/figures_FWIZ/params_fwi_regression.png",
       width = 17, height = 10, units = "cm")


# Steps-FWI-Area predictions (all fires) -----------------------------------

aall <- abind::abind(list(fixef_arr, ranef_arr, steps_arr), along = 3)
draws_arr <- as_draws_array(aall)
summ <- summarise_draws(draws_arr)

steps_summ <- summ[grep("ranef__steps__", summ$variable), ]
steps_summ$fwi <- fwi_all
steps_summ$fwi_ori <- fwi_all * fwi_sd + fwi_mean
steps_summ$area_log <- area_all
steps_summ$ig_location <- c(rep("Known", J1),
                            rep("Unknown", J2))


# fwi prediction
npred <- 150
pred_fwi <- data.frame(fwi = seq(min(fwi_all), max(fwi_all), length.out = npred))
predmat <- matrix(NA, npred, npost)
for(i in 1:npred) {
  predmat[i, ] <-
    as.vector(draws_arr[, , "fixef_steps_a"]) +
    as.vector(draws_arr[, , "fixef_steps_b"]) * pred_fwi$fwi[i]
}
pred_fwi$mu <- rowMeans(predmat)
pred_fwi$mu_lower <- apply(predmat, 1, quantile, prob = 0.025)
pred_fwi$mu_upper <- apply(predmat, 1, quantile, prob = 0.975)
pred_fwi$fwi_ori <- pred_fwi$fwi * fwi_sd + fwi_mean

# add title wierdly
pred_fwi$ig <- "Ignition location"
steps_summ$ig <- "Ignition location"

fig_steps <-
ggplot(steps_summ, aes(fwi_ori, mean, ymin = q5, ymax = q95)) +

  geom_smooth(aes(fwi_ori, mean), inherit.aes = F,
              method = "lm", se = F, linetype = "dashed", color = vir1,
              linewidth = 0.35) +

  geom_ribbon(data = pred_fwi,
              mapping = aes(fwi_ori, mu, ymin = mu_lower, ymax = mu_upper),
              inherit.aes = F, color = NA, alpha = 0.4, fill = vir1) +
  geom_line(data = pred_fwi, color = vir1,
            mapping = aes(fwi_ori, mu),
            inherit.aes = F) +
  geom_linerange(alpha = 0.5, color = vir2) +
  geom_point(alpha = 0.7, shape = 21, color = vir1, fill = vir2, stroke = 0.35) +

  ggh4x::facet_nested_wrap(vars(ig, ig_location),
                           axes = "all", remove_labels = "all") +

  ylab("Steps (log)") +
  xlab("Fire Weather Index") +

  scale_y_continuous(expand = c(0.05, 0.05)) +

  theme(strip.background = element_rect(color = "white", fill = "white"),
        strip.text = element_text(size = 11),
        panel.spacing.x = unit(4, "mm"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.line = element_line(linewidth = 0.3))
fig_steps

# area ~ steps
npred <- 150
pred_area <- data.frame(steps = seq(min(steps_summ$mean),
                                    max(steps_summ$q95),
                                    length.out = npred))
predmat <- matrix(NA, npred, npost)
for(i in 1:npred) {
  predmat[i, ] <-
    as.vector(draws_arr[, , "fixef_area_a"]) +
    as.vector(draws_arr[, , "fixef_area_b"]) * pred_area$steps[i]

}

# predmat contains the normal mean, but as the area is truncated-normal,
# the true mean is another.
ss <- as.vector(draws_arr[, , "fixef_area_s2"]) |> sqrt()
predmat_mu <- predmat
# # Truncated-normal mean:
# for(i in 1:npost) {
#   predmat_mu[, i] <- etruncnorm(a = areaL, mean = predmat[, i], sd = ss[i])
# }


pred_area$mu <- rowMeans(predmat_mu)
pred_area$mu_lower <- apply(predmat_mu, 1, quantile, prob = 0.025)
pred_area$mu_upper <- apply(predmat_mu, 1, quantile, prob = 0.975)

fig_area <-
ggplot() +
  geom_hline(yintercept = log(10), linetype = "dashed", color = "gray") +

  geom_ribbon(data = pred_area,
              mapping = aes(steps, mu, ymin = mu_lower, ymax = mu_upper),
              inherit.aes = F, color = NA, alpha = 0.4, fill = vir1) +
  geom_line(data = pred_area, color = vir1,
            mapping = aes(steps, mu),
            inherit.aes = F) +

  geom_linerange(data = steps_summ, orientation = "y",
                 mapping = aes(x = mean, y = area_log, xmin = q5, xmax = q95,
                               color = ig_location),
                 alpha = 0.5, color = vir2) +
  geom_point(data = steps_summ,
             mapping = aes(mean, area_log),
             alpha = 0.7, shape = 21, color = vir1, fill = vir2) +

  scale_color_viridis(discrete = TRUE, end = 0.6) +

  facet_wrap(vars(ig_location), axes = "all", axis.labels = "margins") +

  xlab("Steps (log)") +
  ylab("Fire size (log ha)") +

  scale_y_continuous(expand = c(0.05, 0.05)) +

  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing.x = unit(4, "mm"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.line = element_line(linewidth = 0.3))
fig_area

fig_steps_area <- ggarrange2(
  fig_steps + theme(plot.margin = margin(b = 7, unit = "mm")),
  fig_area,
  nrow = 2
)
ggsave("spread/figures_FWIZ/steps_fwi_area.png",
       plot = fig_steps_area,
       width = 17, height = 15, units = "cm")


#######################################################################
# CONTINUE HERE -----------------------------------------------------------
#######################################################################


# Spread probability curves ---------------------------------------------

# Four panels, with spread prob as a function of vfi, tfi, slope and wind.
# In all cases, the remaining predictors are fixed at zero (makes sense?),
# so the FWI affects only the focal slope and the intercept.
# From every posterior sample, 3 curves (3 FWI values) could be drawn, which are
# the average across 200 random-effect curves each.

# start by making a design matrix from the (expand-gridded) prediction
# data.frame. Here, use regional data to define the range for predictors.

# output: array curves_samples[npred, 3, npost]
# For every posterior sample (i),
# 01 - get simulated spread parameters for the three FWI values.
#      array[3, n_coef, 200]
# 02 - compute the array of spread-prob linear predictors, with dimension
#      [npred, 3, 200]. then, apply plogis().
# 03 - average curves across random effects obtaining the average curve as
#      curve_mean[npred, 3]
# 04 - store curves_samples[, , i] <- curve_mean

# Load landscapes to find the 95 % percentile of wind speed.
# Then, bear in mind that in landscapes wind speed was divided by the sd (1.41).
# In the case of slope, it was not scaled.

wind_sd <- 1.46 ## check value later
wind_high_kmh <- 90 ## check later
wind_high_mps <- wind_high_kmh / 3.6
wind_high_z <- wind_high_mps / wind_sd

vfi_low <- fi_params$vfi_z_hdi["vfi_lower"]
vfi_high <- fi_params$vfi_z_hdi["vfi_upper"]

tfi_low <- fi_params$tfi_z_hdi["tfi_lower"]
tfi_high <- fi_params$tfi_z_hdi["tfi_upper"]

nseq <- 400

pdspread <- rbind(
  # vfi ___________
  expand.grid(
    vfi = seq(vfi_low, vfi_high, length.out = nseq),
    tfi = 0,
    slope_ang = 0,
    wind_kmh = 0,
    varying_var = "vfi"
  ),
  # tfi ___________
  expand.grid(
    vfi = 0,
    tfi = seq(tfi_low, tfi_high, length.out = nseq),
    slope_ang = 0,
    wind_kmh = 0,
    varying_var = "tfi"
  ),
  # slope ___________
  expand.grid(
    vfi = 0,
    tfi = 0,
    slope_ang = seq(-45, 45, length.out = nseq),
    wind_kmh = 0,
    varying_var = "slope"
  ),
  # wind ___________
  expand.grid(
    vfi = 0,
    tfi = 0,
    slope_ang = 0,
    wind_kmh = seq(-90, 90, length.out = nseq),
    varying_var = "wind"
  )
)

# scale predictors in the way the fire simulator needs them
names_spread <- c("vfi", "tfi", "slope", "wind")
names_plot <- colnames(pdspread)[1:4]

pdspread$wind <- pdspread$wind_kmh / 3.6 / wind_sd
pdspread$slope <- sin(pdspread$slope_ang * pi / 180)
# turn into zero the slope values below zero
pdspread$slope[pdspread$slope < 0] <- 0

# move those values to varying_val
pdspread$varying_val <- NA
for(v in 1:4) {
  rows <- pdspread$varying_var == names_spread[v]
  colname <- names_plot[v]
  pdspread$varying_val[rows] <- pdspread[rows, colname]
}

# to match later:
npred <- nrow(pdspread)
pdspread$row <- 1:npred

# design matrix
Xspread <- model.matrix(~ vfi + tfi + slope + wind, data = pdspread)

# prediction data for FWI:

# fwi_ref <- quantile(fwi_all, prob = c(0.025, 0.5, 0.975), method = 8)
# fwi_ref_ori <- fwi_ref * fwi_sd + fwi_mean
# range(fwi_all) * fwi_sd + fwi_mean
fwi_ref_ori <- c(0, 11, 32)
fwi_ref <- (fwi_ref_ori - fwi_mean) / fwi_sd

Xfwi <- cbind(rep(1, 3), fwi_ref)
npred_mu <- length(fwi_ref)

# Array to store samples of curves
curves_samples <- array(
  NA, dim = c(npred, npred_mu, npost),
  dimnames = list(
    row = 1:npred,
    fwi_level = fwi_ref_ori,
    iter = 1:npost
  )
)

# Compute curves

nr <- 200 # number of random effects to simulate in order to compute mean
ranef_raw <- matrix(rnorm(nr * (n_coef-1)), ncol = n_coef-1)
# random effects array, placeholders
ranef_tmp <- array(NA, dim = c(npred_mu, n_coef-1, nr))
ranef_cons <- array(NA, dim = c(npred_mu, n_coef-1, nr))

fitted_prob <- array(NA, dim = c(npred, npred_mu, nr))

for(i in 1:npost) {
  if(i %% 100 == 0) print(i)

  # mu at unconstrained scale
  mumat <- Xfwi %*% ab_arr[, 1:(n_coef-1), i]

  # Compute choleski factor of vcov matrix for random effects
  sds <- sigma_mat[i, , drop = T]
  rho_vec <- rho_mat[i, , drop = T]
  rr <- matrix(1, n_coef, n_coef)
  rr[lower.tri(rr)] <- rho_vec
  rr <- t(rr)
  rr[lower.tri(rr)] <- rho_vec
  V <- diag(sds) %*% rr %*% diag(sds)
  Vchol_U <- chol(V[1:(n_coef-1), 1:(n_coef-1)])

  # unconstrained centred random effects
  ranef_centred <- ranef_raw %*% Vchol_U

  # unconstrained random effects
  for(j in 1:nr) {
    ranef_tmp[, , j] <- t(t(mumat) + ranef_centred[j, ])
  }

  # constrain (only logit-link ones, leave steps at log)
  for(v in 1:(n_coef-1)) {
    # v = 1
    ranef_cons[, v, ] <- constrain_vec(ranef_tmp[, v, ], v = par_names[v],
                                         support = support)
  }
  ranef_cons3 <- aperm(ranef_cons, c(2, 1, 3))

  # Compute spread probability curves for all raneffs and FWI values
  for(j in 1:nr) {
    fitted_prob[, , j] <- plogis(Xspread %*% ranef_cons3[, , j])
  }

  # average curves across raneffs
  curves_samples[, , i] <- apply(fitted_prob, 1:2, mean)
}

# summarize posterior and longanize.
curves_summ <- apply(curves_samples, 1:2, summarise)
names(dimnames(curves_summ))[1] <- "metric"

curves_summ_sub <- curves_summ[c("mean", "eti_lower_95", "eti_upper_95"), , ]
curves_df1 <- as.data.frame.table(curves_summ_sub, responseName = "probfit")
curves_df <- pivot_wider(curves_df1, names_from = "metric", values_from = "probfit")
curves_df$row <- as.numeric(as.character(curves_df$row))

# merge with predictions data
curves_df <- left_join(curves_df,
                       pdspread[, c("row", "varying_val", "varying_var")],
                       by = "row")

curves_df$varying_var2 <- factor(as.character(curves_df$varying_var),
                                 levels = c("vfi", "tfi", "slope", "wind"),
                                 labels = c("VFI", "TFI", "Slope (°)",
                                            "Wind speed (km/h)"))

saveRDS(curves_df, file.path("files", "hierarchical_model_FWIZ", "curves_df_prediction.rds"))
curves_df <- readRDS(file.path("files", "hierarchical_model_FWIZ", "curves_df_prediction.rds"))

# Plot

ggplot(curves_df, aes(varying_val, mean,
                      ymin = eti_lower_95, ymax = eti_upper_95,
                      color = fwi_level, fill = fwi_level)) +

  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +

  scale_color_viridis(option = "A", end = 0.8, discrete = T,
                      name = "FWI") +
  scale_fill_viridis(option = "A", end = 0.8, discrete = T,
                     name = "FWI") +

  facet_wrap(vars(varying_var2), scales = "free_x", strip.position = "bottom",
             axes = "all", axis.labels = "margins") +

  ylab("Spread probability") +
  # scale_y_continuous(expand = c(0.05, 0.05)) +
  theme(strip.placement = "outside",
        strip.background = element_rect(color = "white", fill = "white"),
        strip.text = element_text(margin = margin(r = 0, l = 0, unit = "mm"),
                                  size = 11),
        panel.spacing.x = unit(0, "mm"),
        panel.spacing.y = unit(6, "mm"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right")

ggsave("spread/figures_FWIZ/spread_prob_curves.png",
       width = 14, height = 12, units = "cm")



# Spread probability curves (raw variables) -------------------------------

# First, plot the relationship between raw variables and flammability indices.
# Then, plot spraed probability as a function of raw variables
# (ndvi * veg, elevation, slope-weighted northing)

veg_levels <- data_summ$ndvi$vegetation
nveg <- length(veg_levels)

nseq <- 300

pd_veg <- do.call("rbind", lapply(1:nveg, function(v) {
  expand.grid(
    ndvi = seq(data_summ$ndvi$hdi_lower_95[v],
               data_summ$ndvi$hdi_upper_95[v],
               length.out = nseq),
    vegnum = v,
    vegetation = veg_levels[v],
    elevation = 0,
    northing = 0,
    varying_var = "ndvi"
  )
}))

pd_topo <- rbind(
  expand.grid(
    ndvi = 0,
    vegnum = 1,
    vegetation = veg_levels[1],
    elevation = seq(data_summ$elevation["hdi_lower_95"],
                    data_summ$elevation["hdi_upper_95"],
                    length.out = nseq),
    northing = 0,
    varying_var = "elevation"
  ),
  expand.grid(
    ndvi = 0,
    vegnum = 1,
    vegetation = veg_levels[1],
    elevation = data_summ$elevation["mean"],
    northing = seq(-1, 1, length.out = nseq),
    varying_var = "northing"
  )
)

pdspread_fi <- rbind(pd_veg, pd_topo)

names_spread <- c("ndvi", "elevation", "northing")
# move values to varying_val
pdspread_fi$varying_val <- NA
for(v in 1:3) {
  rows <- pdspread_fi$varying_var == names_spread[v]
  colname <- names_spread[v]
  pdspread_fi$varying_val[rows] <- pdspread_fi[rows, colname]
}

# Compute flammability indices
pdspread_fi$vfi <- pdspread_fi$tfi <- 0

# VFI
for(v in 1:nveg) {
  id_fill <- pdspread_fi$vegnum == v
  pdspread_fi$vfi[id_fill] <-
    fi_params$a[v] +
    fi_params$b[v] *
    (pdspread_fi$ndvi[id_fill] - fi_params$o[v]) ^ 2
}
# (standardize)
pdspread_fi$vfi <- (pdspread_fi$vfi - fi_params$vfi_mean) / fi_params$vfi_sd

# TFI
pdspread_fi$tfi <-
  fi_params$b_elev_ori * pdspread_fi$elevation +
  fi_params$b_north_ori * pdspread_fi$northing
pdspread_fi$tfi <- (pdspread_fi$tfi - fi_params$tfi_mean) / fi_params$tfi_sd

# Plot indices

# VFI
vfiplot <-
ggplot(pdspread_fi[pdspread_fi$varying_var == "ndvi", ]) +
  geom_line(aes(ndvi, vfi, color = vegetation), linewidth = 0.6) +
  scale_color_viridis(discrete = T, begin = 0, end = 0.9, option = "D",
                      name = "Vegetation\ntype") +
  ylab("Vegetation\nFlammability Index (VFI)") +
  xlab("NDVI") +
  scale_y_continuous(breaks = seq(-1, 1, 1), limits = c(-2, 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right")

# TFI
# elevation
tfiplot1 <-
ggplot(pdspread_fi[pdspread_fi$varying_var == "elevation", ]) +
  geom_line(aes(elevation, tfi), linewidth = 0.6) +
  ylab("Topographic\nFlammability Index (TFI)") +
  xlab("Elevation (m a.s.l.)") +
  scale_y_continuous(breaks = seq(-2, 2, 1), limits = c(-2, 2)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.3))

# northing
tfiplot2 <-
ggplot(pdspread_fi[pdspread_fi$varying_var == "northing", ]) +
  geom_line(aes(northing, tfi), linewidth = 0.6) +
  xlab("Northing") +
  scale_y_continuous(breaks = seq(-2, 2, 1), limits = c(-2, 2)) +
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        axis.title.y = element_blank())

vegleg <- ggpubr::get_legend(vfiplot) |> ggpubr::as_ggplot()

fi_fig <- ggarrange2(
  vfiplot + theme(legend.position = "none"),
  vegleg, tfiplot1, tfiplot2,
  nrow = 2
)

ggsave("spread/figures_FWIZ/flammability_indices.png",
       width = 14, height = 12, units = "cm",
       plot = fi_fig)


## Compute spread probability curves.

# make zero the non-varying index
pdspread_fi$vfi[pdspread_fi$varying_var != "ndvi"] <- 0
pdspread_fi$tfi[pdspread_fi$varying_var == "ndvi"] <- 0

# to match later:
npred <- nrow(pdspread_fi)
pdspread_fi$row <- 1:npred

# design matrix
Xspread_fi <- model.matrix(~ vfi + tfi, data = pdspread_fi)

# prediction data for FWI:

# fwi_ref <- quantile(fwi_all, prob = c(0.025, 0.5, 0.975), method = 8)
# fwi_ref_ori <- fwi_ref * fwi_sd + fwi_mean
# range(fwi_all) * fwi_sd + fwi_mean
fwi_ref_ori <- c(0, 11, 32)
fwi_ref <- (fwi_ref_ori - fwi_mean) / fwi_sd

Xfwi <- cbind(rep(1, 3), fwi_ref)
npred_mu <- length(fwi_ref)

# Array to store samples of curves
curves_samples <- array(
  NA, dim = c(npred, npred_mu, npost),
  dimnames = list(
    row = 1:npred,
    fwi_level = fwi_ref_ori,
    iter = 1:npost
  )
)

# Compute curves

nr <- 200 # number of random effects to simulate in order to compute mean
n_coef_fi <- 3
ranef_raw <- matrix(rnorm(nr * n_coef_fi), ncol = 3) # only simulates the used parameters
# random effects array, placeholders
ranef_tmp <- array(NA, dim = c(npred_mu, n_coef_fi, nr))
ranef_cons <- array(NA, dim = c(npred_mu, n_coef_fi, nr))

fitted_prob <- array(NA, dim = c(npred, npred_mu, nr))

for(i in 1:npost) {
  if(i %% 100 == 0) print(i)

  # mu at unconstrained scale
  mumat <- Xfwi %*% ab_arr[, 1:n_coef_fi, i]

  # Compute choleski factor of vcov matrix for random effects
  sds <- sigma_mat[i, , drop = T]
  rho_vec <- rho_mat[i, , drop = T]
  rr <- matrix(1, n_coef, n_coef)
  rr[lower.tri(rr)] <- rho_vec
  rr <- t(rr)
  rr[lower.tri(rr)] <- rho_vec
  V <- diag(sds) %*% rr %*% diag(sds)
  Vchol_U <- chol(V[1:3, 1:3])

  # unconstrained centred random effects
  ranef_centred <- ranef_raw %*% Vchol_U

  # unconstrained random effects
  for(j in 1:nr) {
    ranef_tmp[, , j] <- t(t(mumat) + ranef_centred[j, ])
  }

  # constrain (only logit-link ones, leave steps at log)
  for(v in 1:n_coef_fi) {
    # v = 1
    ranef_cons[, v, ] <- constrain_vec(ranef_tmp[, v, ], v = par_names[v],
                                         support = support)
  }

  ranef_cons3 <- aperm(ranef_cons, c(2, 1, 3))

  # Compute spread probability curves for all raneffs and FWI values
  for(j in 1:nr) {
    fitted_prob[, , j] <- plogis(Xspread_fi %*% ranef_cons3[, , j])
  }

  # average curves across raneffs
  curves_samples[, , i] <- apply(fitted_prob, 1:2, mean)
}

# summarize posterior and longanize.
curves_summ <- apply(curves_samples, 1:2, summarise)
names(dimnames(curves_summ))[1] <- "metric"

curves_summ_sub <- curves_summ[c("mean", "eti_lower_95", "eti_upper_95"), , ]
curves_df1 <- as.data.frame.table(curves_summ_sub, responseName = "probfit")
curves_df <- pivot_wider(curves_df1, names_from = "metric", values_from = "probfit")
curves_df$row <- as.numeric(as.character(curves_df$row))

# merge with predictions data
curves_df <- left_join(curves_df,
                       pdspread_fi[, c("row", "vegetation", "varying_val", "varying_var")],
                       by = "row")

saveRDS(curves_df, file.path("files", "hierarchical_model_FWIZ", "curves_df_prediction_raw_x.rds"))
curves_df <- readRDS(file.path("files", "hierarchical_model_FWIZ", "curves_df_prediction_raw_x.rds"))


# Plots
bgcol <- "#1a1a1aff"

curves_df$varying_var2 <- as.character(curves_df$varying_var)
curves_df$varying_var2[curves_df$varying_var == "elevation"] <-
  "Elevation (m a.s.l.)"
curves_df$varying_var2[curves_df$varying_var == "northing"] <-
  "Northing"

# NDVI * veg
curves_veg <-
ggplot(curves_df[curves_df$varying_var == "ndvi", ],
       aes(varying_val, mean,
           ymin = eti_lower_95, ymax = eti_upper_95,
           color = fwi_level, fill = fwi_level)) +

  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +

  scale_color_viridis(option = "A", end = 0.8, discrete = T,
                      name = "FWI") +
  scale_fill_viridis(option = "A", end = 0.8, discrete = T,
                     name = "FWI") +

  facet_wrap(vars(vegetation), scales = "fixed", strip.position = "top",
             axes = "all", axis.labels = "margins", nrow = 3) +

  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1),
                     expand = c(0.005, 0.005)) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +

  xlab("NDVI") +
  theme(strip.placement = "outside",
        strip.background = element_rect(color = bgcol, fill = bgcol),
        strip.text = element_text(#margin = margin(r = 0, l = 0, unit = "mm"),
                                  size = 10, color = "white"),
        panel.spacing.x = unit(3, "mm"),
        panel.spacing.y = unit(6, "mm"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        axis.title.y = element_blank(),
        legend.position = c(0.78, 0.13))
curves_veg

# Elevation and topography
curves_topo <-
ggplot(curves_df[curves_df$varying_var != "ndvi", ],
       aes(varying_val, mean,
           ymin = eti_lower_95, ymax = eti_upper_95,
           color = fwi_level, fill = fwi_level)) +

  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +

  scale_color_viridis(option = "A", end = 0.8, discrete = T,
                      name = "FWI") +
  scale_fill_viridis(option = "A", end = 0.8, discrete = T,
                     name = "FWI") +

  facet_wrap(vars(varying_var2), scales = "free_x", strip.position = "bottom",
             axes = "all", axis.labels = "margins", nrow = 1) +

  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1),
                     expand = c(0.005, 0.005)) +
  theme(strip.placement = "outside",
        strip.background = element_rect(color = "white", fill = "white"),
        strip.text = element_text(margin = margin(r = 0, l = 0, unit = "mm"),
                                  size = 11),
        panel.spacing.x = unit(3, "mm"),
        panel.spacing.y = unit(6, "mm"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none")
curves_topo


ylabb <- grid::textGrob("Spread probability",
                        gp = grid::gpar(fontsize = 11, fontface = "plain"),
                        rot = 90)

rawx_fig <-
ggarrange2(curves_veg + theme(plot.margin = margin(b = 5, unit = "mm")),
           curves_topo,
           left = ylabb, #"Spread probability",
           label.args = list(gp = gpar(size = 11)),
           nrow = 2, heights = c(3.9, 1))

ggsave("spread/figures_FWIZ/spread_prob_curves_raw_x.png",
       width = 14, height = 18, units = "cm",
       plot = rawx_fig)


# Assessing model fit -----------------------------------------------------

# Compute overlap distribution and dharma residuals for total fire size and
# by veg type, conditioning on fitted ranef and simulating new ones.

nsim <- 2000 # simulations
ids_sim <- sample(1:npost, nsim, replace = F)

# Prepare parameters _______________

# Tidy fitted random effects
ranef_names <- data.frame(nn = colnames(ranef_mat))
ranef_names <- separate(ranef_names, "nn",
                        into = c("type", "par_name", "fire_num"))

ranef_fit <- array(NA, dim = c(npost, n_coef, J1),
                   dimnames = list(
                     iter = 1:npost,
                     par_name = par_names,
                     fire_num = 1:J1
                   ))
ranef_sim <- ranef_fit # to fill later

for(j in 1:J1) {
  for(p in 1:n_coef) {
    # j = 1; p = 2
    col <-
      ranef_names$par_name == par_names[p] &
      ranef_names$fire_num == j
    ranef_fit[, p, j] <- ranef_mat[, col]
  }
}
# now, constrain!
for(j in 1:J1) {
  ranef_fit[, , j] <- constrain(ranef_fit[, , j], support)
  ss <- ranef_fit[, "steps", j]    # steps below 1 are 1, as zero means unlimited
  ranef_fit[ss < 1, "steps", j] <- 1
}

# Get samples of random effects sds
sigma_mat <- fixef_mat[, grep("s2", colnames(fixef_mat))]
sigma_mat <- sigma_mat[, 1:n_coef] |> sqrt() # remove area and turn into sd

# Get samples of intercepts and slopes for FWI regressions
ab_arr <- abind::abind(fixef_mat[, grep("_a$", colnames(fixef_mat)), drop = T],
                       fixef_mat[, grep("_b$", colnames(fixef_mat)), drop = T],
                       along = 3)
ab_arr <- ab_arr[, 1:n_coef, ] # remove area parameters
dimnames(ab_arr) <- list(
  iter = 1:npost,
  param = par_names,
  coef = c("a", "b")
)
ab_arr <- aperm(ab_arr, perm = c(3, 2, 1))
ab_arr[, , 1]

# tidy correlation matrices
lays <- matrix(NA, n_coef, n_coef)
lays[lower.tri(lays)] <- lays[upper.tri(lays)] <- 1:ncol(combn(1:n_coef, 2))

# raw random effects
ranef_raw <- matrix(rnorm(J1 * n_coef), ncol = n_coef)
colnames(ranef_raw) <- par_names

# Loop to compute means
for(i in 1:npost) {
  # i = 1
  if(i %% 1000 == 0) print(i)

  # mu at unconstrained scale
  mumat <- X_ %*% ab_arr[, , i]

  # Compute choleski factor of vcov matrix for random effects
  sds <- sigma_mat[i, , drop = T]
  rho_vec <- rho_mat[i, , drop = T]
  rr <- matrix(1, n_coef, n_coef)
  rr[lower.tri(rr)] <- rho_vec
  rr <- t(rr)
  rr[lower.tri(rr)] <- rho_vec
  V <- diag(sds) %*% rr %*% diag(sds)
  Vchol_U <- chol(V)

  # Compute actual ranef
  ranef_centred <- ranef_raw %*% Vchol_U
  ranef_unc <- mumat + ranef_centred
  ranef_cons <- constrain(ranef_unc, support)
  ss <- ranef_cons[, "steps"]    # steps below 1 are 1, as zero means unlimited
  ranef_cons[ss < 1, "steps"] <- 1

  ranef_sim[i, , ] <- t(ranef_cons)
}

# subset only nsim posterior samples in both arrays
ranef_fit <- ranef_fit[ids_sim, , ]
ranef_sim <- ranef_sim[ids_sim, , ]

# merge to simplify paralelization
ranef_both <- abind::abind(ranef_fit, ranef_sim, along = 1)

# Matrices to fill with fire metrics
metrics_table <- array(NA, dim = c(nsim * 2, nmet, J1),
                       dimnames = list(
                         iter = c(paste("fit", 1:nsim, sep = "_"),
                                  paste("sim", 1:nsim, sep = "_")),
                         metric = met_names,
                         fire_num = 1:J1
                       ))

size_obs <- matrix(NA, J1, n_veg + 1)
veg_available <- matrix(NA, J1, n_veg) # available veg types by fire, to control
# residuals

gc()
registerDoMC(cores = 12)

# loop over fires
for(j in 1:J1) {
  # j = 56
  fire_name <- fires_data_spread$fire_id[j]
  print(paste(j, ": ", fire_name, sep = ""))
  file_name <- paste(fire_name, ".rds", sep = "")

  full_data <- readRDS(file.path(lands_dir, file_name))

  # subset data needed for spread (to be cloned across workers)
  spread_data <- full_data[c("landscape", "ig_rowcol",
                             "burned_layer", "burned_ids",
                             "counts_veg", "counts_veg_available")]


  size_obs[j, 1] <- sum(spread_data$counts_veg)
  size_obs[j, 2:(n_veg+1)] <- spread_data$counts_veg
  veg_available[j, ] <- spread_data$counts_veg_available

  # metrics_table[, , j] <- simulate_metrics_parallel(ranef_both[, , j], spread_data)
  # rm(full_data, spread_data)
  # gc()
}
# saveRDS(metrics_table, file.path("files", "hierarchical_model_FWIZ", "metrics_table.rds"))
# saveRDS(size_obs, file.path("files", "hierarchical_model_FWIZ", "size_obs.rds"))
# saveRDS(veg_available, file.path("files", "hierarchical_model_FWIZ", "veg_available.rds"))
metrics_table <- readRDS(file.path("files", "hierarchical_model_FWIZ", "metrics_table.rds"))

# Overlap for fitted and simulated random effects
ids_fit <- 1:nsim
ids_sim <- (nsim+1):(nsim*2)

ov_fit_summ <- apply(t(metrics_table[ids_fit, "overlap", ]), 1,
                     summarise) |> t() |> as.data.frame()
ov_sim_summ <- apply(t(metrics_table[ids_sim, "overlap", ]), 1,
                     summarise) |> t() |> as.data.frame()

cc <- c("mean", "hdi_lower_95", "hdi_upper_95")
#cc <- c("mean", "eti_lower_95", "eti_upper_95")
ovfit <- ov_fit_summ[, cc]
ovsim <- ov_sim_summ[, cc]
colnames(ovfit) <- paste(c("mean", "lower", "upper"), "fit", sep = "_")
colnames(ovsim) <- paste(c("mean", "lower", "upper"), "sim", sep = "_")

ovtable <- cbind(ovfit, ovsim)

ggplot(ovtable) +
  geom_linerange(aes(x = mean_sim, y = mean_fit, ymin = lower_fit, ymax = upper_fit),
                 alpha = 0.5) +
  geom_linerange(aes(y = mean_fit, xmin = lower_sim, xmax = upper_sim),
                 orientation = "y", alpha = 0.5) +
  geom_point(aes(mean_sim, mean_fit), size = 2) +
  scale_y_continuous(limits = c(0, 0.9), expand = c(0.005, 0.005)) +
  scale_x_continuous(limits = c(0, 0.6), expand = c(0.01, 0.01)) +
  coord_fixed() +
  nice_theme() +
  theme() +
  ylab("Overlap from fitted parameters") +
  xlab("Overlap from simulated parameters")

ggsave("spread/figures_FWIZ/overlap_fit_sim.png",
       width = 9, height = 12, units = "cm")



# Dharma for size ______________________

size_obs2 <- rbind(size_obs, size_obs)
size_sim <- metrics_table[, -1, ]
str(size_sim)
size_sim <- aperm(size_sim, c(3, 1, 2))
size_sim2 <- abind::abind(size_sim[, ids_fit, ], size_sim[, ids_sim, ],
                          along = 1)
str(size_sim2)
sort(rnorm(10))

rows_fit <- 1:J1
rows_sim <- (J1+1):(J1*2)

veg_levels_all <- c("All vegetation types", veg_levels)
unif_q <- ppoints(J1)

res_table <- do.call("rbind", lapply(1:(nmet-1), function(m) {
  # m = 1
  rr <- createDHARMa(size_sim2[, , m], size_obs2[, m], integerResponse = T)
  out <- data.frame(
    veg_class = veg_levels_all[m],
    q_obs = c(
      sort(rr$scaledResiduals[rows_fit]),
      sort(rr$scaledResiduals[rows_sim])
    ),
    q_exp = c(unif_q, unif_q),
    ranef = rep(c("Fitted parameters", "Simulated parameters"), each = J1)
  )
  return(out)
}))

res_table$veg_class <- factor(as.character(res_table$veg_class),
                              levels = c(veg_levels, veg_levels_all[1]))


bgcol <- "#1a1a1aff"

ggplot(res_table, aes(q_exp, q_obs, color = ranef, fill = ranef)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5, linetype = "dashed") +
  geom_point(shape = 21, size = 2, stroke = 0.5) +
  scale_color_viridis(discrete = T, option = "D", end = 0.5) +
  scale_fill_viridis(discrete = T, option = "D", end = 0.5, alpha = 0.4) +
  facet_wrap(vars(veg_class), nrow = 2, axes = "all", axis.labels = "margins") +
  nice_theme() +
  coord_fixed() +
  ylab("Observed quantiles") +
  xlab("Expected quantiles") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(margin = margin(l = -1, unit = "mm"),
                                   size = 10),
        panel.spacing.y = unit(5, "mm"),
        panel.spacing.x = unit(5, "mm"),
        strip.background = element_rect(color = bgcol, fill = bgcol),
        strip.text = element_text(color = "white", size = 11))

ggsave("spread/figures_FWIZ/dharma_size_fit_sim.png",
       width = 15, height = 13, units = "cm")
# quitar los fuegos menores a 10 ha antes de hacer los dharma?
