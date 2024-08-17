# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(bayesplot)
library(bayestestR) # hdi
library(posterior)
library(deeptime)
library(terra)
library(kdevine)
library(truncnorm)
theme_set(theme_bw())

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

# simulator of random effects, based on the ranef_kde list (a list with a
# kdevine density for each fire.)
ranef_rng <- function(ids) {
  x <- sapply(ids, function(i) {
    rkdevine(1, ranef_kde[[i]])
  })
  return(x)
}

# climatic data ------------------------------------------------------------

fwi_data <- read.csv(file.path("data", "climatic_data_by_fire_fwi-fortnight-cumulative.csv"))

# fire polygons (to get area) ---------------------------------------------

ff <- vect("data/patagonian_fires_spread.shp")
ff$area_ha <- expanse(ff) / 1e4 # turn m2 to ha

# A few constants ---------------------------------------------------------

par_names <- c("wet", "subalpine", "dry", "shrubland", "grassland",
               "slope", "wind", "steps")
n_coef <- length(par_names)
# load file with constants to standardize
ndvi_params <- readRDS(file.path("data", "NDVI_regional_data",
                                 "ndvi_optim_and_proportion.rds"))
slope_sd <- ndvi_params$slope_term_sd

# support for parameters
n_veg <- 5
n_terrain <- 2

ext_alpha <- 50
ext_beta <- 30

params_lower <- c(rep(-ext_alpha, n_veg), rep(0, n_terrain))
params_upper <- c(rep(ext_alpha, n_veg), ext_beta / slope_sd, ext_beta)

support <- rbind(params_lower, params_upper)
colnames(support) <- names(params_lower) <- names(params_upper) <- par_names[-n_coef]
support_width <- apply(support, 2, diff)

# Load and prepare data ---------------------------------------------------

# dir to load files
target_dir <- file.path("files", "overlaps")

samples_files <- list.files(target_dir, pattern = "-posterior_samples.rds")

# import a list with matrix-samples
samples_list <- lapply(samples_files, function(s) {
  rrr <- readRDS(file.path(target_dir, s))
  return(rrr)
})

fire_ids <- lapply(samples_list, function(x) attr(x, "fire_id")) %>% unlist()
names(samples_list) <- fire_ids

samples_df <- do.call("rbind", lapply(samples_list, function(x) {
  fire_id <- attr(x, "fire_id")
  r <- cbind(
    data.frame(fire_id = fire_id),
    as.data.frame(x)
  )
  return(r)
}))

# make list of samples in the unconstrained scale
samples_list_unc <- lapply(samples_list, function(x) {
  x2 <- unconstrain(x, support)
})

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

fwi_mean <- mean(fires_data$fwi_expquad_fortnight)
fwi_sd <- sd(fires_data$fwi_expquad_fortnight)

fires_data$fwi <- (fires_data$fwi_expquad_fortnight - fwi_mean) / fwi_sd
rownames(fires_data) <- fires_data$fire_id

fires_data_spread <- fires_data[fires_data$fire_id %in% fire_ids, ]
fires_data_spread <- fires_data_spread[fire_ids, ]

fires_data_nonspread <- fires_data[!(fires_data$fire_id %in% fire_ids), ]

fire_ids_nonspread <- fires_data_nonspread$fire_id
nfires_nonspread <- length(fire_ids_nonspread)
nfires_spread <- length(fire_ids)

# Get bounds for steps parameters
steps_bounds <- do.call("rbind", lapply(samples_list, function(x) {
  bb <- attr(x, "support")
  return(bb[, "steps"])
}))


# Brute estimates ---------------------------------------------------------

# using samples from the first-step abc, get a point estimate of parameters
# across fires. Those are:

# intercepts, mean and sd
# alpha, beta and sd for log_steps ~ log_area, with
#   predicted steps for non-spread fires (to use as starting values)
# alpha, beta and sd for (slope, wind, steps) ~ fwi

par_start <- vector("list", n_coef + 2)
names(par_start) <- c(par_names[-n_coef],
                      "steps_area", "steps_fwi", "area_steps")
nsim <- 5000

# make matrices for intercepts
for(v in par_names[1:5]) {
  mm <- matrix(NA, nsim, 2)
  colnames(mm) <- c("alpha", "sigma")
  par_start[[v]] <- mm
}

# make matrices for slope, wind and steps
for(v in c("slope", "wind", "steps_fwi", "area_steps")) {
  mm <- matrix(NA, nsim, 3)
  colnames(mm) <- c("alpha", "beta", "sigma")
  par_start[[v]] <- mm
}

# make matrices for steps ~ area
for(v in "steps_area") {
  mm1 <- matrix(NA, nsim, 3)
  colnames(mm1) <- c("alpha", "beta", "sigma")

  mm2 <- matrix(NA, nsim, nfires_nonspread)
  colnames(mm2) <- fire_ids_nonspread

  par_start[[v]] <- list(fixef = mm1, ranef = mm2)
}

# Turn samples into array
samples_temp <- lapply(samples_list_unc, function(x) x[1:nsim, ])
samples_arr <- abind::abind(samples_temp, along = 3)

# covariates with all fires
fwi_full <- c(fires_data_spread$fwi, fires_data_nonspread$fwi)
area_full <- c(fires_data_spread$area_ha_log, fires_data_nonspread$area_ha_log)

# # Compute estimates
# for(i in 1:nsim) {
#   # i = 1
#   print(i)
#
#   # intercepts
#   for(v in par_names[1:5]) {
#     par_start[[v]][i, "alpha"] <- mean(samples_arr[i, v, ])
#     par_start[[v]][i, "sigma"] <- sd(samples_arr[i, v, ])
#   }
#
#   # slope, wind
#   for(v in c("slope", "wind")) {
#     parvals <- samples_arr[i, v, ]
#     mod <- lm(parvals ~ fwi, data = fires_data_spread)
#     par_start[[v]][i, ] <- c(coef(mod), sigma(mod))
#   }
#
#   # steps ~ area model
#   # used to provide initial values for the non-spread steps parameters
#   ss <- samples_arr[i, "steps", ]
#   mod <- lm(ss ~ area_ha_log, data = fires_data_spread)
#   par_start[["steps_area"]][["fixef"]][i, ] <- c(coef(mod), sigma(mod)) # will not be used
#   # predict steps in non-spread fires
#   steps_pred <- predict(mod, fires_data_nonspread)
#   par_start[["steps_area"]][["ranef"]][i, ] <- steps_pred
#
#   # estimate steps ~ fwi and area ~ steps using all fires
#   ss_full <- c(ss, steps_pred)
#   # steps ~ fwi
#   mod <- lm(ss_full ~ fwi_full)
#   par_start[["steps_fwi"]][i, ] <- c(coef(mod), sigma(mod))
#   # area ~ steps
#   mod <- lm(ss_full ~ area_full)
#   par_start[["area_steps"]][i, ] <- c(coef(mod), sigma(mod))
# }
# saveRDS(par_start, file.path("files", "hierarchical_model", "par_start.rds"))
par_start <- readRDS(file.path("files", "hierarchical_model", "par_start.rds"))


# Compute kdes for random effects (spread parameters) ----------------------

# ranef_kde <- lapply(1:nfires_spread, function(i) {
#   print(i)
#   sup_steps_log <- attr(samples_list[[i]], "support")[, "steps"] %>% log
#   xx <- samples_temp[[i]] # unconstrained
#   kde <- kdevine(
#     xx,
#     xmin = c(rep(-Inf, n_coef-1), sup_steps_log[1]),
#     xmac = c(rep(Inf, n_coef-1), sup_steps_log[2]),
#     copula.type = "kde",
#     cores = 6
#   )
#   return(kde)
# })
# saveRDS(ranef_kde, file.path("files", "hierarchical_model", "ranef_kde.rds"))
ranef_kde <- readRDS(file.path("files", "hierarchical_model", "ranef_kde.rds"))

# Visualize starting values -----------------------------------------------

# plot(density(par_start$subalpine[, "alpha"]))
# lines(density(par_start$shrubland[, "alpha"]), col = 2)
# lines(density(par_start$wet[, "alpha"]), col = 3)
# lines(density(par_start$dry[, "alpha"]), col = 4)
# lines(density(par_start$grassland[, "alpha"]), col = 5)

# MCMC for a single parameter ---------------------------------------------

# (alpha subalpine)

kde_list <- lapply(samples_temp, function(x) {
  kk <- kde1d(x[, "subalpine"])
})

# Copying code from Hooten and Hefley (2019) book.

# setup variables
nsim <- 500
J <- length(kde_list)
mu_hat = rep(0, nsim)
sigma_hat = rep(0, nsim)
muj_hat = matrix(0, J, nsim)

# prior parameters
q = 1       # for sigma
r = 1000
mu_0 = 0    # for mu
s2_0 = 5^2

# initial values
muj = sapply(kde_list, function(x) rkde1d(1, x))
mu = mean(sapply(kde_list, function(x) rkde1d(1, x)))
s2 = var(sapply(kde_list, function(x) rkde1d(1, x)))

# mcmc loop

for(k in 1:nsim ) {
  print(k)

  #### Sample s2 from the conditional distribution
  q_tmp = J / 2 + q
  r_tmp = 1 / (sum((muj - mu) ^ 2) / 2 + 1 / r)
  s2 = 1 / rgamma(1, q_tmp, 1 / r_tmp) # r = 1/beta = scale

  #### sample mu from the conditional distribution
  tmp_var = 1 / (J / s2 + 1 / s2_0)
  tmp_mn = tmp.var * (sum(muj) / s2 + mu_0 / s2_0)
  mu = rnorm(1, tmp_mn, sqrt(tmp_var))

  ### Sample mu_j (all at once!)
  muj_star = sapply(kde_list, function(x) rkde1d(1, x))
    # proposal taken from previous density, fitted with kde

  mh1 = dnorm(muj_star, mu, sqrt(s2), log = TRUE) +
          dlogis(muj, 0, 1, log = TRUE) # prior from stage 1 for previous sample
  mh2 = dnorm(muj, mu, sqrt(s2), log = TRUE) +
          dlogis(muj_star, 0, 1, log = TRUE) # prior from stage 1 for proposal
  keep_idx = exp(mh1 - mh2) > runif(J)

  muj[keep_idx] = muj_star[keep_idx]

  # save samples
  mu_hat[k] = mu
  sigma_hat[k] = sqrt(s2)
  muj_hat[, k] = muj
}
# va lentazoooooo!!!!

matplot(1:nsim, t(muj_hat[1:5, ]), type = "l")
plot(mu_hat, type = "l")
plot(sigma_hat, type = "l")

# MCMC for FWI-controlled parameter (wind) --------------------------------

kde_list <- lapply(samples_temp, function(x) {
  kk <- kde1d(x[, "wind"])
})

# Copying code from Hooten and Hefley (2019) book.

# setup variables
nsim <- 500
J <- length(kde_list)
alpha_hat = numeric(nsim)
beta_hat = numeric(nsim)
sigma_hat = numeric(nsim)
muj_hat = matrix(0, J, nsim)

# prior parameters
q = 1       # for sigma
r = 1000

a_mu_prior = 0    # for alpha
a_sd_prior = 5 ^ 2

b_mu_prior = 0    # for beta
b_sd_prior = 3 ^ 2

# proposal scales for a and b
a_sd_jump <- sd(par_start$wind[, "alpha"])
b_sd_jump <- sd(par_start$wind[, "beta"])

# initial values
muj_hat[, 1] = sapply(kde_list, function(x) rkde1d(1, x))
id_start <- sample(1:5000, 1)
bb_start <- par_start$wind[id_start, ]
alpha_hat[1] = bb_start["alpha"]
beta_hat[1] = bb_start["beta"]
sigma_hat[1] = bb_start["sigma"]

# mcmc loop

for(k in 2:nsim) {
  # k = 2
  print(k)

  #### Sample s2 from the conditional distribution

  # compute mu_fitted, from the fwi, to get error terms
  mu_fitted <- alpha_hat[k-1] + beta_hat[k-1] * fires_data_spread$fwi

  q_tmp = J / 2 + q
  r_tmp = 1 / (sum((muj_hat[, k-1] - mu_fitted) ^ 2) / 2 + 1 / r)
  s2 = 1 / rgamma(1, q_tmp, 1 / r_tmp) # r = 1/beta = scale

  #### sample alpha and beta using metropolis update
  alpha_star = rnorm(1, alpha_hat[k-1], a_sd_jump)
  mu_star = alpha_star + beta_hat[k-1] * fires_data_spread$fwi
  lp1 <- sum(dnorm(muj_hat[, k-1], mu_star, sqrt(s2), log = T)) + # likelihood of random effects
         dnorm(alpha_star, a_mu_prior, a_sd_prior, log = T)       # prior
  lp2 <- sum(dnorm(muj_hat[, k-1], mu_fitted, sqrt(s2), log = T)) + # likelihood of random effects
         dnorm(alpha_hat[k-1], a_mu_prior, a_sd_prior, log = T)   # prior
  alpha <- ifelse(exp(lp1 - lp2) > runif(1), alpha_star, alpha_hat[k-1])
  # update mu with new alpha
  mu_fitted <- alpha + beta_hat[k-1] * fires_data_spread$fwi

  beta_star <- rnorm(1, beta_hat[k-1], b_sd_jump)
  mu_star = alpha + beta_star * fires_data_spread$fwi
  lp1 <- sum(dnorm(muj_hat[, k-1], mu_star, sqrt(s2), log = T)) + # likelihood of random effects
         dnorm(beta_star, b_mu_prior, b_sd_prior, log = T)
  lp2 <- sum(dnorm(muj_hat[, k-1], mu_fitted, sqrt(s2), log = T)) + # likelihood of random effects
         dnorm(beta_hat[k-1], b_mu_prior, b_sd_prior, log = T)
  beta <- ifelse(exp(lp1 - lp2) > runif(1), beta_star, beta_hat[k-1])
  # update mu_fitted with new beta
  mu_fitted <- alpha + beta * fires_data_spread$fwi

  ### Sample proposal for mu_j (vectorized)
  muj_star = sapply(kde_list, function(x) rkde1d(1, x))
  # proposals sampled from stage1 empirical density.
  # (the metropolis dividends are likelihood * other's_prior)
  lp1 = dnorm(muj_star, mu_fitted, sqrt(s2), log = TRUE) +
        dlogis(muj_hat[, k-1], 0, 1, log = TRUE) # prior from stage 1 for previous sample
  lp2 = dnorm(muj_hat[, k-1], mu_fitted, sqrt(s2), log = TRUE) +
        dlogis(muj_star, 0, 1, log = TRUE) # prior from stage 1 for proposal
  keep_idx = exp(lp1 - lp2) > runif(J)

  # save samples
  muj = muj_hat[, k-1]
  muj[keep_idx] = muj_star[keep_idx]

  muj_hat[, k] = muj
  alpha_hat[k] = alpha
  beta_hat[k] = beta
  sigma_hat[k] = sqrt(s2)
}
# va lentazoooooo!!!!

matplot(1:nsim, t(muj_hat[1:5, ]), type = "l")
plot(alpha_hat, type = "l")
plot(beta_hat, type = "l")
plot(sigma_hat, type = "l")


# exp-uniform lpdf for the prior -------------------------------------------

# In the first stage, steps were assigned a flat prior, between 5 and S, with S
# varying among fires. As we aim to estimate the steps at the log scale, we
# need to know which was the prior distribution at the log scale, which is not
# a uniform, but an exponential-uniform distribution with pdf:

dexpunif <- function(x, l = 1, u = 10) { # for x \in [log(l), log(u)]
  d <- ifelse(x >= log(l) & x <= log(u),
              exp(x) / (u-l), 0)
  return(d)
}

expunif_lpdf <- function(x, l = 1, u = 10) { # for x \in [log(l), log(u)]
  ld <- ifelse(x >= log(l) & x <= log(u),
               x, -Inf)
  return(ld)
}

# Old code below in this section

curve(dexpunif(x, 2, 6), from = log(2-1), to = log(6+1))

p1 <- runif(10000, 5, 38)
p2 <- log(p1)
plot(density(p1, from = 5, to = 38))
plot(density(p2, from = log(5), to = log(38)))

ktest <- kde1d(p2, xmin = log(5), xmax = log(38), mult = 3)
curve(dkde1d(x, ktest), from = log(5), to = log(38))

steps_log_prior_kde <- lapply(samples_list, function(x) {
  # x = samples_list[[1]]
  sup <- attr(x, "support")
  bounds <- sup[, "steps"]
  log_bounds <- log(bounds)
  unif_samples <- runif(10000, bounds[1], bounds[2])
  log_samples <- log(unif_samples)
  kde <- kde1d(log_samples, xmin = log_bounds[1], xmax = log_bounds[2], mult = 3)
  attr(kde, "fire_id") <- attr(x, "fire_id")
  # curve(dkde1d(x, kde), from = log_bounds[1], to = log_bounds[2],
  #       main = attr(x, "fire_id"))
  return(kde)
})

steps_log_prior_d <- function(x) {
  ll <- length(x)
  ld <- numeric(ll)
  for(i in 1:ll) {
    ld[i] <- dkde1d(x[i], steps_log_prior_kde[[i]]) %>% log
  }
  return(ld)
}

# prior to evaluate single value, i indicates the spread-fire id.
steps_prior_lpdf <- function(x, i) {
  ld <- dkde1d(x, steps_log_prior_kde[[i]]) %>% log
  return(ld)
}

steps_kde <- lapply(samples_temp, function(x) {
  kk <- kde1d(x[, "steps"])
})

steps_rng <- function() {
  ss <- numeric(nfires_spread)
  for(i in 1:nfires_spread) {
    prior_dens <- 0
    while(prior_dens == 0) {
      stry <- rkde1d(1, steps_kde[[i]])
      prior_dens <- dkde1d(stry, steps_log_prior_kde[[i]])
    }
    ss[i] <- stry
  }
  return(ss)
}

## Test steps rng and prior density
# all(is.finite(steps_log_prior_d(steps_rng())))
# replicate(100, all(is.finite(steps_log_prior_d(steps_rng()))))

# MCMC for steps and area -----------------------------------------------

# parameters are named according to the variable they control.
# The model for area is
#   area_all[i] ~ Normal(area_mu[i], area_s)T[0, )
#   area_mu[i] = area_a + area_b * steps[i]
# The model for steps is
#   steps[i] ~ Normal(steps_mu[i], steps_s)
#   steps_mu[i] = steps_a + steps_b * fwi_all[i]


# lower bound for fire area (10 ha)
areaL <- log(10)

# Copying code from Hooten and Hefley (2019) book.

# setup variables
nsim <- 500
J1 <- length(steps_kde)
J2 <- nfires_nonspread
J <- J1 + J2
ids1 <- 1:J1
ids2 <- (J1+1):J

fwi_all <- c(fires_data_spread$fwi, fires_data_nonspread$fwi)
area_all <- c(fires_data_spread$area_ha_log, fires_data_nonspread$area_ha_log)

# steps parameters
steps_a = numeric(nsim)
steps_b = numeric(nsim)
steps_s = numeric(nsim)
steps = matrix(NA, J, nsim)

# area parameters
area_a = numeric(nsim)
area_b = numeric(nsim)
area_s = numeric(nsim)

# prior parameters
q = 1       # for sigma, all sigmas, non-informative, conjugate
r = 1000

# steps regression (steps ~ fwi)
steps_a_mu_prior = mean(par_start$steps_fwi[, "alpha"])  # for alpha
steps_a_sd_prior = 5
steps_b_mu_prior = 0    # for beta
steps_b_sd_prior = 5

# area regression (area ~ steps)
area_a_mu_prior = mean(par_start$area_steps[, "alpha"])    # for alpha
area_a_sd_prior = 5
area_b_mu_prior = 0    # for beta
area_b_sd_prior = 5
area_s_mu_prior = 0    # for sigma
area_s_sd_prior = 5    # updated with M-H, for the likelihood is truncated-
                       # normal.

# proposal scales for a and b
steps_a_sd_jump <- sd(par_start$steps_fwi[, "alpha"]) * 2
steps_b_sd_jump <- sd(par_start$steps_fwi[, "beta"]) * 2
steps_sd_jump <- apply(par_start$steps_area$ranef, 2, sd) * 2

# for area
area_a_sd_jump <- sd(par_start$area_steps[, "alpha"]) * 2
area_b_sd_jump <- sd(par_start$area_steps[, "beta"]) * 2
area_s_sd_jump <- sd(par_start$area_steps[, "sigma"]) * 2

# initial values
id_start <- sample(1:5000, 1)

steps[ids1, 1] <- steps_rng()
steps[ids2, 1] <- par_start$steps_area$ranef[id_start, ]

steps_coef_start <- par_start$steps_fwi[id_start, ]
steps_a[1] = steps_coef_start["alpha"]
steps_b[1] = steps_coef_start["beta"]
steps_s[1] = steps_coef_start["sigma"]
# for area
area_coef_start <- par_start$area_steps[id_start, ]
area_a[1] = area_coef_start["alpha"]
area_b[1] = area_coef_start["beta"]
area_s[1] = area_coef_start["sigma"]

# mcmc loop

# "_" suffix indicates parameters already updated.
# "_try" identifies proposals to be evaluated.

for(k in 2:nsim) {
  # k = 2
  print(k)

  #### Update parameters for area ~ steps regression

  # As the likelihood is truncated-normal, the conditional distribution for
  # area_s (sigma) is unknown. Use M-H
  area_mu <- area_a[k-1] + area_b[k-1] * steps[, k-1]
  area_s_try <- rtruncnorm(1, a = 0, mean = area_s[k-1], sd = area_s_sd_jump)

  lp1 <-
    sum(dtruncnorm(area_all, a = areaL,
                   mean = area_mu, sd = area_s_try) %>% log) +         # likelihood
    dtruncnorm(area_s_try, a = 0,
               mean = area_s_mu_prior, sd = area_s_sd_prior) %>% log + # prior
    dtruncnorm(area_s[k-1], a = 0,
               mean = area_s_try, sd = area_s_sd_jump) %>% log         # jump
  # this terms takes into account the asymmetric proposal (truncnorm)

  lp2 <-
    sum(dtruncnorm(area_all, a = areaL,
                   mean = area_mu, sd = area_s[k-1]) %>% log) +        # likelihood
    dtruncnorm(area_s[k-1], a = 0,
               mean = area_s_mu_prior, sd = area_s_sd_prior) %>% log + # prior
    dtruncnorm(area_s_try, a = 0,
               mean = area_s[k-1], sd = area_s_sd_jump) %>% log        # jump

  area_s_ <- ifelse(exp(lp1 - lp2) > runif(1), area_s_try, area_s[k-1])

  # a and b using metropolis
  area_a_try <- rnorm(1, area_a[k-1], area_a_sd_jump)
  area_mu_try <- area_a_try + area_b[k-1] * steps[, k-1]
  lp1 <-
    sum(dnorm(area_all, area_mu_try, area_s_, log = T)) +          # likelihood
    dnorm(area_a_try, area_a_mu_prior, area_a_sd_prior, log = T)  # prior
  lp2 <- sum(dnorm(area_all, area_mu, area_s_, log = T)) +              # likelihood
    dnorm(area_a[k-1], area_a_mu_prior, area_a_sd_prior, log = T) # prior
  area_a_ <- ifelse(exp(lp1 - lp2) > runif(1), area_a_try, area_a[k-1])
  # update area_mu with new a
  area_mu <- area_a_ + area_b[k-1] * steps[, k-1]

  area_b_try <- rnorm(1, area_b[k-1], area_b_sd_jump)
  area_mu_try <- area_a_ + area_b_try * steps[, k-1]
  lp1 <- sum(dnorm(area_all, area_mu_try, area_s_, log = T)) +          # likelihood
    dnorm(area_b_try, area_b_mu_prior, area_b_sd_prior, log = T)  # prior
  lp2 <- sum(dnorm(area_all, area_mu, area_s_, log = T)) +              # likelihood
    dnorm(area_b[k-1], area_b_mu_prior, area_b_sd_prior, log = T) # prior
  area_b_ <- ifelse(exp(lp1 - lp2) > runif(1), area_b_try, area_b[k-1])
  # update area_mu with new a and b
  area_mu <- area_a_ + area_b_ * steps[, k-1]

  #### Update parameters for steps ~ fwi regression

  # sample steps_s2 from the conditional distribution (gibbs step)
  # (compute steps_fitted, from the fwi, to get error terms)
  steps_mu <- steps_a[k-1] + steps_b[k-1] * fwi_all
  q_tmp <- J / 2 + q
  r_tmp <- 1 / (sum((steps[, k-1] - steps_mu) ^ 2) / 2 + 1 / r) # prior is inv-gamma
  steps_s_ <- sqrt(1 / rgamma(1, q_tmp, 1 / r_tmp)) # r = 1/b = scale

  # a and b using metropolis
  steps_a_try <- rnorm(1, steps_a[k-1], steps_a_sd_jump)
  steps_mu_try <- steps_a_try + steps_b[k-1] * fwi_all
  lp1 <- sum(dnorm(steps[, k-1], steps_mu_try, steps_s_, log = T)) + # likelihood of random effects
         dnorm(steps_a_try, steps_a_mu_prior, steps_a_sd_prior, log = T) # prior
  lp2 <- sum(dnorm(steps[, k-1], steps_mu, steps_s_, log = T)) + # likelihood of random effects
         dnorm(steps_a[k-1], steps_a_mu_prior, steps_a_sd_prior, log = T)   # prior
  steps_a_ <- ifelse(exp(lp1 - lp2) > runif(1), steps_a_try, steps_a[k-1])
  # update steps_fitted with new a
  steps_mu <- steps_a_ + steps_b[k-1] * fwi_all

  steps_b_try <- rnorm(1, steps_b[k-1], steps_b_sd_jump)
  steps_mu_try <- steps_a_ + steps_b_try * fwi_all
  lp1 <- sum(dnorm(steps[, k-1], steps_mu_try, steps_s_, log = T)) + # likelihood of random effects
         dnorm(steps_b_try, steps_b_mu_prior, steps_b_sd_prior, log = T) # prior
  lp2 <- sum(dnorm(steps[, k-1], steps_mu, steps_s_, log = T)) + # likelihood of random effects
         dnorm(steps_b[k-1], steps_b_mu_prior, steps_b_sd_prior, log = T)   # prior
  steps_b_ <- ifelse(exp(lp1 - lp2) > runif(1), steps_b_try, steps_b[k-1])
  # update steps_mu with new a and b
  steps_mu <- steps_a_ + steps_b_ * fwi_all

  #### Update steps (random effects)

  # first, steps used for fire spread

  # proposals sampled from stage 1 empirical density.
  # log_posterior = area_likelihood + hierarchical_prior - stage1_prior
  steps_try1 <- steps_rng()

  # compute area_mu for the proposed steps
  area_mu_try1 <- area_a_ + area_b_ * steps_try1

  lp1 <-
    dtruncnorm(area_full[ids1], a = areaL,
                mean = area_mu_try1, sd = area_s_) %>% log +
    dnorm(steps_try1, steps_mu[ids1], steps_s_, log = TRUE) -
    steps_log_prior_d(steps_try1)

  lp2 <-
    dtruncnorm(area_full[ids1], a = areaL,
               mean = area_mu[ids1], sd = area_s_) %>% log +
    dnorm(steps[ids1, k-1], steps_mu[ids1], steps_s_, log = TRUE) -
    steps_log_prior_d(steps[ids1, k-1])

  # lp1
  # lp2
  # Inf-Inf
  #
  # print(exp(lp1 - lp2))
  keep_idx <- exp(lp1 - lp2) > runif(J1)
  # keep_idx[is.na(keep_idx)] <- FALSE # override -Inf prior
  steps1_ <- steps[ids1, k-1]
  steps1_[keep_idx] <- steps_try1[keep_idx]

  # second, steps not used for fire spread
  steps_try2 <- rnorm(J2, steps[ids2, k-1], steps_sd_jump)

  # log_posterior = area_likelihood + hierarchical_prior
  area_mu_try2 <- area_a_ + area_b_ * steps_try2

  lp1 <-
    dtruncnorm(area_full[ids2], a = areaL,
               mean = area_mu_try2, sd = area_s_) %>% log +
    dnorm(steps_try2, steps_mu[ids2], steps_s_, log = TRUE)

  lp2 <-
    dtruncnorm(area_full[ids2], a = areaL,
               mean = area_mu[ids2], sd = area_s_) %>% log +
    dnorm(steps[ids2, k-1], steps_mu[ids2], steps_s_, log = TRUE)

  keep_idx <- exp(lp1 - lp2) > runif(J2)
  steps2_ <- steps[ids2, k-1]
  steps2_[keep_idx] <- steps_try2[keep_idx]

  # save samples
  steps[, k] <- c(steps1_, steps2_)
  steps_a[k] <- steps_a_
  steps_b[k] <- steps_b_
  steps_s[k] <- steps_s_

  area_a[k] <- area_a_
  area_b[k] <- area_b_
  area_s[k] <- area_s_
}

par(mfrow = c(2, 3))
plot(steps_a, type = "l")
plot(steps_b, type = "l")
plot(steps_s, type = "l")
plot(area_a, type = "l")
plot(area_b, type = "l")
plot(area_s, type = "l")
par(mfrow = c(1, 1))

par(mfrow = c(2, 3))
plot(steps[1, ], type = "l")
plot(steps[2, ], type = "l")
plot(steps[3, ], type = "l")
plot(steps[J1+1, ], type = "l")
plot(steps[J1+2, ], type = "l")
plot(steps[J1+3, ], type = "l")
par(mfrow = c(1, 1))


# Sample all parameters ---------------------------------------------------

# nsim: number of iterations to run.
# sd_jump: list containing the proposal sd for each parameter that is updated
#   through m-h, using either normal or truncated normal in the case of sigma.
# start: list with initial values

# returns a list with 3 arrays of samples: fixef, ranef and steps_extra.
spread_mcmc <- function(nsim, sd_jump, start, jump_factor = 4) {

  #### Allocate memory to save samples, using three arrays

  # fixef:       [parname, partype (a, b, s), iter]
  # ranef:       [parname, fire_id,           iter]
  # steps_extra: [1,       fire_id,           iter]

  # b param is NA for intercepts, and fixef include the area ~ steps
  # regression.

  fixef <- array(NA, dim = c(n_coef + 1, 3, nsim),
                  dimnames = list(
                    par_names = c(par_names, "area"),
                    par_class = c("a", "b", "s"),
                    iter = 1:nsim
                  ))

  ranef <- array(NA, dim = c(n_coef, nfires_spread, nsim),
                  dimnames = list(
                    par_names = par_names,
                    fire_id = fire_ids,
                    iter = 1:nsim
                  ))

  steps_extra <- array(NA, dim = c(nfires_nonspread, nsim),
                       dimnames = list(
                         fire_id = fire_ids_nonspread,
                         iter = 1:nsim
                       ))

  #### Define proposals sd (if missing)
  if(missing(sd_jump)) {

    # intercepts' alpha and sigma are updated through gibbs
    fixef_jump <- array(NA, dim = c(n_coef + 1 - 5, 3), # no intercepts
                        dimnames = list(
                          par_names = c(par_names[-(1:5)], "area"),
                          par_class = c("a", "b", "s")
                        ))

    fixef_jump["slope", ] <- apply(par_start$slope, 2, sd) * jump_factor
    fixef_jump["wind", ] <- apply(par_start$wind, 2, sd) * jump_factor
    fixef_jump["steps", ] <- apply(par_start$steps_fwi, 2, sd) * jump_factor
    fixef_jump["area", ] <- apply(par_start$area_steps, 2, sd) * jump_factor

    steps_extra_jump <- apply(par_start$steps_area$ranef, 2, sd) * jump_factor
  } else {
    fixef_jump <- sd_jump$fixef
    steps_extra_jump <- sd_jump$steps_extra
  }

  #### Define initial values (if missing)
  if(missing(start)) {

    is <- sample(1:5000, size = 1)

    fixef_start <- rbind(
      c(par_start$wet[is, 1], NA, par_start$wet[is, 2]),
      c(par_start$subalpine[is, 1], NA, par_start$subalpine[is, 2]),
      c(par_start$dry[is, 1], NA, par_start$dry[is, 2]),
      c(par_start$shrubland[is, 1], NA, par_start$shrubland[is, 2]),
      c(par_start$grassland[is, 1], NA, par_start$grassland[is, 2]),

      par_start$slope[is, ],
      par_start$wind[is, ],
      par_start$steps_fwi[is, ],
      par_start$area_steps[is, ]
    )
    colnames(fixef_start) <- c("a", "b", "s")
    rownames(fixef_start) <- c(par_names, "area")

    ranef_start <- ranef_rng(1:J)
    steps_extra_start <- par_start$steps_area$ranef[is, ]
  } else {
    fixef_start <- start$fixef
    ranef_start <- start$ranef
    steps_extra_start <- start$steps_extra
  }

  #### Initialize chain
  fixef[, , 1] <- fixef_start
  ranef[, , 1] <- ranef_start
  steps_extra[, 1] <- steps_extra_start

  #### Transient matrices to store intermediate results
  mu_fitted_keep <- array(NA, dim = c(nfires_spread, 2),
                          dimnames = list(
                            fire_id = fire_ids,
                            param = c("slope", "wind")
                          ))

  #### Run MCMC

  for(k in 2:nsim) {

    ### Fixed effects:
    # intercepts (a, s) updated through gibbs;
    # slope, wind, steps, use mh for a and b, gibbs for s;
    # area updated with mh because it uses truncated-normal likelihood.

    ## Intercepts
    for(v in 1:5) {
      # # prior parameters
      # prior_s_q = 1       # for sigma
      # prior_s_r = 1000
      # prior_a_mu = 0    # for mu
      # prior_a_sd = 5

      # sample s from the conditional distribution
      q_tmp <- J / 2 + q
      r_tmp <- 1 / (sum((ranef[v, , k-1] - fixef[v, "a"]) ^ 2) / 2 + 1 / r)
      s2 <- sqrt(1 / rgamma(1, q_tmp, 1 / r_tmp)) # r = 1/beta = scale
      fixef[v, "s", k] <- sqrt(s2)

      # sample mu (alpha) from the conditional distribution
      tmp_var <- 1 / (J / s2 + 1 / prior_a_sd^2)
      tmp_mn <- tmp.var * (sum(ranef[v, , k-1]) / s2 + mu_0 / prior_a_sd^2)
      fixef[v, "a", k] = rnorm(1, tmp_mn, sqrt(tmp_var))
    }
    # Acá arriba, usar matriz con previas


    ## Slope and wind
    for(v in c("slope", "wind")) {
      # sample s from the conditional distribution
      # (compute mu_fitted, from the fwi, to get error terms)
      mu_fitted <- fixef[v, "a", k-1] + fixef[v, "b", k-1] * fwi_sub
      q_tmp <- J / 2 + q
      r_tmp <- 1 / (sum((ranef[v, , k-1] - mu_fitted) ^ 2) / 2 + 1 / r)
      s2 <- sqrt(1 / rgamma(1, q_tmp, 1 / r_tmp)) # r = 1/beta = scale
      fixef[v, "s", k] <- sqrt(s2)

      # a and b using metropolis
      a_try <- rnorm(1, fixef[v, "a", k-1], fixef_jump[v, "a"])
      mu_try <- a_try + fixef[v, "b", k-1] * fwi_sub
      lp_new <-
        sum(dnorm(ranef[v, , k-1], mu_try, fixef[v, "s", k], log = T)) + # likelihood of random effects
        dnorm(a_try,
              fixef_prior[v, "a", "mu"], fixef_prior[v, "a", "sd"], log = T) # prior
      lp_old <-
        sum(dnorm(ranef[v, , k-1], mu_fitted, fixef[v, "s", k], log = T)) + # likelihood of random effects
        dnorm(fixef[v, "a", k-1],
              fixef_prior[v, "a", "mu"], fixef_prior[v, "a", "sd"], log = T) # prior
      fixef[v, "a", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                                 a_try, fixef[v, "a", k-1])
      # update mu with updated a
      mu_fitted <- fixef[v, "a", k] + fixef[v, "b", k-1] * fwi_sub

      b_try <- rnorm(1, fixef[v, "b", k-1], fixef_jump[v, "b"])
      mu_try <- fixef[v, "a", k] + b_try * fwi_sub
      lp_new <-
        sum(dnorm(ranef[v, , k-1], mu_try, fixef[v, "s", k], log = T)) + # likelihood of random effects
        dnorm(b_try,
              fixef_prior[v, "b", "mu"], fixef_prior[v, "b", "sd"], log = T) # prior
      lp_old <-
        sum(dnorm(ranef[v, , k-1], mu_fitted, fixef[v, "s", k], log = T)) + # likelihood of random effects
        dnorm(fixef[v, "b", k-1],
              fixef_prior[v, "b", "mu"], fixef_prior[v, "b", "sd"], log = T) # prior
      fixef[v, "b", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                                 b_try, fixef[v, "b", k-1])
      # update mu with updated b, and save (will be used to update random effects)
      mu_fitted_keep[, v] <- fixef[v, "a", k] + fixef[v, "b", k] * fwi_sub
    }

    ## Area parameters (~ steps)
    # As the likelihood is truncated-normal, the conditional distribution for
    # area_s (sigma) is unknown. Use M-H

    v <- "area" # ease subsetting

    steps_all <- c(ranef["steps", , k-1], steps_exta[, k-1])
    mu_fitted <- fixef[v, "a", k-1] + fixef[v, "b", k-1] * steps_all
    area_s_try <- rtruncnorm(1, a = 0, mean = fixef[v, "s", k-1],
                             sd = fixef_jump[v, "s"])

    lp_new <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_fitted, sd = area_s_try) %>% log) +     # likelihood
      dtruncnorm(area_s_try, a = 0,
                 mean = fixef_prior[v, "b", "mu"],
                 sd = fixef_prior[v, "s", "sd"]) %>% log +           # prior
      dtruncnorm(fixef[v, "s", k-1], a = 0,
                 mean = area_s_try, sd = fixef_jump[v, "s"]) %>% log # jump
    # this terms takes into account the asymmetric proposal (truncnorm)

    lp_old <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_fitted, sd = fixef[v, "s", k-1]) %>% log) +     # likelihood
      dtruncnorm(fixef[v, "s", k-1], a = 0,
                 mean = fixef_prior[v, "b", "mu"],
                 sd = fixef_prior[v, "s", "sd"]) %>% log +                   # prior
      dtruncnorm(area_s_try, a = 0,
                 mean = fixef[v, "s", k-1], sd = fixef_jump[v, "s"]) %>% log # jump

    fixef[v, "s", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               area_s_try, fixef[v, "s", k-1])

    # a and b using metropolis
    a_try <- rnorm(1, fixef[v, "a", k-1], fixef_jump[v, "a"])
    mu_try <- a_try + fixef[v, "b", k-1] * steps_all
    lp_new <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_try, sd = fixef[v, "s", k]) %>% log) +      # likelihood
      dnorm(a_try,
            fixef_prior[v, "a", "mu"], fixef_prior[v, "a", "sd"], log = T) # prior
    lp_old <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_fitted, sd = fixef[v, "s", k]) %>% log) +   # likelihood
      dnorm(fixef[v, "a", k-1],
            fixef_prior[v, "a", "mu"], fixef_prior[v, "a", "sd"], log = T) # prior
    fixef[v, "a", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               a_try, fixef[v, "a", k-1])
    # update mu with updated a
    mu_fitted <- fixef[v, "a", k] + fixef[v, "b", k-1] * steps_all

    b_try <- rnorm(1, fixef[v, "b", k-1], fixef_jump[v, "b"])
    mu_try <- fixef[v, "a", k] + b_try * steps_all
    lp_new <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_try, sd = fixef[v, "s", k]) %>% log) +      # likelihood
      dnorm(b_try,
            fixef_prior[v, "b", "mu"], fixef_prior[v, "b", "sd"], log = T) # prior
    lp_old <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_fitted, sd = fixef[v, "s", k]) %>% log) +   # likelihood
      dnorm(fixef[v, "b", k-1],
            fixef_prior[v, "b", "mu"], fixef_prior[v, "b", "sd"], log = T) # prior
    fixef[v, "b", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               b_try, fixef[v, "b", k-1])

    # compute mu to compute the likelihood when updating steps
    area_mu <- fixef[v, "a", k] + fixef[v, "b", k] * steps_all

    ## Steps parameters (~ fwi)

    v <- "steps" # ease subsetting
    mu_fitted <- fixef[v, "a", k-1] + fixef[v, "b", k-1] * fwi_all
    q_tmp <- J / 2 + q
    r_tmp <- 1 / (sum((ranef[v, , k-1] - mu_fitted) ^ 2) / 2 + 1 / r) # prior is inv-gamma
    fixef[v, "s", k] <- sqrt(1 / rgamma(1, q_tmp, 1 / r_tmp)) # r = 1/b = scale

    # a and b using metropolis
    a_try <- rnorm(1, fixef[v, "a", k-1], fixef_jump[v, "a"])
    mu_try <- a_try + fixef[v, "b", k-1] * fwi_all
    lp_new <-
      sum(dnorm(steps_all, mu_try, fixef[v, "s", k], log = T)) +           # likelihood for ranef
      dnorm(a_try,
            fixef_prior[v, "a", "mu"], fixef_prior[v, "a", "sd"], log = T) # prior
    lp_old <-
      sum(dnorm(steps_all, mu_fitted, fixef[v, "s", k], log = T)) +        # likelihood for ranef
      dnorm(fixef[v, "a", k-1],
            fixef_prior[v, "a", "mu"], fixef_prior[v, "a", "sd"], log = T) # prior
    fixef[v, "a", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               a_try, fixef[v, "a", k-1])
    # update mu with updated a
    mu_fitted <- fixef[v, "a", k] + fixef[v, "b", k-1] * fwi_all

    b_try <- rnorm(1, fixef[v, "b", k-1], fixef_jump[v, "b"])
    mu_try <- fixef[v, "a", k] + b_try * fwi_all
    lp_new <-
      sum(dnorm(steps_all, mu_try, fixef[v, "s", k], log = T)) +      # likelihood
      dnorm(b_try,
            fixef_prior[v, "b", "mu"], fixef_prior[v, "b", "sd"], log = T) # prior
    lp_old <-
      sum(dnorm(steps_all, mu_fitted, fixef[v, "s", k], log = T)) +   # likelihood
      dnorm(fixef[v, "b", k-1],
            fixef_prior[v, "b", "mu"], fixef_prior[v, "b", "sd"], log = T) # prior
    fixef[v, "b", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               b_try, fixef[v, "b", k-1])

    # keep steps mu to update random effects
    steps_fitted <- fixef[v, "a", k] + fixef[v, "b", k] * fwi_all

    ## Random effects (Lunn method)
    # simulate new random effects (proposals from stage1 kdes)
    for(f in 1:nfires_spread) {
      # simulate new random effect, handling possible Inf in the log-steps prior
      log_prior_steps <- -Inf
      while(!is.finite(log_prior_steps)) {
        ranef_try <- rkdevine(1, ranef_kde[[f]])
        names(ranef_try) <- par_names
        log_prior_steps <- expunif_lpdf(ranef_try[n_coef],
                                        l = steps_bounds[f, 1],
                                        u = steps_bounds[f, 2])
                           #steps_prior_lpdf(ranef_try[n_coef], f)
      }

      # compute area mu for the proposed steps
      area_mu_try <- fixef["area", "a", k] +
                     fixef["area", "b", k] * ranef_try["steps"]

      lp_new <-
        # steps log-posterior: area likelihood + hierarchical prior
        dtruncnorm(area_full[f], a = areaL,
                   mean = area_mu_try, sd = fixef["area", "s", k]) %>% log +
        dnorm(ranef_try["steps"], steps_fitted[f], fixef["steps", "s", k], log = TRUE) +

        # slope and wind log-prior (hierarchical), based on varying-mu
        dnorm(ranef_try["slope"],
              mu_fitted_keep[f, "slope"], fixef["slope", "s", k], log = TRUE) +
        dnorm(ranef_try["wind"],
              mu_fitted_keep[f, "wind"], fixef["wind", "s", k], log = TRUE) +

        # intercepts hierarchical prior, based on alpha (a)
        sum(dnorm(ranef_try[1:5],
                  fixef[1:5, "a", k], fixef[1:5, "s", k], log = TRUE)) -

        # substract the stage1 prior
        (
          log_prior_steps +
          dlogis(ranef_try[1:(n_coef-1)], log = T)
        )

      lp_old <-
        # steps log-posterior: area likelihood + hierarchical prior
        dtruncnorm(area_full[f], a = areaL,
                   mean = area_mu[f], sd = fixef["area", "s", k]) %>% log +
        dnorm(ranef["steps", f, k-1],
              steps_fitted[f], fixef["steps", "s", k], log = TRUE) +

        # slope and wind log-prior (hierarchical), based on varying-mu
        dnorm(ranef["slope", f, k-1],
              mu_fitted_keep[f, "slope"], fixef["slope", "s", k], log = TRUE) +
        dnorm(ranef["wind", f, k-1],
              mu_fitted_keep[f, "wind"], fixef["wind", "s", k], log = TRUE) +

        # intercepts hierarchical prior, based on alpha (a)
        sum(dnorm(ranef[1:5, f, k-1],
                  fixef[1:5, "a", k], fixef[1:5, "s", k], log = TRUE)) -

        # substract the stage1
        (
          steps_prior_lpdf(ranef["steps", f, k-1], f) +
          dlogis(ranef[1:(n_coef-1), f, k-1], log = T)
        )

      ranef[, f, k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                              ranef_try, ranef[, f, k-1])
    }

    ## Steps random effects not involved in fire spread
    steps_try <- rnorm(J2, steps_extra[, k-1], steps_extra_jump)

    # log_posterior = area_likelihood + hierarchical_prior
    area_mu_try2 <- fixef["area", "a", k] + fixef["area", "b", k] * steps_try

    lp_new <-
      dtruncnorm(area_full[ids2], a = areaL,
                 mean = area_mu_try2, sd = fixef["area", "s", k]) %>% log +
      dnorm(steps_try,
            steps_fitted[ids2], fixef["steps", "s", k], log = TRUE)

    lp_old <-
      dtruncnorm(area_full[ids2], a = areaL,
                 mean = area_mu[ids2], sd = fixef["area", "s", k]) %>% log +
      dnorm(steps_extra[, k-1],
            steps_fitted[ids2], fixef["steps", "s", k], log = TRUE)

    keep_idx <- exp(lp1 - lp2) > runif(J2)
    steps_extra[keep_idx, k] <- steps_try[keep_idx]
    steps_extra[!keep_idx, k] <- steps_extra[!keep_idx, k-1]
  }

  #### Merge samples
  out <- list(
    fixef = fixef,
    ranef = ranef,
    steps_extra = steps_extra
  )

  return(out)
}


# returns a list with 3 arrays of samples: fixef, ranef and steps_extra.
# It also returns the acceptance ratio of a run
spread_mcmc_adapt <- function(nsim, sd_jump, start, jump_factor = 4) {

  #### Allocate memory to save samples, using three arrays

  # fixef:       [parname, partype (a, b, s), iter]
  # ranef:       [parname, fire_id,           iter]
  # steps_extra: [1,       fire_id,           iter]

  # b param is NA for intercepts, and fixef include the area ~ steps
  # regression.

  fixef <- array(NA, dim = c(n_coef + 1, 3, nsim),
                 dimnames = list(
                   par_names = c(par_names, "area"),
                   par_class = c("a", "b", "s"),
                   iter = 1:nsim
                 ))

  ranef <- array(NA, dim = c(n_coef, nfires_spread, nsim),
                 dimnames = list(
                   par_names = par_names,
                   fire_id = fire_ids,
                   iter = 1:nsim
                 ))

  steps_extra <- array(NA, dim = c(nfires_nonspread, nsim),
                       dimnames = list(
                         fire_id = fire_ids_nonspread,
                         iter = 1:nsim
                       ))

  # Create the equivalent arrays to compute the acceptance rate
  fixef_accept <- fixef[, , 1]
  fixef_accept[, ] <- 0

  ranef_accept <- ranef[, , 1]
  ranef_accept[, ] <- 0

  steps_extra_accept <- steps_extra[1, ]
  steps_extra_accept[] <- 0

  #### Define proposals sd (if missing)
  if(missing(sd_jump)) {

    # intercepts' alpha and sigma are updated through gibbs
    fixef_jump <- array(NA, dim = c(n_coef + 1 - 5, 3), # no intercepts
                        dimnames = list(
                          par_names = c(par_names[-(1:5)], "area"),
                          par_class = c("a", "b", "s")
                        ))

    fixef_jump["slope", ] <- apply(par_start$slope, 2, sd) * jump_factor
    fixef_jump["wind", ] <- apply(par_start$wind, 2, sd) * jump_factor
    fixef_jump["steps", ] <- apply(par_start$steps_fwi, 2, sd) * jump_factor
    fixef_jump["area", ] <- apply(par_start$area_steps, 2, sd) * jump_factor

    steps_extra_jump <- apply(par_start$steps_area$ranef, 2, sd) * jump_factor
  } else {
    fixef_jump <- sd_jump$fixef
    steps_extra_jump <- sd_jump$steps_extra
  }

  #### Define initial values (if missing)
  if(missing(start)) {

    is <- sample(1:5000, size = 1)

    fixef_start <- rbind(
      c(par_start$wet[is, 1], NA, par_start$wet[is, 2]),
      c(par_start$subalpine[is, 1], NA, par_start$subalpine[is, 2]),
      c(par_start$dry[is, 1], NA, par_start$dry[is, 2]),
      c(par_start$shrubland[is, 1], NA, par_start$shrubland[is, 2]),
      c(par_start$grassland[is, 1], NA, par_start$grassland[is, 2]),

      par_start$slope[is, ],
      par_start$wind[is, ],
      par_start$steps_fwi[is, ],
      par_start$area_steps[is, ]
    )
    colnames(fixef_start) <- c("a", "b", "s")
    rownames(fixef_start) <- c(par_names, "area")

    ranef_start <- ranef_rng()
    steps_extra_start <- par_start$steps_area$ranef[is, ]
  } else {
    fixef_start <- start$fixef
    ranef_start <- start$ranef
    steps_extra_start <- start$steps_extra
  }

  #### Initialize chain
  fixef[, , 1] <- fixef_start
  ranef[, , 1] <- ranef_start
  steps_extra[, 1] <- steps_extra_start

  #### Transient matrices to store intermediate results
  mu_fitted_keep <- array(NA, dim = c(nfires_spread, 2),
                          dimnames = list(
                            fire_id = fire_ids,
                            param = c("slope", "wind")
                          ))

  #### Run MCMC, counting acceptance

  for(k in 2:nsim) {

    ### Fixed effects:
    # intercepts (a, s) updated through gibbs;
    # slope, wind, steps, use mh for a and b, gibbs for s;
    # area updated with mh because it uses truncated-normal likelihood.

    ## Intercepts
    for(v in 1:5) {
      # # prior parameters
      # prior_s_q = 1       # for sigma
      # prior_s_r = 1000
      # prior_a_mu = 0    # for mu
      # prior_a_sd = 5

      # sample s from the conditional distribution
      q_tmp <- J / 2 + q
      r_tmp <- 1 / (sum((ranef[v, , k-1] - fixef[v, "a"]) ^ 2) / 2 + 1 / r)
      s2 <- sqrt(1 / rgamma(1, q_tmp, 1 / r_tmp)) # r = 1/beta = scale
      fixef[v, "s", k] <- sqrt(s2)

      # sample mu (alpha) from the conditional distribution
      tmp_var <- 1 / (J / s2 + 1 / prior_a_sd^2)
      tmp_mn <- tmp.var * (sum(ranef[v, , k-1]) / s2 + mu_0 / prior_a_sd^2)
      fixef[v, "a", k] = rnorm(1, tmp_mn, sqrt(tmp_var))
    }
    # Acá arriba, usar matriz con previas


    ## Slope and wind
    for(v in c("slope", "wind")) {
      # sample s from the conditional distribution
      # (compute mu_fitted, from the fwi, to get error terms)
      mu_fitted <- fixef[v, "a", k-1] + fixef[v, "b", k-1] * fwi_sub
      q_tmp <- J / 2 + q
      r_tmp <- 1 / (sum((ranef[v, , k-1] - mu_fitted) ^ 2) / 2 + 1 / r)
      s2 <- sqrt(1 / rgamma(1, q_tmp, 1 / r_tmp)) # r = 1/beta = scale
      fixef[v, "s", k] <- sqrt(s2)

      # a and b using metropolis
      a_try <- rnorm(1, fixef[v, "a", k-1], fixef_jump[v, "a"])
      mu_try <- a_try + fixef[v, "b", k-1] * fwi_sub
      lp_new <-
        sum(dnorm(ranef[v, , k-1], mu_try, fixef[v, "s", k], log = T)) + # likelihood of random effects
        dnorm(a_try,
              fixef_prior[v, "a", "mu"], fixef_prior[v, "a", "sd"], log = T) # prior
      lp_old <-
        sum(dnorm(ranef[v, , k-1], mu_fitted, fixef[v, "s", k], log = T)) + # likelihood of random effects
        dnorm(fixef[v, "a", k-1],
              fixef_prior[v, "a", "mu"], fixef_prior[v, "a", "sd"], log = T) # prior
      fixef[v, "a", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                                 a_try, fixef[v, "a", k-1])
      # update mu with updated a
      mu_fitted <- fixef[v, "a", k] + fixef[v, "b", k-1] * fwi_sub

      b_try <- rnorm(1, fixef[v, "b", k-1], fixef_jump[v, "b"])
      mu_try <- fixef[v, "a", k] + b_try * fwi_sub
      lp_new <-
        sum(dnorm(ranef[v, , k-1], mu_try, fixef[v, "s", k], log = T)) + # likelihood of random effects
        dnorm(b_try,
              fixef_prior[v, "b", "mu"], fixef_prior[v, "b", "sd"], log = T) # prior
      lp_old <-
        sum(dnorm(ranef[v, , k-1], mu_fitted, fixef[v, "s", k], log = T)) + # likelihood of random effects
        dnorm(fixef[v, "b", k-1],
              fixef_prior[v, "b", "mu"], fixef_prior[v, "b", "sd"], log = T) # prior
      fixef[v, "b", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                                 b_try, fixef[v, "b", k-1])
      # update mu with updated b, and save (will be used to update random effects)
      mu_fitted_keep[, v] <- fixef[v, "a", k] + fixef[v, "b", k] * fwi_sub
    }

    ## Area parameters (~ steps)
    # As the likelihood is truncated-normal, the conditional distribution for
    # area_s (sigma) is unknown. Use M-H

    v <- "area" # ease subsetting

    steps_all <- c(ranef["steps", , k-1], steps_exta[, k-1])
    mu_fitted <- fixef[v, "a", k-1] + fixef[v, "b", k-1] * steps_all
    area_s_try <- rtruncnorm(1, a = 0, mean = fixef[v, "s", k-1],
                             sd = fixef_jump[v, "s"])

    lp_new <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_fitted, sd = area_s_try) %>% log) +     # likelihood
      dtruncnorm(area_s_try, a = 0,
                 mean = fixef_prior[v, "b", "mu"],
                 sd = fixef_prior[v, "s", "sd"]) %>% log +           # prior
      dtruncnorm(fixef[v, "s", k-1], a = 0,
                 mean = area_s_try, sd = fixef_jump[v, "s"]) %>% log # jump
    # this terms takes into account the asymmetric proposal (truncnorm)

    lp_old <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_fitted, sd = fixef[v, "s", k-1]) %>% log) +     # likelihood
      dtruncnorm(fixef[v, "s", k-1], a = 0,
                 mean = fixef_prior[v, "b", "mu"],
                 sd = fixef_prior[v, "s", "sd"]) %>% log +                   # prior
      dtruncnorm(area_s_try, a = 0,
                 mean = fixef[v, "s", k-1], sd = fixef_jump[v, "s"]) %>% log # jump

    fixef[v, "s", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               area_s_try, fixef[v, "s", k-1])

    # a and b using metropolis
    a_try <- rnorm(1, fixef[v, "a", k-1], fixef_jump[v, "a"])
    mu_try <- a_try + fixef[v, "b", k-1] * steps_all
    lp_new <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_try, sd = fixef[v, "s", k]) %>% log) +      # likelihood
      dnorm(a_try,
            fixef_prior[v, "a", "mu"], fixef_prior[v, "a", "sd"], log = T) # prior
    lp_old <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_fitted, sd = fixef[v, "s", k]) %>% log) +   # likelihood
      dnorm(fixef[v, "a", k-1],
            fixef_prior[v, "a", "mu"], fixef_prior[v, "a", "sd"], log = T) # prior
    fixef[v, "a", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               a_try, fixef[v, "a", k-1])
    # update mu with updated a
    mu_fitted <- fixef[v, "a", k] + fixef[v, "b", k-1] * steps_all

    b_try <- rnorm(1, fixef[v, "b", k-1], fixef_jump[v, "b"])
    mu_try <- fixef[v, "a", k] + b_try * steps_all
    lp_new <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_try, sd = fixef[v, "s", k]) %>% log) +      # likelihood
      dnorm(b_try,
            fixef_prior[v, "b", "mu"], fixef_prior[v, "b", "sd"], log = T) # prior
    lp_old <-
      sum(dtruncnorm(area_all, a = areaL,
                     mean = mu_fitted, sd = fixef[v, "s", k]) %>% log) +   # likelihood
      dnorm(fixef[v, "b", k-1],
            fixef_prior[v, "b", "mu"], fixef_prior[v, "b", "sd"], log = T) # prior
    fixef[v, "b", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               b_try, fixef[v, "b", k-1])

    # compute mu to compute the likelihood when updating steps
    area_mu <- fixef[v, "a", k] + fixef[v, "b", k] * steps_all

    ## Steps parameters (~ fwi)

    v <- "steps" # ease subsetting
    mu_fitted <- fixef[v, "a", k-1] + fixef[v, "b", k-1] * fwi_all
    q_tmp <- J / 2 + q
    r_tmp <- 1 / (sum((ranef[v, , k-1] - mu_fitted) ^ 2) / 2 + 1 / r) # prior is inv-gamma
    fixef[v, "s", k] <- sqrt(1 / rgamma(1, q_tmp, 1 / r_tmp)) # r = 1/b = scale

    # a and b using metropolis
    a_try <- rnorm(1, fixef[v, "a", k-1], fixef_jump[v, "a"])
    mu_try <- a_try + fixef[v, "b", k-1] * fwi_all
    lp_new <-
      sum(dnorm(steps_all, mu_try, fixef[v, "s", k], log = T)) +           # likelihood for ranef
      dnorm(a_try,
            fixef_prior[v, "a", "mu"], fixef_prior[v, "a", "sd"], log = T) # prior
    lp_old <-
      sum(dnorm(steps_all, mu_fitted, fixef[v, "s", k], log = T)) +        # likelihood for ranef
      dnorm(fixef[v, "a", k-1],
            fixef_prior[v, "a", "mu"], fixef_prior[v, "a", "sd"], log = T) # prior
    fixef[v, "a", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               a_try, fixef[v, "a", k-1])
    # update mu with updated a
    mu_fitted <- fixef[v, "a", k] + fixef[v, "b", k-1] * fwi_all

    b_try <- rnorm(1, fixef[v, "b", k-1], fixef_jump[v, "b"])
    mu_try <- fixef[v, "a", k] + b_try * fwi_all
    lp_new <-
      sum(dnorm(steps_all, mu_try, fixef[v, "s", k], log = T)) +      # likelihood
      dnorm(b_try,
            fixef_prior[v, "b", "mu"], fixef_prior[v, "b", "sd"], log = T) # prior
    lp_old <-
      sum(dnorm(steps_all, mu_fitted, fixef[v, "s", k], log = T)) +   # likelihood
      dnorm(fixef[v, "b", k-1],
            fixef_prior[v, "b", "mu"], fixef_prior[v, "b", "sd"], log = T) # prior
    fixef[v, "b", k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                               b_try, fixef[v, "b", k-1])

    # keep steps mu to update random effects
    steps_fitted <- fixef[v, "a", k] + fixef[v, "b", k] * fwi_all

    ## Random effects (Lunn method)
    # simulate new random effects (proposals from stage1 kdes)
    for(f in 1:nfires_spread) {
      # simulate new random effect, handling possible Inf in the log-steps prior
      log_prior_steps <- Inf
      while(!is.finite(log_prior_steps)) {
        ranef_try <- rkdevine(1, ranef_kde[[f]])
        names(ranef_try) <- par_names
        log_prior_steps <- steps_prior_lpdf(ranef_try[n_coef], f)
      }

      # compute area mu for the proposed steps
      area_mu_try <- fixef["area", "a", k] +
        fixef["area", "b", k] * ranef_try["steps"]

      lp_new <-
        # steps log-posterior: area likelihood + hierarchical prior
        dtruncnorm(area_full[f], a = areaL,
                   mean = area_mu_try, sd = fixef["area", "s", k]) %>% log +
        dnorm(ranef_try["steps"], steps_fitted[f], fixef["steps", "s", k], log = TRUE) +

        # slope and wind log-prior (hierarchical), based on varying-mu
        dnorm(ranef_try["slope"],
              mu_fitted_keep[f, "slope"], fixef["slope", "s", k], log = TRUE) +
        dnorm(ranef_try["wind"],
              mu_fitted_keep[f, "wind"], fixef["wind", "s", k], log = TRUE) +

        # intercepts hierarchical prior, based on alpha (a)
        sum(dnorm(ranef_try[1:5],
                  fixef[1:5, "a", k], fixef[1:5, "s", k], log = TRUE)) -

        # substract the stage1 prior
        (
          log_prior_steps +
            dlogis(ranef_try[1:(n_coef-1)], log = T)
        )

      lp_old <-
        # steps log-posterior: area likelihood + hierarchical prior
        dtruncnorm(area_full[f], a = areaL,
                   mean = area_mu[f], sd = fixef["area", "s", k]) %>% log +
        dnorm(ranef["steps", f, k-1],
              steps_fitted[f], fixef["steps", "s", k], log = TRUE) +

        # slope and wind log-prior (hierarchical), based on varying-mu
        dnorm(ranef["slope", f, k-1],
              mu_fitted_keep[f, "slope"], fixef["slope", "s", k], log = TRUE) +
        dnorm(ranef["wind", f, k-1],
              mu_fitted_keep[f, "wind"], fixef["wind", "s", k], log = TRUE) +

        # intercepts hierarchical prior, based on alpha (a)
        sum(dnorm(ranef[1:5, f, k-1],
                  fixef[1:5, "a", k], fixef[1:5, "s", k], log = TRUE)) -

        # substract the stage1
        (
          steps_prior_lpdf(ranef["steps", f, k-1], f) +
            dlogis(ranef[1:(n_coef-1), f, k-1], log = T)
        )

      ranef[, f, k] <- ifelse(exp(lp_new - lp_old) > runif(1),
                              ranef_try, ranef[, f, k-1])
    }

    ## Steps random effects not involved in fire spread
    steps_try <- rnorm(J2, steps_extra[, k-1], steps_extra_jump)

    # log_posterior = area_likelihood + hierarchical_prior
    area_mu_try2 <- fixef["area", "a", k] + fixef["area", "b", k] * steps_try

    lp_new <-
      dtruncnorm(area_full[ids2], a = areaL,
                 mean = area_mu_try2, sd = fixef["area", "s", k]) %>% log +
      dnorm(steps_try,
            steps_fitted[ids2], fixef["steps", "s", k], log = TRUE)

    lp_old <-
      dtruncnorm(area_full[ids2], a = areaL,
                 mean = area_mu[ids2], sd = fixef["area", "s", k]) %>% log +
      dnorm(steps_extra[, k-1],
            steps_fitted[ids2], fixef["steps", "s", k], log = TRUE)

    keep_idx <- exp(lp1 - lp2) > runif(J2)
    steps_extra[keep_idx, k] <- steps_try[keep_idx]
    steps_extra[!keep_idx, k] <- steps_extra[!keep_idx, k-1]
  }

  #### Merge samples
  out <- list(
    fixef = fixef,
    ranef = ranef,
    steps_extra = steps_extra
  )

  return(out)
}

##########################################################
######  TAREAS:
# CREAR MATRIZ DE PREVIAS.
# CREAR MATRIZ DE MU_FITTED INTERNA para guardar los mu de slope, wind y steps,
#   y usarlos al ACTUALIZAR RANDOM EFFECTS.
# create FWI_SUB (wind and slope) and FWI_ALL.

# Hacer el conteo de aceptaciones para cada parámetro que requiera M-H,
# así se puede usar esta func para tunear las proposals.
