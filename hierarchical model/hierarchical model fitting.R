# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(bayesplot)
library(bayestestR) # hdi
library(posterior)
library(deeptime)
library(terra)
library(kdevine)
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



# Load posteriors ---------------------------------------------------------

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
fire_ids_nonspread <- fires_data_nonspread$fire_id
nfires_nonspread <- length(fire_ids_nonspread)

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


# log-uniform priors ------------------------------------------------------

# In the first stage, steps were assigned a flat prior, between 5 and S, with S
# varying among fires. As we aim to estimate the steps at the log scale, we
# need to know which was the prior distribution at the log scale, which is not
# a uniform:

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
    steps_log_prior_kde
  }
}

# MCMC for steps and area -----------------------------------------------

steps_kde <- lapply(samples_temp, function(x) {
  kk <- kde1d(x[, "steps"])
})

# Copying code from Hooten and Hefley (2019) book.

# setup variables
nsim <- 500
J1 <- length(kde_list)
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
steps_j = matrix(NA, J, nsim)

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

# proposal scales for a and b
steps_a_sd_jump <- sd(par_start$steps_fwi[, "alpha"]) * 2
steps_b_sd_jump <- sd(par_start$steps_fwi[, "beta"]) * 2
steps_j_sd_jump <- apply(par_start$steps_area$ranef, 2, sd) * 2

# for area
area_a_sd_jump <- sd(par_start$area_steps[, "alpha"]) * 2
area_b_sd_jump <- sd(par_start$area_steps[, "beta"]) * 2

# initial values
id_start <- sample(1:5000, 1)

steps_j[ids1, 1] = sapply(kde_list, function(x) rkde1d(1, x))
steps_j[ids2, 1] = par_start$steps_area$ranef[id_start, ]

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

for(k in 2:nsim) {
  # k = 2
  print(k)

  #### Update parameters for steps ~ fwi regression

  # sample steps_s2 from the conditional distribution
  # (compute steps_fitted, from the fwi, to get error terms)
  steps_mu <- steps_a[k-1] + steps_b[k-1] * fwi_all
  q_tmp = J / 2 + q
  r_tmp = 1 / (sum((steps_j[, k-1] - steps_mu) ^ 2) / 2 + 1 / r)
  steps_s = sqrt(1 / rgamma(1, q_tmp, 1 / r_tmp)) # r = 1/b = scale

  # a and b using metropolis
  steps_a_try = rnorm(1, steps_a[k-1], steps_a_sd_jump)
  steps_mu_try = steps_a_try + steps_b[k-1] * fwi_all
  lp1 <- sum(dnorm(steps_j[, k-1], steps_mu_try, steps_s, log = T)) + # likelihood of random effects
         dnorm(steps_a_try, steps_a_mu_prior, steps_a_sd_prior, log = T) # prior
  lp2 <- sum(dnorm(steps_j[, k-1], steps_mu, steps_s, log = T)) + # likelihood of random effects
         dnorm(steps_a[k-1], steps_a_mu_prior, steps_a_sd_prior, log = T)   # prior
  steps_a <- ifelse(exp(lp1 - lp2) > runif(1), steps_a_try, steps_a[k-1])
  # update steps_fitted with new a
  steps_mu <- steps_a + steps_b[k-1] * fwi_all

  steps_b_try = rnorm(1, steps_b[k-1], steps_b_sd_jump)
  steps_mu_try = steps_a + steps_b_try * fwi_all
  lp1 <- sum(dnorm(steps_j[, k-1], steps_mu_try, steps_s, log = T)) + # likelihood of random effects
         dnorm(steps_b_try, steps_b_mu_prior, steps_b_sd_prior, log = T) # prior
  lp2 <- sum(dnorm(steps_j[, k-1], steps_mu, steps_s, log = T)) + # likelihood of random effects
         dnorm(steps_b[k-1], steps_b_mu_prior, steps_b_sd_prior, log = T)   # prior
  steps_b <- ifelse(exp(lp1 - lp2) > runif(1), steps_b_try, steps_b[k-1])
  # update steps_mu with new a and b
  steps_mu <- steps_a + steps_b * fwi_all

  #### Update parameters for area ~ steps regression

  # sample area_s2 from the conditional distribution
  # (compute area_mu, from the steps, to get error terms)
  area_mu <- area_a[k-1] + area_b[k-1] * steps_j[, k-1]
  q_tmp = J / 2 + q
  r_tmp = 1 / (sum((area_j[, k-1] - area_mu) ^ 2) / 2 + 1 / r)
  area_s = sqrt(1 / rgamma(1, q_tmp, 1 / r_tmp)) # r = 1/b = scale

  # a and b using metropolis
  area_a_try = rnorm(1, area_a[k-1], area_a_sd_jump)
  area_mu_try = area_a_try + area_b[k-1] * steps_j[, k-1]
  lp1 <- sum(dnorm(area_all, area_mu_try, area_s, log = T)) +          # likelihood
         dnorm(area_a_try, area_a_mu_prior, area_a_sd_prior, log = T)  # prior
  lp2 <- sum(dnorm(area_all, area_mu, area_s, log = T)) +              # likelihood
         dnorm(area_a[k-1], area_a_mu_prior, area_a_sd_prior, log = T) # prior
  area_a <- ifelse(exp(lp1 - lp2) > runif(1), area_a_try, area_a[k-1])
  # update area_mu with new a
  area_mu <- area_a + area_b[k-1] * steps_j[, k-1]

  area_b_try = rnorm(1, area_b[k-1], area_b_sd_jump)
  area_mu_try = area_a + area_b_try * steps_j[, k-1]
  lp1 <- sum(dnorm(area_all, area_mu_try, area_s, log = T)) +          # likelihood
         dnorm(area_b_try, area_b_mu_prior, area_b_sd_prior, log = T)  # prior
  lp2 <- sum(dnorm(area_all, area_mu, area_s, log = T)) +              # likelihood
         dnorm(area_b[k-1], area_b_mu_prior, area_b_sd_prior, log = T) # prior
  area_b <- ifelse(exp(lp1 - lp2) > runif(1), area_b_try, area_b[k-1])
  # update area_mu with new a and b
  area_mu <- area_a + area_b * steps_j[, k-1]


  ### Update steps (random effects) seguiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiirrrrrrrrrrrrrrr

  # first, steps used for fire spread
  steps_star1 = sapply(steps_kde, function(x) rkde1d(1, x))
  # proposals sampled from stage1 empirical density.
  # (the metropolis dividends are likelihood * other's_prior)
  lp1 = dnorm(muj_star, mu_fitted, sqrt(s2), log = TRUE) +
        sapply(steps_log_prior_kde, )
  lp2 = dnorm(muj_hat[, k-1], mu_fitted, sqrt(s2), log = TRUE) +
    dlogis(muj_star, 0, 1, log = TRUE) # prior from stage 1 for proposal
  keep_idx = exp(lp1 - lp2) > runif(J)

  # save samples
  muj = muj_hat[, k-1]
  muj[keep_idx] = muj_star[keep_idx]

  muj_hat[, k] = muj
  a_hat[k] = a
  b_hat[k] = b
  s_hat[k] = sqrt(s2)
}
