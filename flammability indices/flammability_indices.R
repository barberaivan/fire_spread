# Assume known the relative differences in logit-burn-probability related
# to vegetation and topography, estimated by regional logistic regressions.

# Vegetation flammability index (VFI):
#   intercepts[v] + b_ndvi[v] * (ndvi[i] - opt[v]) ^ 2
#
# Considering random effects for the three parameters, with 5 vegetation types:
# {
#   1: wet forest,
#   2: subalpine forest,
#   3: dry forest,
#   4: shrubland,
#   5: grassland
# }

# The topographic flammability index (TFI) is defined by elevation and northing:
#   b_elev * elev[i] + b_north * north[i]

# Slope is included in the model to remove its effect,
# but will not be used to compute the TFI.

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_minimal())
library(terra)
library(rstan)
library(posterior)
library(mgcv)
library(DHARMa)
library(HDInterval)
library(loo)

# Functions ---------------------------------------------------------------

# to compute mean and 95 % CI from samples.
mean_ci <- function(x, name = "p_") {
  qs <- quantile(x, probs = c(0.025, 0.975), method = 8)
  tmp <- c(mean(x), qs) %>% unname
  names(tmp) <- paste(name, c("mean", "lower", "upper"), sep = "")
  return(tmp)
}

mean_hdi <- function(x, name = "p_", prob = 0.95) {
  hh <- hdi(x, prob)
  tmp <- c(mean(x), hh["lower"], hh["upper"])
  names(tmp) <- paste(name, c("mean", "lower", "upper"), sep = "")
  return(tmp)
}

# Prepare data -------------------------------------------------------------

# model to detrend ndvi
mdetrend <- readRDS(file.path("data", "flammability indices",
                              "ndvi_detrender_model.rds"))

v0 <- vect(file.path("data", "flammability indices",
                     "ndvi_regional_points.shp"))

dveg <- readxl::read_excel("/home/ivan/Insync/Mapa vegetación WWF - Lara et al. 1999/clases de vegetacion y equivalencias.xlsx",
                           sheet = "Sheet2")

names(dveg) <- c("vegnum1", "vegfac1", "vegnum2", "vegfac2")
dveg$vegfac2 <- factor(dveg$vegfac2,
                       levels = c("Wet forest", "Subalpine forest", "Dry forest",
                                  "Shrubland", "Grassland", "Non burnable"))

names(v0)[names(v0) == "veg"] <- "vegnum1"
names(v0)[names(v0) == "elev"] <- "elevation"

# remove the non-burnable
data <- as.data.frame(v0)
data <- left_join(data, dveg, by = "vegnum1")
data <- data[data$vegnum1 < 7, ] # urban, water, high andean, and ice and snow are out

# compute slope-weighted northing
# (see code <northing importance function.R>)
data$northing <- cos(data$aspect * pi / 180) * plogis(-5 + 0.35 * data$slope)

# Data to fit model -------------------------------------------------------

rows_fit <- complete.cases(data[, c("burned",
                                    "vegfac2",
                                    "vegnum2",
                                    "ndvi_dyn",
                                    "year_dyn",
                                    "ndvi_22",
                                    "elevation",
                                    "slope",
                                    "northing")])
dfit <- data[rows_fit, ]; nrow(dfit); nrow(data)# OK

veg_levels <- unique(dveg$vegfac2[dveg$vegnum2 < 6] %>% as.character)

dfit$vegfac2 <- factor(dfit$vegfac2 %>% as.character,
                       levels = veg_levels)

dfit <- dfit[complete.cases(dfit), ]

veg <- model.matrix(~ vegfac2 - 1, dfit)
veg_num <- dfit$vegnum2
V <- ncol(veg)
veg_classes <- levels(dfit$vegfac2)

dfit$elev_z <- scale(dfit$elevation) %>% as.numeric
dfit$slope_z <- scale(dfit$slope) %>% as.numeric
dfit$north_z <- scale(dfit$northing) %>% as.numeric

# detrend ndvi
dpred_ndvi <- data.frame(
  ndvi_dyn_logit = qlogis((dfit$ndvi_dyn + 1) / 2),
  ndvi01_22 = (dfit$ndvi_22 + 1) / 2,
  year = dfit$year
)

dpred_ndvi$diff_logit <- predict(mdetrend, dpred_ndvi, se.fit = F)
dpred_ndvi$ndvi_dt <- plogis(dpred_ndvi$ndvi_dyn_logit - dpred_ndvi$diff_logit) * 2 - 1
dfit$ndvi_dt <- dpred_ndvi$ndvi_dt

# plot(ndvi_dyn ~ ndvi_dt, dfit)
# abline(0, 1, col = "red")


# Summarize predictors and export -----------------------------------------

ndvi_summ <- aggregate(ndvi_dt ~ vegfac2, dfit, mean_hdi)
ndvi_summ <- cbind(vegetation = ndvi_summ$vegfac2,
                   as.data.frame(ndvi_summ$ndvi_dt))
colnames(ndvi_summ) <- c("vegetation", "mean", "hdi_lower_95", "hdi_upper_95")

elev_summ <- mean_hdi(dfit$elevation, name = "")
names(elev_summ) <- c("mean", "hdi_lower_95", "hdi_upper_95")

data_summ <- list(ndvi = ndvi_summ,
                  elevation = elev_summ)

saveRDS(data_summ, file.path("data", "flammability indices",
                             "ndvi_elevation_summary.rds"))

# Prior check for b_ndvi --------------------------------------------------

# model it at log scale, with normal random effect.

# check for mu
priorscale_b_ndvi_mu_log = 2
priormean_b_ndvi_mu_log = log(100)
b_ndvi_mu_log <- rnorm(1, priormean_b_ndvi_mu_log, priorscale_b_ndvi_mu_log)
b_ndvi_mu <- exp(b_ndvi_mu_log) * (-1)

curve(plogis(b_ndvi_mu * (x - 0.5) ^ 2), ylim = c(0, 1), n = 301,
      col = rgb(0, 0, 0, 0.08))

for(i in 1:300) {
  b_ndvi_mu_log <- rnorm(1, priormean_b_ndvi_mu_log, priorscale_b_ndvi_mu_log)
  b_ndvi_mu <- exp(b_ndvi_mu_log) * (-1)
  curve(plogis(b_ndvi_mu * (x - 0.5) ^ 2), add = T, n = 301,
        col = rgb(0, 0, 0, 0.08))
}


# Compile stan model ------------------------------------------------------

smodel <- stan_model("flammability indices/flammability_indices.stan")

# Model fit ----------------------------------------------------------------

sdata <- list(
  y = dfit$burned, N = nrow(dfit), V = V,
  veg = veg,
  veg_num = dfit$vegnum2,
  ndvi = dfit$ndvi_dt,
  elev = dfit$elevation,
  slope = dfit$slope,
  north = dfit$northing,

  # For NDVI hyper-parameters
  priorscale_a_mu = 10,
  priorscale_a_sigma = 3,

  # b at log scale
  priormean_b_mu = log(100),
  priorscale_b_mu = 2,
  priorscale_b_sigma = 1.75,

  # optim at logit scale
  priorscale_o_mu = 3,
  priorscale_o_sigma = 1.5,

  # prior scale for the slope of topographic variables
  priorscale_b_topo = 10
)

mod <- sampling(
  smodel, sdata, refresh = 100, seed = 32425,
  # cores = 1, chains = 1, iter = 5,
  cores = 6, chains = 6, iter = 3000, warmup = 1000, thin = 1,
  control = list(adapt_delta = 0.95)
)
# 1472.29 / 60 = 24 min
# 4 div transitions after warmup (ignore them)
saveRDS(mod, file.path("data", "flammability indices",
                       "flammability_indices_samples.rds"))
mod <- readRDS(file.path("data", "flammability indices",
                         "flammability_indices_samples.rds"))

smod <- summary(mod)[[1]]
min(smod[, "n_eff"]) # 2728.666
max(smod[, "Rhat"])  # 1.001616

# extract parameters
samples <- as_draws_matrix(mod)
parnames <- colnames(samples)
samples_mean <- colMeans(samples)
names(samples_mean) <- parnames
best_iter <- which.max(samples[, "lp__"])

alpha <- samples[, grep("a\\[", parnames)]
optim_ndvi <- samples[, grep("o\\[", parnames)]
b_ndvi <- samples[, grep("b\\[", parnames)]
b_topo_z <- samples[, c("b_elev", "b_slope", "b_north")]
b_topo <- samples[, c("b_elev_ori", "b_slope_ori", "b_north_ori")]
intercepts <- samples[, grep("intercepts\\[", parnames)]

# Explore posterior densities -------------------------------------------------

fplot <- function(x, fac = veg_classes) {
  xvec <- as.vector(x)
  df <- data.frame(x = xvec, veg = rep(fac, each = nrow(x)))
  return(df)
}

ggplot(fplot(optim_ndvi, fac = veg_classes), aes(x)) +
  geom_density() +
  facet_wrap(vars(veg)) +
  xlab("optim ndvi")

ggplot(fplot(b_ndvi), aes(x)) +
  geom_density(adjust = 1.5) +
  facet_wrap(vars(veg)) +
  xlab("b ndvi props")

ggplot(fplot(b_ndvi, fac = colnames(b_ndvi)), aes(x)) +
  geom_density() +
  facet_wrap(vars(veg)) +
  xlab("b ndvi")

ggplot(fplot(b_topo_z, fac = colnames(b_topo_z)), aes(x)) +
  geom_density() +
  facet_wrap(vars(veg), scales = "free_y", nrow = 2) +
  xlab("b topo")

ggplot(fplot(alpha), aes(x)) +
  geom_density() +
  facet_wrap(vars(veg)) +
  xlab("alpha")


# DHARMa residuals --------------------------------------------------------

S <- nrow(alpha)
veg_mat <- model.matrix(~ vegfac2 - 1, dfit)
topo_mat <- model.matrix(~ elev_z + slope_z + north_z - 1, dfit)

int_term <- veg_mat %*% t(alpha)
topo_term <- topo_mat %*% t(b_topo_z)
ndvi_term <- matrix(NA, nrow(dfit), S)

for(i in 1:S) {
  print(i)
  for(v in 1:V) {
    # i = 1; v = 1
    rows <- veg_num == v
    ndvi_term[rows, i] <-
      (dfit$ndvi_dt[rows] - optim_ndvi[i, v, drop = T]) ^ 2 * b_ndvi[i, v, drop = T]
  }
}

pfit <- plogis(int_term + topo_term + ndvi_term)

ysim <- sapply(1:S, function(i) {
  rbinom(nrow(dfit), size = 1, prob = pfit[, i])
})

res4 <- createDHARMa(simulatedResponse = ysim,
                     observedResponse = dfit$burned,
                     fittedPredictedResponse = rowMeans(pfit))
plot(res4, rank = F)
dfit$res4 <- res4$scaledResiduals

ggplot(dfit, aes(ndvi_dt, res4)) +
  geom_smooth(method = "gam", method.args = list(family = betar()),
              formula = y ~ s(x, bs = "cr", k = 10)) +
  facet_wrap(vars(vegfac2)) +
  ylim(0, 1)

ggplot(dfit, aes(elevation, res4)) +
  geom_smooth(method = "gam", method.args = list(family = betar()),
              formula = y ~ s(x, bs = "cr", k = 10)) +
  ylim(0, 1)

ggplot(dfit, aes(slope, res4)) +
  geom_smooth(method = "gam", method.args = list(family = betar()),
              formula = y ~ s(x, bs = "cr", k = 10)) +
  ylim(0, 1)

ggplot(dfit, aes(northing, res4)) +
  geom_smooth(method = "gam", method.args = list(family = betar()),
              formula = y ~ s(x, bs = "cr", k = 10)) +
  ylim(0, 1)

# biutiful

# burn prob by veg type, observed and fitted
dfit$pfit <- rowMeans(pfit)
aggveg <- aggregate(cbind(burned, pfit) ~ vegfac2, dfit, mean)
aggveg$vegfac2 <- factor(aggveg$vegfac2, levels = veg_levels)
ggplot(aggveg) +
  geom_point(aes(vegfac2, pfit), color = "red", size = 3) +
  geom_point(aes(vegfac2, burned), color = "black", size = 3, alpha = 0.6)  +
  ylim(0, 0.2) +
  ylab("burn probability") +
  xlab("veg type")
# perfect.

aggveg <- aggregate(cbind(burned, pfit) ~ vegfac1, dfit, mean)

ggplot(aggveg) +
  geom_point(aes(vegfac1, pfit), color = "red", size = 3) +
  geom_point(aes(vegfac1, burned), color = "black", size = 3, alpha = 0.6)  +
  ylim(0, 0.2)
# gooood


# Exports -----------------------------------------------------------------

# export the mean and sd of the term that will be multiplied by beta in the
# spread model, so the landscape layer is passed standardized to the simulator:


# VFI[i] = a[v][i] + b[v][i] * (ndvi[i] - o[v]) ^ 2

a_mean <- samples_mean[grep("a\\[", parnames)]
b_mean <- samples_mean[grep("b\\[", parnames)]
o_mean <- samples_mean[grep("o\\[", parnames)]
names(a_mean) <- names(b_mean) <- names(o_mean) <- veg_classes

vfi_ori <- numeric(nrow(dfit))
for(i in 1:nrow(dfit)) {
  vfi_ori[i] <-
    a_mean[veg_num[i]] +
    b_mean[veg_num[i]] *
    (dfit$ndvi_dt[i] - o_mean[veg_num[i]]) ^ 2
}

# compute whole mean and sd
vfi_mean <- mean(vfi_ori)
vfi_sd <- sd(vfi_ori)

# TFI[i] = b_elev_ori * elevation + b_north_ori * north
b_elev_ori_mean <- samples_mean[grep("b_elev_ori", parnames)]
b_north_ori_mean <- samples_mean[grep("b_north_ori", parnames)]

tfi_ori <- b_elev_ori_mean * dfit$elevation +
           b_north_ori_mean * dfit$northing

# compute whole mean and sd
tfi_mean <- mean(tfi_ori)
tfi_sd <- sd(tfi_ori)

# standardized indices in the study area, to get extreme values
vfi_z <- (vfi_ori - vfi_mean) / vfi_sd
tfi_z <- (tfi_ori - tfi_mean) / tfi_sd

theta_hat <- list(
  a = a_mean,
  b = b_mean,
  o = o_mean,

  b_elev_ori = b_elev_ori_mean,
  b_north_ori = b_north_ori_mean,

  # mean and sd of vfi and tfi in the landscape, to standardize
  vfi_mean = vfi_mean,
  vfi_sd = vfi_sd,
  tfi_mean = tfi_mean,
  tfi_sd = tfi_sd,

  # 95 % HDI of the standardized indices
  vfi_z_hdi = mean_hdi(vfi_z, name = "vfi_"),
  tfi_z_hdi = mean_hdi(tfi_z, name = "tfi_"),

  # predictors mean and sds (maybe not used)
  northing_mean = mean(dfit$northing),
  northing_sd = sd(dfit$northing),

  elevation_mean = mean(dfit$elevation),
  elevation_sd = sd(dfit$elevation),

  slope_term_mean = mean(sin(dfit$slope * pi / 180)),
  slope_term_sd = sd(sin(dfit$slope * pi / 180)),

  notes =
    "Posterior means of the estimates from a regional logistic regression,
similar as the model from Barberá et al. 2024 (simplified).
Mean NDVI from the previous summer was used, detrending and at 2022 scale.
Means and sds of predictors/terms are exported to standardize landscape arrays,
so spread parameters (betas) are in similar scales.
See <flammability_indices.R> for details."
)

saveRDS(theta_hat, file.path("data", "flammability indices",
                             "flammability_indices.rds"))
