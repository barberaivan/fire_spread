# The ideal vegetation terms to include should be
# eta[v] =
#   a[v] +
#   b[v] * (ndvi - optim[v]) ^ 2

# But that requires too many parameters to fit. A reduction could be like this:
# v in {1 = forest, 2 = shrubland, 3 = grassland},

# To reduce dimensionality even more, some parameters are estimated by a
# logistic regression for fire ocurrence in 24 years at the regional scale.

# b_ndvi could be defined as follows:
# b[2] = b[1] * pi_shrubland,
# b[3] = b[1] * pi_grassland,
# with pi in [0, 1].
# With b[1] estimated in the fire spread model, but pi_ by the logistic regression.
# In addition, optim[v] is also estimated with the regional logistic regression.

# In this code we estimate the regional model to obtain
# optim[v] and pi_shrubland and pi_grassland.

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
mdetrend <- readRDS(file.path("..", "data", "NDVI_regional_data",
                              "ndvi_detrender_model.rds"))

v0 <- vect(file.path("..", "data", "NDVI_regional_data",
                     "ndvi_regional_points.shp"))

dveg <- readxl::read_excel("/home/ivan/Insync/Mapa vegetación WWF - Lara et al. 1999/clases de vegetacion y equivalencias.xlsx",
                           sheet = "Sheet2")

names(dveg) <- c("vegnum1", "vegfac1", "vegnum2", "vegfac2")
dveg$vegfac3 <- dveg$vegfac2
dveg$vegfac3[grep("forest", dveg$vegfac2)] <- "Forest"
dveg$vegfac3 <- factor(dveg$vegfac3, levels = c("Forest", "Shrubland", "Grassland", "Non burnable"))
dveg$vegnum3 <- as.numeric(dveg$vegfac3)

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
                                    "vegfac3",
                                    "vegnum3",
                                    "ndvi_dyn",
                                    "year_dyn",
                                    "ndvi_22",
                                    "elevation",
                                    "slope",
                                    "northing")])
dfit <- data[rows_fit, ]; nrow(dfit); nrow(data)# OK

veg_levels <- unique(dveg$vegfac2[dveg$vegnum2 < 6] %>% as.character)

dfit$vegfac3 <- factor(dfit$vegfac3 %>% as.character,
                       levels = c("Forest", "Shrubland", "Grassland"))

dfit <- dfit[complete.cases(dfit), ]

veg <- model.matrix(~ vegfac3 - 1, dfit)
veg_num <- dfit$vegnum3
V <- ncol(veg)
veg_classes <- levels(dfit$vegfac3)

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

smodel <- stan_model("ndvi_effects.stan")

# Model fit ----------------------------------------------------------------

sdata <- list(
  y = dfit$burned, N = nrow(dfit), V = V,
  veg = veg,
  veg_num = dfit$vegnum3,
  ndvi = dfit$ndvi_dt,
  elev = dfit$elevation,
  slope = dfit$slope,
  north = dfit$northing,

  # prior parameters
  priorscale_alpha = 10,
  priorscale_b_ndvi_log = 2,
  priormean_b_ndvi_log = log(100),
  priorscale_b = 10
)

mod <- sampling(
  smodel, sdata, refresh = 100, seed = 32425,
  # cores = 1, chains = 1, iter = 5,
  cores = 6, chains = 6, iter = 3000, warmup = 1000, thin = 1,
  control = list(adapt_delta = 0.9)
)
saveRDS(mod, file.path("..", "data", "NDVI_regional_data",
                       "ndvi_effects_samples.rds"))
# mod <- readRDS(file.path("..", "data", "NDVI_regional_data",
#                          "ndvi_effects_samples.rds"))

smod <- summary(mod)[[1]]
min(smod[, "n_eff"]) # 2788.108
max(smod[, "Rhat"])  # 1.001979

pairs(mod, pars = c("optim_ndvi", "b_ndvi"))
pairs(mod, pars = c("optim_ndvi[3]", "pi_ndvi_grass", "optim_ndvi[2]", "pi_ndvi_shrub"))

# extract parameters
samples <- as_draws_matrix(mod)
parnames <- colnames(samples)
best_iter <- which.max(samples[, "lp__"])

optim_ndvi <- samples[, grep("optim_ndvi", parnames)]
b_ndvi_pi <- samples[, c("pi_ndvi_grass", "pi_ndvi_shrub")]
b_ndvi <- samples[, grep("b_ndvi\\[", parnames)]
b_topo_z <- samples[, c("b_elev", "b_slope", "b_north")]
b_topo <- samples[, c("b_elev_ori", "b_slope_ori", "b_north_ori")]
intercepts <- samples[, grep("intercepts\\[", parnames)]
alpha <- samples[, grep("alpha\\[", parnames)]

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

ggplot(fplot(b_ndvi_pi, fac = colnames(b_ndvi_pi)), aes(x)) +
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
veg_mat <- model.matrix(~ vegfac3 - 1, dfit)
topo_mat <- model.matrix(~ elev_z + slope_z + north_z - 1, dfit)

int_term <- veg_mat %*% t(alpha)
topo_term <- topo_mat %*% t(b_topo_z)
ndvi_term <- matrix(NA, nrow(dfit), S)
optim_reps <- outer(rep(1, nrow(dfit)), optim_ndvi)
squared_diff <- (dfit$ndvi_max - optim_reps) ^ 2
for(i in 1:S) {
  ndvi_term[, i] <- (squared_diff[, i, ] * veg_mat) %*% b_ndvi[i, , drop = T]
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
  facet_wrap(vars(vegfac3)) +
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

b_ndvi_for <- samples[, c("b_ndvi_for")]
pi_ndvi <- samples[, c("pi_ndvi_shrub", "pi_ndvi_grass")]

optim_ndvi_mean <- colMeans(optim_ndvi)
names(optim_ndvi_mean) <- veg_classes

pi_ndvi_mean <- colMeans(pi_ndvi)
names(pi_ndvi_mean) <- veg_classes[2:3]

oprep <- outer(rep(1, nrow(dfit)), optim_ndvi_mean)
ndvi_diff_quad <- rowSums((oprep - c(dfit$ndvi_dt)) ^ 2 * veg)
ndvi_dq_mean <- mean(ndvi_diff_quad)
ndvi_dq_sd <- sd(ndvi_diff_quad)

hist(b_ndvi_for * ndvi_dq_sd) # andaría por -1 si usara el término estandarizado

theta_hat <- list(
  optim_ndvi_mean = optim_ndvi_mean,
  pi_ndvi_mean = pi_ndvi_mean,

  ndvi_diffquad_mean = ndvi_dq_mean,
  ndvi_diffquad_sd = ndvi_dq_sd,

  elevation_mean = mean(dfit$elevation),
  elevation_sd = sd(dfit$elevation),

  notes =
"Posterior means of the estimates from a regional logistic regression,
similar as the model from Barberá et al. 2024 (simplified).
Mean NDVI from the previous summer was used, detrending and at 2022 scale.
dq stands for squared difference.
See <NDVI effects at the regional scale.R> for details."
)

saveRDS(theta_hat, file.path("..", "data", "NDVI_regional_data",
                             "ndvi_optim_and_proportion.rds"))
