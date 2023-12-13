# To reduce the dimensionality of predictors in the spread model we build
# vegetation and topographic flammability indices using a logistic regression
# that predicts burned vs unburned pixels as a function of these covariates.
# We use the same data as in the regional-scale model for the fire drivers
# interactions study, and fit a simpler model, not including FWI.

# Fitting a logistic model for burn probability, similar as in
# fire_drivers_interactions

# ndvi_max is ndvi_max_max, the maximum of all summer maximums.
# ndvi_mean is ndvi_mean_max, the maximum of all summer means.

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_minimal())
library(terra)
library(rstan)
library(mgcv)
library(DHARMa)

# Prepare data ------------------------------------------------------------

v <- vect(file.path("data", "points_burned-unburned_logistic_static_variables_2.shp"))

d <- as.data.frame(values(v))
d <- rename(d, ndvi_mean = ndvi_mean_, ndvi_max = ndvi_max_m)
names(d)
# remove unburnable
d <- d[d$vegetation > 1, ]
nrow(d) # 11943

# recode vegetation and remove unuseful categories
burned_veg <- read.csv("/home/ivan/Insync/Fire drivers interactions/fire_drivers_interactions/data/data_burned_and_available_area_by_vegetation_dryforest2.csv")
d$vegetation_code <- d$vegetation

# code dry forest A as subalpine forest
d$vegetation_code[d$vegetation_code == 4] <- 2
# plantation as dry forest B
d$vegetation_code[d$vegetation_code == 9] <- 5
# anthrop as shrubland
d$vegetation_code[d$vegetation_code == 8] <- 6

data <- left_join(d, burned_veg[, c("vegetation_code", "vegetation_class")],
               by = "vegetation_code")

data$vegetation_class[data$vegetation_class == "Dry forest B"] <- "Dry forest"

unique(data$vegetation_class)

veg_levels <- c(
  "Wet forest",
  "Subalpine forest",
  "Dry forest",
  "Shrubland",
  "Steppe and grassland"
)

veg_labels <- c(
  "Wet forest",
  "Subalpine\nforest",
  "Dry forest",
  "Shrubland",
  "Steppe and\ngrassland"
)

data$vegetation_class <- factor(data$vegetation_class, levels = veg_levels,
                             labels = veg_labels)

# compute slope-weighted northing
# (see code <northing importance function.R>)
data$northing <- cos(data$aspect * pi / 180) * plogis(-5 + 0.35 * data$slope)


# NDVI max and mean by veg type -------------------------------------------

# par(mfrow = c(3, 2))
# for(v in 1:length(veg_levels)) {
#   dd <- data[data$vegetation_class == veg_labels[v], ]
#   
#   dmax <- density(dd$ndvi_max, from = 0, to = 1)
#   dmean <- density(dd$ndvi_mean, from = 0, to = 1)
#   
#   plot(dmax, col = 1, ylim = c(0, max(dmax$y, dmean$y)),
#        main = veg_levels[v])
#   lines(dmean, col = 2)
# }
# par(mfrow = c(1, 1))

# Model scales and prior checks -------------------------------------------

# Our aim is to get parameters applicable to the predictors in the original
# scale, so that one does not have to transform them. But the real parameters
# (raw) will be sampled from uniform or std_normal distributions to improve
# the sampling efficiency.

# Prior checks

# elevation # [201, 2253]
elev_optim <- 2000
b <- 25 / sd((data$elevation - mean(data$elevation)) ^ 2)
curve(plogis(-b * (x - elev_optim) ^ 2),
      from = min(data$elevation), to = max(data$elevation),
      ylim = c(0, 1), n = 300)
# b_elev_raw ~ std_normal(); upper = 0
# b_elev <- b_elev_raw * 25 / sd((elevation - mean(elevation)) ^ 2)
# elev_optim_raw ~ unif(0, 1)
# elev_optim <- elev_optim_raw * 2500

# tpi # [0, 1]
tpi_optim <- 0.5
b <- 25 / sd((data$TPI1k - mean(data$TPI1k)) ^ 2)
curve(plogis(-b * (x - tpi_optim) ^ 2),
      from = min(data$TPI1k), to = max(data$TPI1k),
      ylim = c(0, 1), n = 300)
# b_tpi_raw ~ std_normal(); upper = 0
# b_tpi <- b_tpi_raw * 25 / sd((tpi - mean(TPI1k)) ^ 2)
# tpi_optim ~ unif(0, 1)

# ndvi # [0, 1]
ndvi_optim <- 0.7
b <- 0.5 / sd((data$ndvi_max - mean(data$ndvi_max)) ^ 2)
curve(plogis(-b * (x - ndvi_optim) ^ 2),
      from = min(data$ndvi_max), to = max(data$ndvi_max),
      ylim = c(0, 1), n = 300)
# b_ndvi_raw ~ std_normal(); upper = 0
# b_ndvi <- b_ndvi_raw * 25 / sd((ndvi - mean(ndvi)) ^ 2)
# ndvi_optim ~ unif(0, 1)

# pp # [381, 2175]
pp_optim <- 1200
b <- 0.5 / sd((data$pp - mean(data$pp)) ^ 2)
curve(plogis(-b * (x - pp_optim) ^ 2),
      from = min(data$pp), to = max(data$pp),
      ylim = c(0, 1), n = 300)
# b_pp_raw ~ std_normal(); upper = 0
# b_pp <- b_pp_raw * 25 / sd((pp - mean(pp)) ^ 2)
# ndvi_optim ~ unif(300, 2200)

# northing # [-1, 1]
b <- 5 / sd(data$northing) # sd = 3 is regularizer, >0
curve(plogis(b * (x - mean(data$northing))),
      from = min(data$northing), to = max(data$northing),
      ylim = c(0, 1))
# b_north_raw ~ std_normal(), lower = 0
# b_north = b_north_raw * 5 / sd(data$northing)

# slope # [0, 73.55016]
b <- 5 / sd(data$slope)
curve(plogis(b * (x - mean(data$slope))),
      from = min(data$slope), to = max(data$slope),
      ylim = c(0, 1))
# b_north_raw ~ std_normal(), lower = 0
# b_north = b_north_raw * 5 / sd(data$slope)


# Data to fit models -------------------------------------------------------

rows_fit <- complete.cases(data[, c("burned", "vegetation_class", 
                                    "ndvi_max", "ndvi_mean",
                                    "elevation", "TPI1k", "northing", "slope",
                                    "pp", "dist_human", "dist_roads")])
dfit <- data[rows_fit, ]; nrow(dfit); nrow(data)# OK

veg <- model.matrix(~ vegetation_class - 1, dfit)
V <- ncol(veg)

# smodel <- stan_model("landscape_flammability.stan")

# Model fit (ndvi_max) ----------------------------------------------------

ndvi_mean <- aggregate(ndvi_max ~ vegetation_class, dfit, mean)[, 2]
# ndvi_scale is sd((ndvi - ndvi_mean) ^ 2) by veg_type
ndvi_scale <- numeric(V)
for(v in 1:V) {
  nn <- dfit$ndvi_max[veg[, v] == 1]
  ndvi_scale[v] <- sd((nn - ndvi_mean[v]) ^ 2)
}

sdata <- list(
  y = dfit$burned, N = nrow(dfit), V = V,
  
  veg = veg, ndvi = dfit$ndvi_max,
  elev = dfit$elevation, tpi = dfit$TPI1k, 
  north = dfit$northing, slope = dfit$slope,
  
  ndvi_mean = ndvi_mean, ndvi_scale = ndvi_scale,
  
  prior_sd_int =  5, 
  prior_sd_linear = 5, 
  prior_sd_quad = 15
)

# tried to optimize, but did not work well

# m_max <- sampling(
#   smodel, sdata, refresh = 10, seed = 32425,
#   cores = 4, 
#   chains = 4, 
#   iter = 2000, warmup = 500,
# )
# saveRDS(m_max, "landscape_flammability_model_samples_ndvi-max.rds")
# takes 40 min
m_max <- readRDS("landscape_flammability_model_samples_ndvi-max.rds")
sm1 <- summary(m_max)[[1]]

# elevation-tpi competition?
pairs(m_max, pars = c("optim_elev", "optim_tpi", "b_elev", "b_tpi"))

# extract parameters
best_iter <- which.max(as.matrix(m_max, "lp__") %>% as.numeric)

b_elev <- as.matrix(m_max, "b_elev") %>% as.numeric
optim_elev <- as.matrix(m_max, "optim_elev") %>% as.numeric

b_tpi <- as.matrix(m_max, "b_tpi") %>% as.numeric
optim_tpi <- as.matrix(m_max, "optim_tpi") %>% as.numeric

b_ndvi <- as.matrix(m_max, "b_ndvi")
optim_ndvi <- as.matrix(m_max, "optim_ndvi")

b_north <- as.matrix(m_max, "b_north") %>% as.numeric
b_slope <- as.matrix(m_max, "b_slope") %>% as.numeric

intercepts <- as.matrix(m_max, "intercepts")

S <- length(b_elev)

# Partial predictions (ndvi_max) -------------------------------------------

par(mfrow = c(2, 2))

# elevation
plot(1000, 0.5, col = "white", ylim = c(0, 1), xlim = c(200, 2200),
     xlab = "Elevation", ylab = "Burn probability")
abline(v = quantile(dfit$elevation, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_elev[s] * (x - optim_elev[s]) ^ 2), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_elev[best_iter] * (x - optim_elev[best_iter]) ^ 2), add = TRUE, 
      col = rgb(0, 1, 1, 1), lwd = 2)

# tpi # positive effect, maybe because most vegetation reaches up to 0.7?
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$TPI1k),
     xlab = "Topographic position", ylab = "Burn probability")
abline(v = quantile(dfit$TPI1k, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_tpi[s] * (x - optim_tpi[s]) ^ 2), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_tpi[best_iter] * (x - optim_tpi[best_iter]) ^ 2), add = TRUE, 
      col = rgb(0, 1, 1, 1), lwd = 2)


# north (le había puesto upper = 0, pero debería ser lower = 0)
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$northing),
     ylab = "Burn probability", xlab = "Northing")
abline(v = quantile(dfit$northing, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_north[s] * x), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_north[best_iter] * x), add = TRUE,
      col = rgb(0, 1, 1, 1), lwd = 2)

# slope (le había puesto upper = 0, pero debería ser lower = 0)
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$slope),
     ylab = "Burn probability", xlab = "Slope")
abline(v = quantile(dfit$slope, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_slope[s] * x), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_slope[best_iter] * x), add = TRUE,
      col = rgb(0, 1, 1, 1), lwd = 2)

par(mfrow = c(1, 1))


# ndvi
samples <- sample(1:S, 500)
par(mfrow = c(2, 3))
for(v in 1:V) {
  plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$ndvi_max),
       main = veg_levels[v],
       ylab = "Burn probability", xlab = "NDVI")
  abline(v = quantile(dfit$ndvi_max[dfit$vegetation_class == veg_labels[v]],
                      prob = c(0.05, 0.95)),
         lty = 2)
  for(s in samples) {
    curve(plogis(b_ndvi[s, v] * (x - optim_ndvi[s, v]) ^ 2), add = TRUE, 
          col = rgb(0, 0, 0, 0.05))
  }
  curve(plogis(b_ndvi[best_iter, v] * (x - optim_ndvi[best_iter, v]) ^ 2), add = TRUE, 
        col = rgb(0, 1, 1, 1), lwd = 2)
  
}
par(mfrow = c(1, 1))

# intercepts
int_df <- data.frame(intercept = as.numeric(intercepts),
                     veg = rep(veg_levels, each = nrow(intercepts)))
boxplot(intercept ~ veg, int_df)

# Ahora los efectos son mucho mayores!


# Model fit (ndvi_mean) ----------------------------------------------------

ndvi_mean <- aggregate(ndvi_mean ~ vegetation_class, dfit, mean)[, 2]
# ndvi_scale is sd((ndvi - ndvi_mean) ^ 2) by veg_type
ndvi_scale <- numeric(V)
for(v in 1:V) {
  nn <- dfit$ndvi_mean[veg[, v] == 1]
  ndvi_scale[v] <- sd((nn - ndvi_mean[v]) ^ 2)
}

sdata <- list(
  y = dfit$burned, N = nrow(dfit), V = V,
  
  veg = veg, ndvi = dfit$ndvi_mean,
  elev = dfit$elevation, tpi = dfit$TPI1k, 
  north = dfit$northing, slope = dfit$slope,
  
  ndvi_mean = ndvi_mean, ndvi_scale = ndvi_scale,
  
  prior_sd_int =  5, 
  prior_sd_linear = 5, 
  prior_sd_quad = 15
)

# m_mean <- sampling(
#   smodel, 
#   sdata, refresh = 20, seed = 32425,
#   cores = 4,
#   chains = 4,
#   iter = 2000, warmup = 500,
# )
# saveRDS(m_mean, "landscape_flammability_model_samples_ndvi-mean.rds")
# 2360.59 / 60 = 39 min 
m_mean <- readRDS("landscape_flammability_model_samples_ndvi-mean.rds")

sm1 <- summary(m_mean)[[1]]

# elevation-tpi competition?
pairs(m_mean, pars = c("optim_elev", "optim_tpi", "b_elev", "b_tpi"))

# extract parameters
best_iter <- which.max(as.matrix(m_mean, "lp__") %>% as.numeric)

b_elev <- as.matrix(m_mean, "b_elev") %>% as.numeric
optim_elev <- as.matrix(m_mean, "optim_elev") %>% as.numeric

b_tpi <- as.matrix(m_mean, "b_tpi") %>% as.numeric
optim_tpi <- as.matrix(m_mean, "optim_tpi") %>% as.numeric

b_ndvi <- as.matrix(m_mean, "b_ndvi")
optim_ndvi <- as.matrix(m_mean, "optim_ndvi")

b_north <- as.matrix(m_mean, "b_north") %>% as.numeric
b_slope <- as.matrix(m_mean, "b_slope") %>% as.numeric

intercepts <- as.matrix(m_mean, "intercepts")

S <- length(b_elev)

# Partial predictions (ndvi_mean) -------------------------------------------

par(mfrow = c(2, 2))

# elevation
plot(1000, 0.5, col = "white", ylim = c(0, 1), xlim = c(200, 2200),
     xlab = "Elevation", ylab = "Burn probability")
abline(v = quantile(dfit$elevation, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_elev[s] * (x - optim_elev[s]) ^ 2), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_elev[best_iter] * (x - optim_elev[best_iter]) ^ 2), add = TRUE, 
      col = rgb(0, 1, 1, 1), lwd = 2)

# tpi # positive effect, maybe because most vegetation reaches up to 0.7?
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$TPI1k),
     xlab = "Topographic position", ylab = "Burn probability")
abline(v = quantile(dfit$TPI1k, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_tpi[s] * (x - optim_tpi[s]) ^ 2), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_tpi[best_iter] * (x - optim_tpi[best_iter]) ^ 2), add = TRUE, 
      col = rgb(0, 1, 1, 1), lwd = 2)


# north (le había puesto upper = 0, pero debería ser lower = 0)
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$northing),
     ylab = "Burn probability", xlab = "Northing")
abline(v = quantile(dfit$northing, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_north[s] * x), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_north[best_iter] * x), add = TRUE,
      col = rgb(0, 1, 1, 1), lwd = 2)

# slope (le había puesto upper = 0, pero debería ser lower = 0)
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$slope),
     ylab = "Burn probability", xlab = "Slope")
abline(v = quantile(dfit$slope, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_slope[s] * x), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_slope[best_iter] * x), add = TRUE,
      col = rgb(0, 1, 1, 1), lwd = 2)

par(mfrow = c(1, 1))


# ndvi
samples <- sample(1:S, 4000)
par(mfrow = c(2, 3))
for(v in 1:V) {
  plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$ndvi_mean),
       main = veg_levels[v],
       ylab = "Burn probability", xlab = "NDVI")
  abline(v = quantile(dfit$ndvi_mean[dfit$vegetation_class == veg_labels[v]],
                      prob = c(0.05, 0.95)),
         lty = 2)
  for(s in samples) {
    curve(plogis(b_ndvi[s, v] * (x - optim_ndvi[s, v]) ^ 2), add = TRUE, 
          col = rgb(0, 0, 0, 0.01))
  }
  curve(plogis(b_ndvi[best_iter, v] * (x - optim_ndvi[best_iter, v]) ^ 2), add = TRUE, 
        col = rgb(0, 1, 1, 1), lwd = 2)
  
}
par(mfrow = c(1, 1))

# intercepts
int_df <- data.frame(intercept = as.numeric(intercepts),
                     veg = rep(veg_levels, each = nrow(intercepts)))
boxplot(intercept ~ veg, int_df)



# DHARMa (ndvi_mean) ------------------------------------------------------

b_ndvi[best_iter, ]; optim_ndvi[best_iter, ]

b_elev[best_iter]; optim_elev[best_iter]
b_tpi[best_iter]; optim_tpi[best_iter]

b_north[best_iter]
b_slope[best_iter]

intercepts[best_iter, ]

# topographic logit
topo_logit <- b_elev[best_iter] * (dfit$elevation - optim_elev[best_iter]) ^ 2 +
              b_tpi[best_iter] *  (dfit$TPI1k - optim_tpi[best_iter]) ^ 2 +
              b_north[best_iter] * dfit$northing +
              b_slope[best_iter] * dfit$slope

# veg logits
veg_alpha <- sdata$veg %*% intercepts[best_iter, ]
ndvi_rep <- matrix(rep(dfit$ndvi_mean, V), nrow(dfit), V)
optims_rep <- outer(rep(1, nrow(dfit)), optim_ndvi[best_iter, ])
ndvi_quad <- (ndvi_rep - optims_rep) ^ 2 * sdata$veg
ndvi_logit <- ndvi_quad %*% b_ndvi[best_iter, ]

# fitted p
fit <- plogis(veg_alpha + ndvi_logit + topo_logit) %>% as.numeric()
hist(fit)
mean(fit)
mean(dfit$burned)

# Con la parametrización del simpler anda bien, pero con el otro los intercepts 
# estaban mal

dfit$pfit <- fit

ysim <- matrix(rbinom(nrow(dfit) * S, size = 1, prob = fit),
               nrow(dfit), S)

res <- createDHARMa(ysim, dfit$burned, fittedPredictedResponse = fit)
plot(res)
plotResiduals(res, form = dfit$vegetation_class)

# da hermosoooo
dfit$res <- res$scaledResiduals

ggplot(dfit, aes(ndvi_mean, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(elevation, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(TPI1k, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(slope, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(northing, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

# The ndvi_mean model fits perfectly

# R2
var(fit) / (var(fit) + mean(fit * (1 - fit))) # 0.17

# compare to a very good R2
bgud <- 3.3
curve(plogis(bgud * x), ylim = c(0, 1), from = -0.5, to = 0.5)
xs <- seq(-0.5, 0.5, length.out = 10000)
ff <- plogis(bgud * xs)
var(ff) / (var(ff) + mean(ff * (1 - ff))) # 0.17 con bgud = 3.3
# O sea, es análogo a esto.

# intercepts by veg type --------------------------------------------------

boxplot(fit ~ dfit$vegetation_class)
intercepts[best_iter, ] %>% plogis
# intercepts and mean burn prob are very different, because 
# important predictors vary with vegetation class.

# NDVI model --------------------------------------------------------------

# is ndvi explained by topography and climate?

ndvi_model <- gam(
  ndvi_mean ~ 
    vegetation_class +
    s(TPI1k, k = 5, bs = "cr") + s(elevation, k = 5, bs = "cr") +
    s(slope, k = 5, bs = "cr") + s(northing, k = 5, bs = "cr") +
    s(pp, k = 5, bs = "cr") + s(temp, k = 5, bs = "cr"),
  family = betar(), method = "REML",
  data = dfit
)

plot(ndvi_model)
summary(ndvi_model)
plot(dfit$ndvi_mean ~ fitted(ndvi_model))
# R2 ~ 0.5




# Model fit (ndvi_mean, pp dist) -----------------------------------------

smodel_ppdist <- stan_model("landscape flammability indices/landscape_flammability_pp_dist.stan")

ndvi_mean <- aggregate(ndvi_mean ~ vegetation_class, dfit, mean)[, 2]
# ndvi_scale is sd((ndvi - ndvi_mean) ^ 2) by veg_type
ndvi_scale <- numeric(V)
for(v in 1:V) {
  nn <- dfit$ndvi_mean[veg[, v] == 1]
  ndvi_scale[v] <- sd((nn - ndvi_mean[v]) ^ 2)
}

sdata <- list(
  y = dfit$burned, N = nrow(dfit), V = V,
  
  veg = veg, ndvi = dfit$ndvi_mean,
  pp = dfit$pp,
  
  elev = dfit$elevation, tpi = dfit$TPI1k, 
  north = dfit$northing, slope = dfit$slope,
  
  dist_human = as.numeric(scale(dfit$dist_human)),
  dist_roads = as.numeric(scale(dfit$dist_roads)),
  
  ndvi_mean = ndvi_mean, ndvi_scale = ndvi_scale,
  
  prior_sd_int =  5, 
  prior_sd_linear = 5, 
  prior_sd_quad = 15
)

# m_mean <- sampling(
#   smodel_ppdist,
#   sdata, refresh = 20, seed = 32425,
#   cores = 4,
#   chains = 4,
#   iter = 2000, warmup = 500,
#   # chains = 1, iter = 10, cores = 1
# )
# saveRDS(
#   m_mean, 
#   file.path("landscape flammability indices",
#             "landscape_flammability_model_samples_ndvi-mean_pp-dist.rds")
# )
# 3344.42 / 60 = 56 min 
m_mean <- readRDS(
  file.path("landscape flammability indices",
            "landscape_flammability_model_samples_ndvi-mean_pp-dist.rds")
)

sm1 <- summary(m_mean)[[1]]
# Low n_eff for optim_ndvi[4]

# correlation between parameters
pairs(m_mean, pars = c("optim_elev", "optim_tpi", "b_elev", "b_tpi"))
pairs(m_mean, pars = c("optim_ndvi", "b_pp", "optim_pp"))
pairs(m_mean, pars = c("b_ndvi", "b_pp", "optim_pp"))
pairs(m_mean, pars = c("b_pp", "optim_pp", "b_human", "b_roads"))
pairs(m_mean, pars = c("b_elev", "optim_elev", "b_human", "b_roads"))
pairs(m_mean, pars = c("b_pp", "optim_pp", "b_elev", "optim_elev"))
# no correlation between 

# extract parameters
best_iter <- which.max(as.matrix(m_mean, "lp__") %>% as.numeric)

b_elev <- as.matrix(m_mean, "b_elev") %>% as.numeric
optim_elev <- as.matrix(m_mean, "optim_elev") %>% as.numeric

b_tpi <- as.matrix(m_mean, "b_tpi") %>% as.numeric
optim_tpi <- as.matrix(m_mean, "optim_tpi") %>% as.numeric

b_ndvi <- as.matrix(m_mean, "b_ndvi")
optim_ndvi <- as.matrix(m_mean, "optim_ndvi")

b_pp <- as.matrix(m_mean, "b_pp")
optim_pp <- as.matrix(m_mean, "optim_pp")

b_north <- as.matrix(m_mean, "b_north") %>% as.numeric
b_slope <- as.matrix(m_mean, "b_slope") %>% as.numeric

b_human <- as.matrix(m_mean, "b_human") %>% as.numeric
b_roads <- as.matrix(m_mean, "b_roads") %>% as.numeric

intercepts <- as.matrix(m_mean, "intercepts")

S <- length(b_elev)

# Partial predictions (ndvi_mean, pp dist) --------------------------------

par(mfrow = c(2, 2))
# elevation
plot(1000, 0.5, col = "white", ylim = c(0, 1), xlim = c(200, 2200),
     xlab = "Elevation", ylab = "Burn probability")
abline(v = quantile(dfit$elevation, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_elev[s] * (x - optim_elev[s]) ^ 2), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_elev[best_iter] * (x - optim_elev[best_iter]) ^ 2), add = TRUE, 
      col = rgb(0, 1, 1, 1), lwd = 2)

# tpi # positive effect, maybe because most vegetation reaches up to 0.7?
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$TPI1k),
     xlab = "Topographic position", ylab = "Burn probability")
abline(v = quantile(dfit$TPI1k, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_tpi[s] * (x - optim_tpi[s]) ^ 2), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_tpi[best_iter] * (x - optim_tpi[best_iter]) ^ 2), add = TRUE, 
      col = rgb(0, 1, 1, 1), lwd = 2)


# north (le había puesto upper = 0, pero debería ser lower = 0)
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$northing),
     ylab = "Burn probability", xlab = "Northing")
abline(v = quantile(dfit$northing, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_north[s] * x), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_north[best_iter] * x), add = TRUE,
      col = rgb(0, 1, 1, 1), lwd = 2)

# slope (le había puesto upper = 0, pero debería ser lower = 0)
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$slope),
     ylab = "Burn probability", xlab = "Slope")
abline(v = quantile(dfit$slope, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_slope[s] * x), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_slope[best_iter] * x), add = TRUE,
      col = rgb(0, 1, 1, 1), lwd = 2)

par(mfrow = c(1, 1))


# ndvi
samples <- sample(1:S, 1000)
par(mfrow = c(2, 3))
# ndvi
for(v in 1:V) {
  plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$ndvi_mean),
       main = veg_levels[v],
       ylab = "Burn probability", xlab = "NDVI")
  abline(v = quantile(dfit$ndvi_mean[dfit$vegetation_class == veg_labels[v]],
                      prob = c(0.05, 0.95)),
         lty = 2)
  for(s in samples) {
    curve(plogis(b_ndvi[s, v] * (x - optim_ndvi[s, v]) ^ 2), add = TRUE, 
          col = rgb(0, 0, 0, 0.01))
  }
  curve(plogis(b_ndvi[best_iter, v] * (x - optim_ndvi[best_iter, v]) ^ 2), add = TRUE, 
        col = rgb(0, 1, 1, 1), lwd = 2)
  
}

# pp
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$pp),
     xlab = "Annual precipitation (mm)", ylab = "Burn probability")
abline(v = quantile(dfit$pp, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_pp[s] * (x - optim_pp[s]) ^ 2), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_pp[best_iter] * (x - optim_pp[best_iter]) ^ 2), add = TRUE, 
      col = rgb(0, 1, 1, 1), lwd = 2)
par(mfrow = c(1, 1))


## human occupation
plot(density(b_human)) # burns more far from human settlements?
plot(density(b_roads))

hu_seq <- seq(min(dfit$dist_human), max(dfit$dist_human), length.out = 150)
hu_seq_z <- (hu_seq - mean(dfit$dist_human)) / sd(dfit$dist_human)

ro_seq <- seq(min(dfit$dist_roads), max(dfit$dist_roads), length.out = 150)
ro_seq_z <- (ro_seq - mean(dfit$dist_roads)) / sd(dfit$dist_roads)

par(mfrow = c(1, 2))

plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$dist_human),
     ylab = "Burn probability", xlab = "Distance from human settlements (m)")
abline(v = quantile(dfit$dist_human, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  lines(plogis(b_human[s] * hu_seq_z) ~ hu_seq, 
        col = rgb(0, 0, 0, 0.05))
}
lines(plogis(b_human[best_iter] * hu_seq_z) ~ hu_seq, 
      col = rgb(0, 1, 1, 1), lwd = 2)

plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$dist_roads),
     ylab = "Burn probability", xlab = "Distance from roads (m)")
abline(v = quantile(dfit$dist_roads, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  lines(plogis(b_roads[s] * ro_seq_z) ~ ro_seq, 
        col = rgb(0, 0, 0, 0.05))
}
lines(plogis(b_roads[best_iter] * ro_seq_z) ~ ro_seq, 
      col = rgb(0, 1, 1, 1), lwd = 2)

par(mfrow = c(1, 1))


# intercepts
int_df <- data.frame(intercept = as.numeric(intercepts),
                     veg = rep(veg_levels, each = nrow(intercepts)))
boxplot(intercept ~ veg, int_df)



# DHARMa (ndvi_mean, pp dist) ---------------------------------------------

b_ndvi[best_iter, ]; optim_ndvi[best_iter, ]

b_elev[best_iter]; optim_elev[best_iter]
b_tpi[best_iter]; optim_tpi[best_iter]

b_north[best_iter]
b_slope[best_iter]

intercepts[best_iter, ]

# topographic logit
topo_logit <- b_elev[best_iter] * (dfit$elevation - optim_elev[best_iter]) ^ 2 +
  b_tpi[best_iter] *  (dfit$TPI1k - optim_tpi[best_iter]) ^ 2 +
  b_north[best_iter] * dfit$northing +
  b_slope[best_iter] * dfit$slope

# veg logits
veg_alpha <- sdata$veg %*% intercepts[best_iter, ]
ndvi_rep <- matrix(rep(dfit$ndvi_mean, V), nrow(dfit), V)
optims_rep <- outer(rep(1, nrow(dfit)), optim_ndvi[best_iter, ])
ndvi_quad <- (ndvi_rep - optims_rep) ^ 2 * sdata$veg
ndvi_logit <- ndvi_quad %*% b_ndvi[best_iter, ]

# pp
pp_logit <- b_pp[best_iter] *  (dfit$pp - optim_pp[best_iter]) ^ 2

# distance
dist_logit <- b_human[best_iter] * sdata$dist_human +
              b_roads[best_iter] * sdata$dist_roads

# fitted p
fit <- plogis(veg_alpha + ndvi_logit + topo_logit + pp_logit + dist_logit) %>% as.numeric()
hist(fit)
mean(fit)
mean(dfit$burned)

# Con la parametrización del simpler anda bien, pero con el otro los intercepts 
# estaban mal

dfit$pfit <- fit

ysim <- matrix(rbinom(nrow(dfit) * S, size = 1, prob = fit),
               nrow(dfit), S)

res <- createDHARMa(ysim, dfit$burned, fittedPredictedResponse = fit)
plot(res)

# da hermosoooo
dfit$res <- res$scaledResiduals

ggplot(dfit, aes(ndvi_mean, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(pp, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(elevation, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(TPI1k, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(slope, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(northing, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(dist_human, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(dist_roads, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")


# R2
var(fit) / (var(fit) + mean(fit * (1 - fit))) # 0.1869553, 
# before adding pp and distances it was 0.17
                                            

# Model fit (ndvi_mean, pp) -----------------------------------------

smodel_pp <- stan_model("landscape flammability indices/landscape_flammability_pp.stan")

ndvi_mean <- aggregate(ndvi_mean ~ vegetation_class, dfit, mean)[, 2]
# ndvi_scale is sd((ndvi - ndvi_mean) ^ 2) by veg_type
ndvi_scale <- numeric(V)
for(v in 1:V) {
  nn <- dfit$ndvi_mean[veg[, v] == 1]
  ndvi_scale[v] <- sd((nn - ndvi_mean[v]) ^ 2)
}

sdata <- list(
  y = dfit$burned, N = nrow(dfit), V = V,
  
  veg = veg, ndvi = dfit$ndvi_mean,
  pp = dfit$pp,
  
  elev = dfit$elevation, tpi = dfit$TPI1k, 
  north = dfit$northing, slope = dfit$slope,

  ndvi_mean = ndvi_mean, ndvi_scale = ndvi_scale,
  
  prior_sd_int =  5, 
  prior_sd_linear = 5, 
  prior_sd_quad = 15
)

# m_mean <- sampling(
#   smodel_pp,
#   sdata, refresh = 50, seed = 32425,
#   cores = 10, # to increase n_eff in shrubland ndvi parameters
#   chains = 10,
#   iter = 2000, warmup = 500
# )
# saveRDS(
#   m_mean,
#   file.path("landscape flammability indices",
#             "landscape_flammability_model_samples_ndvi-mean_pp.rds")
# )
# 6498.92 / 60 = 108 min 
m_mean <- readRDS(
  file.path("landscape flammability indices",
            "landscape_flammability_model_samples_ndvi-mean_pp.rds")
)

sm1 <- summary(m_mean)[[1]]
# Low n_eff for optim_ndvi[4]

# correlation between parameters
pairs(m_mean, pars = c("optim_elev", "optim_tpi", "b_elev", "b_tpi"))
pairs(m_mean, pars = c("optim_ndvi", "b_pp", "optim_pp"))
pairs(m_mean, pars = c("b_ndvi", "b_pp", "optim_pp"))
pairs(m_mean, pars = c("b_pp", "optim_pp", "b_elev", "optim_elev"))
# no correlation between ndvi and pp

# extract parameters
best_iter <- which.max(as.matrix(m_mean, "lp__") %>% as.numeric)

b_elev <- as.matrix(m_mean, "b_elev") %>% as.numeric
optim_elev <- as.matrix(m_mean, "optim_elev") %>% as.numeric

b_tpi <- as.matrix(m_mean, "b_tpi") %>% as.numeric
optim_tpi <- as.matrix(m_mean, "optim_tpi") %>% as.numeric

b_ndvi <- as.matrix(m_mean, "b_ndvi")
optim_ndvi <- as.matrix(m_mean, "optim_ndvi")

b_pp <- as.matrix(m_mean, "b_pp")
optim_pp <- as.matrix(m_mean, "optim_pp")

b_north <- as.matrix(m_mean, "b_north") %>% as.numeric
b_slope <- as.matrix(m_mean, "b_slope") %>% as.numeric

intercepts <- as.matrix(m_mean, "intercepts")

S <- length(b_elev)

# Partial predictions (ndvi_mean, pp) --------------------------------

par(mfrow = c(2, 2))
# elevation
plot(1000, 0.5, col = "white", ylim = c(0, 1), xlim = c(200, 2200),
     xlab = "Elevation", ylab = "Burn probability")
abline(v = quantile(dfit$elevation, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_elev[s] * (x - optim_elev[s]) ^ 2), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_elev[best_iter] * (x - optim_elev[best_iter]) ^ 2), add = TRUE, 
      col = rgb(0, 1, 1, 1), lwd = 2)

# tpi # positive effect, maybe because most vegetation reaches up to 0.7?
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$TPI1k),
     xlab = "Topographic position", ylab = "Burn probability")
abline(v = quantile(dfit$TPI1k, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_tpi[s] * (x - optim_tpi[s]) ^ 2), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_tpi[best_iter] * (x - optim_tpi[best_iter]) ^ 2), add = TRUE, 
      col = rgb(0, 1, 1, 1), lwd = 2)


# north (le había puesto upper = 0, pero debería ser lower = 0)
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$northing),
     ylab = "Burn probability", xlab = "Northing")
abline(v = quantile(dfit$northing, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_north[s] * x), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_north[best_iter] * x), add = TRUE,
      col = rgb(0, 1, 1, 1), lwd = 2)

# slope (le había puesto upper = 0, pero debería ser lower = 0)
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$slope),
     ylab = "Burn probability", xlab = "Slope")
abline(v = quantile(dfit$slope, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_slope[s] * x), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_slope[best_iter] * x), add = TRUE,
      col = rgb(0, 1, 1, 1), lwd = 2)

par(mfrow = c(1, 1))


# ndvi
samples <- sample(1:S, 1000)
par(mfrow = c(2, 3))
# ndvi
for(v in 1:V) {
  plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$ndvi_mean),
       main = veg_levels[v],
       ylab = "Burn probability", xlab = "NDVI")
  abline(v = quantile(dfit$ndvi_mean[dfit$vegetation_class == veg_labels[v]],
                      prob = c(0.05, 0.95)),
         lty = 2)
  for(s in samples) {
    curve(plogis(b_ndvi[s, v] * (x - optim_ndvi[s, v]) ^ 2), add = TRUE, 
          col = rgb(0, 0, 0, 0.01))
  }
  curve(plogis(b_ndvi[best_iter, v] * (x - optim_ndvi[best_iter, v]) ^ 2), add = TRUE, 
        col = rgb(0, 1, 1, 1), lwd = 2)
  
}

# pp
plot(1, type = "n", ylim = c(0, 1), xlim = range(dfit$pp),
     xlab = "Annual precipitation (mm)", ylab = "Burn probability")
abline(v = quantile(dfit$pp, prob = c(0.05, 0.95)), lty = 2)
for(s in sample(1:S, 500)) {
  curve(plogis(b_pp[s] * (x - optim_pp[s]) ^ 2), add = TRUE, 
        col = rgb(0, 0, 0, 0.05))
}
curve(plogis(b_pp[best_iter] * (x - optim_pp[best_iter]) ^ 2), add = TRUE, 
      col = rgb(0, 1, 1, 1), lwd = 2)
par(mfrow = c(1, 1))


# intercepts
int_df <- data.frame(intercept = as.numeric(intercepts),
                     veg = rep(veg_levels, each = nrow(intercepts)))
boxplot(intercept ~ veg, int_df)



# DHARMa (ndvi_mean, pp) ---------------------------------------------

# topographic logit
topo_logit <- b_elev[best_iter] * (dfit$elevation - optim_elev[best_iter]) ^ 2 +
  b_tpi[best_iter] *  (dfit$TPI1k - optim_tpi[best_iter]) ^ 2 +
  b_north[best_iter] * dfit$northing +
  b_slope[best_iter] * dfit$slope

# veg logits
veg_alpha <- sdata$veg %*% intercepts[best_iter, ]
ndvi_rep <- matrix(rep(dfit$ndvi_mean, V), nrow(dfit), V)
optims_rep <- outer(rep(1, nrow(dfit)), optim_ndvi[best_iter, ])
ndvi_quad <- (ndvi_rep - optims_rep) ^ 2 * sdata$veg
ndvi_logit <- ndvi_quad %*% b_ndvi[best_iter, ]

# pp
pp_logit <- b_pp[best_iter] *  (dfit$pp - optim_pp[best_iter]) ^ 2

# fitted p
fit <- plogis(veg_alpha + ndvi_logit + topo_logit + pp_logit) %>% as.numeric()
hist(fit)
mean(fit)
mean(dfit$burned)

# Con la parametrización del simpler anda bien, pero con el otro los intercepts 
# estaban mal

dfit$pfit <- fit

ysim <- matrix(rbinom(nrow(dfit) * S, size = 1, prob = fit),
               nrow(dfit), S)

res <- createDHARMa(ysim, dfit$burned, fittedPredictedResponse = fit)
plot(res)

# da hermosoooo
dfit$res <- res$scaledResiduals

ggplot(dfit, aes(ndvi_mean, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(pp, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(elevation, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(TPI1k, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(slope, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

ggplot(dfit, aes(northing, res)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "gam", method.args = list(family = betar())) +
  facet_wrap(vars(vegetation_class), scales = "free_x")

# R2
var(fit) / (var(fit) + mean(fit * (1 - fit))) # 0.1816227 



# Exports -----------------------------------------------------------------

b_ndvi_export <- b_ndvi[best_iter, ]
optim_ndvi_export <- optim_ndvi[best_iter, ]
intercept_export <- intercepts[best_iter, ]

names(b_ndvi_export) <- names(optim_ndvi_export) <- names(intercept_export) <- veg_levels

b_elev[best_iter]; optim_elev[best_iter]
b_tpi[best_iter]; optim_tpi[best_iter]

b_north[best_iter]
b_slope[best_iter]

intercepts[best_iter, ]

theta_hat <- list(
  intercept = intercept_export,
  
  b_ndvi = b_ndvi_export,
  optim_ndvi = optim_ndvi_export,
  
  b_pp = b_pp[best_iter], 
  optim_pp = optim_pp[best_iter],
  
  b_elev = b_elev[best_iter], 
  optim_elev = optim_elev[best_iter],
  
  b_tpi = b_tpi[best_iter], 
  optim_tpi = optim_tpi[best_iter],
  
  b_north = b_north[best_iter],
  b_slope = b_slope[best_iter],
  
  notes = 
"MAP estimate for a regional logistic regression.
NDVI is the mean of the summers maximums.
See <landscape_flammability.R> for details.
Predictors are at their original scale."
)

cat(theta_hat$notes)
# saveRDS(theta_hat, "landscape flammability indices/landscape_flammability_model_point_estimates.rds")

# To standardize the variables, it was computed in the whole burnable area:

# vfi_mean: -0.9551693163497065
# vfi_sd:    1.187131140133366
# tfi_mean: -1.0654265335115882
# tfi_sd:    1.1019703840457942

theta_hat <- readRDS("landscape_flammability_model_point_estimates.rds")
theta_hat



# TO DO -------------------------------------------------------------------

# clean the code. Use only the ndvi_mean model. ndvi_max fits similarly, but
# the ndvi_max time series is noisy because it increases since 2014, when more
# satellites were added. In turn, the mean is not so sensitive to the higher 
# image availability. The ndvi_mean varies more in space than in time, so it's
# reasonable to consider it a spatial predictor. 

# Use only the "simpler" parameterization for the stan model, and then delete
# the other. It takes twice the time as the other, but it gets the correct intercept. 
# Centering quadratic terms is not easy.

