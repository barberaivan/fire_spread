# To increase the FWI temporal resolution we need to define a temporal window
# in which it is relevant for fire spread, e.g., one week, one month, or 3 months
# before the ignition. However, a more countinuos approach is to consider that the
# FWI effect on fire spread decreases as a function of the temporal lag. Hence, 
# we can fit a cumulative effects model. 

# But the fires dates are uncertain, and some are not so uncertain but spread
# over long time periods (1 or 2 months in large fires), so defining a single
# date in which the fire spread is hard. In an ideal world, we would know the
# daily burned area by fire, but that's not happening. So, for simplicity, 
# we could assume that most of the fire spread in one day, and treat the day
# as an extra parameter to estimate (maybe to be marginalized). Maybe we could 
# also use something like poisson additivity (see Marchal et al. papers), but 
# that would make sense with predictors that affect fire spread in a much smaller
# temporal scale, like wind speed. If we suppose that FWI affects fire spread 
# in a longer time frame, that doesn't make sense. 

# We model the fire size as a function of each fire's FWI. We assume a Gamma
# distribution, with link log, where the log-mean fire size (mu) depends on the
# fwi in the previous 90 days:

# log(mu[d]) = alpha + sum_{l = 0}^L b[l] * FWI[d-l]
# b[l] = b0 * exp(- l ^ 2 / scale)

# where d is the day, L = 90 and l is the lag index in days.
# alpha, b0 and scale are parameters to estimate. The larger the scale, 
# the more lasting is the FWI effect. Bear in mind, for interpretation, that
# FWI is already a cumulative variable. However, it presents considerable 
# daily variation.

# For a faster decline in the FWI effect, it would be sensible to compare this
# model with an exponential decay:
# b[l] = b0 * rho ^ (-l)

# We will fit the model in Stan, and will consider spread date uncertainty.
# The fires database records dale_lower and date_upper, obtained by MODIS or 
# Landsat. In the Landsat case, 
# spread_lower = date_lower + 1, and
# spread_upper = date_upper - 1, 
# because with Landsat we recorded the last day when the NBR was high (lower) 
# and the first day when it was already low (upper).
# With MODIS we recorded the first-last date with heat signal, so it's the proper
# lower or upper date.

# To consider spread date uncertainty we need the FWI time series in 
# [spread_lower - 90, spread_upper]
# for the centroid in every fire. To use a data matrix equal for all fires
# we will select the most uncertain usable fire and get the same range length
# of dates for all. In many cases, we will not use all the data.
# (This is a programmatic issue related to the lack of rugged data structures
# in Stan).

# Thomas will improve the date resolution of the fires he knows.


# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(terra)
library(DHARMa)
library(lubridate)
library(rstan)
library(mgcv)

# FWI raster ts -----------------------------------------------------------

# r <- rast("data/fwi_daily_1998-2022/fwi_daily_19980101_20230630.tif")
# v <- values(r)
# aa <- acf(v[1, ], lag.max = 30)
 
# tt <- ts(v[1, ], start = 1998, end = 2023, frequency = 365)
# plot(tt)
# dd <- decompose(tt)
# plot(dd)

# dd2 <- decompose(tt, "multiplicative")
# plot(dd2)

# It has considerable temp correlation


# Functions ---------------------------------------------------------------

# mean, median and maximum
mmm <- function(x) {
  return(
    c("mean" = mean(x),
      "median" = median(x),
      "max" = max(x),
      "days_focal" = length(x))
  )
}

normalize <- function(x) x / sum(x)

# Data --------------------------------------------------------------------

fires_vec <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp")
fires <- as.data.frame(fires_vec)
# remove fires with uncertain year
uncertain <- which(fires$obs == "uncertain_year")
fires <- fires[-uncertain, ]

fwi_data <- read.csv("data/climatic_data_by_fire_fwi.csv")

# Turn dates into dates
dd <- c("date_start", "date_end", "date_l")
for(d in dd) fwi_data[, d] <- as.Date(fwi_data[, d])

dd <- c("date", "date_l", "date_u")
for(d in dd) fires[, d] <- as.Date(fires[, d])

# Correct date_l and date_u for Landsat-based images, so they are included in the
# focal period. For landsat, date_l and date_u are the dates between which the
# fire occurred, but not included. For modis, they are inclusive. Turn all into
# inclusive. Do this only for date_diff larger than 2

ll <- (fires$datesrc %in% c("landsat", "Mermoz")) &
      (fires$date_u - fires$date_l > 2)  
fires$date_l[ll] <- fires$date_l[ll] + 1
fires$date_u[ll] <- fires$date_u[ll] - 1


# add day of year, centred at July first
yytemp <- yday(fires$date)
# yytemp <- yday(seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = 1)) # test
yy2 <- yytemp
midday <- yday(as.Date("2019-07-01"))
maxday <- yday(as.Date("2019-12-31"))
yy2[yytemp >= midday] <- yytemp[yytemp >= midday] - midday + 1
yy2[yytemp < midday] <- yytemp[yytemp < midday] + (maxday - midday) + 1

fires$doy <- yy2

# fire size as a function of FWI summary (mean, median, max) --------------

# Use mean, median, maximum of the period?
fwi_focal <- do.call("rbind", lapply(fires$fire_id, function(f) {
  # f = fires$fire_id[1]
  
  # get fwi values
  row <- which(fires$fire_id == f)
  dtemp <- fwi_data[fwi_data$fire_id == f, , drop = F]
  fwi_vals <- dtemp[1, grep("fwi", names(fwi_data))] %>% as.numeric
  # subset focal dates
  date_seq_av <- seq(dtemp$date_start, dtemp$date_end, 1)
  date_seq_focal <- seq(fires$date_l[row], fires$date_u[row], 1)
  
  if(length(date_seq_av) != length(fwi_vals)) stop("Date sequence does not match FWI values.")
  
  days_use <- which(date_seq_av %in% date_seq_focal)
  vvv <- mmm(fwi_vals[days_use])
  return(vvv)
})) %>% as.data.frame()

barplot(table(fwi_focal$days_focal)) # most are 1s

# to subset modis fires, use this:
modis_ids <- grep("modis", fires$datesrc)

# merge and longanize
fires_clim <- cbind(fires, fwi_focal)

fires_clim_long <- pivot_longer(
  fires_clim, all_of(which(names(fires_clim) %in% c("mean", "median", "max"))),
  names_to = "fwi_summ", values_to = "fwi"
)

fires_clim_long <- fires_clim_long[fires_clim_long$days_focal < 20, ]

# plots
ggplot(fires_clim_long, aes(x = fwi, y = area_ha)) +
  geom_smooth(method = "glm", method.args = list(family = Gamma(link = "log")),
              se = F) +
  geom_point() +
  facet_wrap(vars(fwi_summ), nrow = 2)

ggplot(fires_clim_long, aes(x = fwi, y = area_ha)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_wrap(vars(fwi_summ), ncol = 2)

ggplot(fires_clim_long, aes(x = fwi, y = area_ha)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_wrap(vars(fwi_summ), ncol = 2) + 
  scale_y_continuous(trans = "log10")


# check 3 models at log scale with DHARMa
m1 <- lm(log(area_ha) ~ mean, fires_clim)
r1 <- simulateResiduals(m1)
plot(r1)
plotResiduals(r1, form = fires_clim$doy, rank = F)

m2 <- lm(log(area_ha) ~ median, fires_clim)
r2 <- simulateResiduals(m2)
plot(r2)
plotResiduals(r2, form = fires_clim$doy, rank = F)

m3 <- lm(log(area_ha) ~ max, fires_clim)
r3 <- simulateResiduals(m3)
plot(r3)
plotResiduals(r3, form = fires_clim$doy, rank = F)

# there seems to be no trend with day of year.

# The area-FWI relationship is noisy when the during-fire FWI is considered, 
# and this doesn't change if we use the modis dataset or the whole dataset.
# Even if we consider only the fires that burned in 1 or 2 days, the noise is
# high. So, trying to estimate the burning date would not make sense, since 
# fire can burn in more than one day and the relationship with fire area is
# weak.
# Try cumulative models, combined with mean and max FWI in the fire-period.

# Cumulative models -------------------------------------------------------

# Assume a log-normal distribution for the burned area by fire. At the log scale,
# the predictors will be
# fwi_focal: mean or max at the focal period, and
# fwi_previous: weighted average of FWI in the previous 120 days, with an estimated
#   weighting function (expquad and exp).
# Probably a FWI effect over sigma will be needed.

# get fwi_previous, in days -1 : -120, where 0 is date_l, the lower limit for
# fire occurrence.

fwi_previous <- do.call("rbind", lapply(fires$fire_id, function(f) {
  # f = fires$fire_id[1]
  
  # get fwi values
  row <- which(fires$fire_id == f)
  dtemp <- fwi_data[fwi_data$fire_id == f, , drop = F]
  fwi_vals <- dtemp[1, grep("fwi", names(fwi_data))] %>% as.numeric
  # subset previous dates
  date_seq_av <- seq(dtemp$date_start, dtemp$date_end, 1)
  date_seq_prev <- seq(fires$date_l[row] - 120, fires$date_l[row] - 1, 1)
  
  if(length(date_seq_av) != length(fwi_vals)) stop("Date sequence does not match FWI values.")
  days_use <- which(date_seq_av %in% date_seq_prev)
  
  dd <- data.frame(fwi = fwi_vals, date = date_seq_av)
  dd <- dd[order(dd$date, decreasing = TRUE), ]
  dd_use <- dd[dd$date %in% date_seq_prev, ]
  
  m <- matrix(dd_use$fwi, 1)
  colnames(m) <- as.character(seq(-1, -120, -1))
  
  return(m)
}))

rownames(fwi_previous) <- fires$fire



# Prior check -------------------------------------------------------------

# expquad model
ls <- 20
curve(exp(-(x / ls) ^ 2), from = 0, to = 120)
for(i in 1:500) {
  ls <- abs(rnorm(1, 0, sd = 150))
  # ls <- runif(1, 0, 500)
  curve(exp(-(x / ls) ^ 2), add = TRUE, col = rgb(0, 0, 0, 0.05))
}
# ls ~ normal(0, 150) allows for no decrease in importance

# exp model
rho <- 0.3
unit <- 15
curve(rho ^ (x / unit), from = 0, to = 120, ylim = c(0, 1))
for(i in 1:500) {
  rho <- runif(1, 0, 1)
  curve(rho ^ (x / unit), add = TRUE, col = rgb(0, 0, 0, 0.05))
}
# rho ~ unif(0, 1) allows for no decrease in importance, with more likelihood of
#                  fast decreases, when unit = 15





# Data for Stan ------------------------------------------------------
rows_fit <- fires_clim$days_focal < 20

# center fwi_previous based on ls = 20 or rho = 0.3 (unit = 15)
days_seq <- 0:(120-1)
fwi_previous_center_exp <- mean(fwi_previous[rows_fit, ] %*% 
                                  normalize(0.3 ^ (days_seq / 15)))
fwi_previous_center_expquad <- mean(fwi_previous[rows_fit, ] %*% 
                                      normalize(exp(-(days_seq / 20) ^ 2)))

sdata <- list(
  y = log(fires$area_ha[rows_fit]),
  fwi_previous = fwi_previous[rows_fit, ],
  N_days = 120,
  days_distance = 0:(120-1), # day zero for previous period is day -1 wrt date_l
  N = sum(rows_fit),
  
  fwi_focal_mean = fires_clim$mean[rows_fit],
  fwi_focal_max = fires_clim$max[rows_fit],
  
  rho_reference = 0.3,
  ls_reference = 20, 
  use_mean = 1,
  unit = 15, # days unit for exponential model
  
  prior_sd_int = 20, 
  prior_mean_int = log(1000),
  prior_sd_b = 10,
  prior_sd_ls = 150,
  prior_sd_sigma = 8
)

# Compile both stan models
model_expquad <- stan_model("temporal scale for climatic predictors/FWI model expquad.stan")

# Model 1: mean fwi, expquad ----------------------------------------------

sdata$use_mean == 0
m1 <- sampling(
  model_expquad, data = sdata, refresh = 100, 
  chains = 4, 
  cores = 4, 
  iter = 2000,
  pars = c("alpha", "intercept", "b_focal", "b_previous", "ls", "sigma", "mu")
)
# una luz.
sm1 <- summary(m1)[[1]]

# para chequear que mu y mu2 son lo mismo hay que pedirle que guarde a mu2
# mu1 <- as.matrix(m1, "mu")
# mu2 <- as.matrix(m1, "mu2")
# hist( as.vector(abs(mu1 - mu2)) )
# son lo mismo, pero R tiene más precisión que Stan.

pairs(m1, pars = c("alpha", "b_focal", "b_previous", "ls", "sigma"))

# legthscale
ls_hat <- as.matrix(m1, "ls") %>% as.numeric
npost <- length(ls_hat)
curve(exp(-(x / ls_hat[1]) ^ 2), from = 0, to = 120, col = rgb(0, 0, 0, 0.05),
      ylim = c(0, 1))
for(i in 2:npost) {
  ls <- ls_hat[i]
  # ls <- runif(1, 0, 500)
  curve(exp(-(x / ls) ^ 2), add = TRUE, col = rgb(0, 0, 0, 0.05))
}

# prior-posterior
plot(density(ls_hat, from = 0))
curve(dnorm(x, 0, sdata$prior_sd_ls) * 2, add = TRUE, col = 2)
# aprende poco.

# DHARMa
mu_hat <- as.matrix(m1, "mu") %>% t
sigma_hat <- as.matrix(m1, "sigma") %>% as.numeric
ysim <- matrix(rnorm(sdata$N * npost), sdata$N, npost) * sigma_hat + mu_hat

# en log
m1_res <- createDHARMa(ysim, sdata$y)
plot(m1_res)

# en exp
m1_res_exp <- createDHARMa(exp(ysim), exp(sdata$y))
plot(m1_res_exp)

# se ve lo mismo que sin FWI acumulado.

# ajuste:

# en log
mu_var <- apply(mu_hat, 2, var)
r2log <- mu_var / (mu_var + sigma_hat ^ 2)
plot(density(r2log, from = 0, to = 1)) # re bajo
abline(v = mean(r2log), lty = 2, col = 2)

# en exp
mu_hat_exp <- exp(mu_hat + 0.5 * sigma_hat ^ 2)
mu_var_exp <- apply(mu_hat_exp, 2, var)
sigma2_exp <- (exp(sigma_hat ^ 2) - 1) * exp(2 * mu_hat + sigma_hat ^ 2)
r2exp <- mu_var_exp / (mu_var_exp + sigma2_exp)
plot(density(r2exp, from = 0, to = 1)) # re bajo
abline(v = mean(r2exp), lty = 2, col = 2)

# comparar con modelo del máximo, y también ver los pesos de los datos, a
# ver si alguno anda muy mal.
plot(sdata$y ~ rowMeans(mu_hat)); abline(0, 1)

# usando fwi_mean y fwi_max da lo mismo, quizás en parte porque 
# son valores muy parecidos. 
plot(fires_clim$mean[rows_fit] ~ fires_clim$max[rows_fit]); abline(0, 1)

# Y si usamos el acumulado hasta date?


# FWI previous from central date ------------------------------------------

fires_sub <- fires[fires_clim$days_focal < 20, ]

fwi_previous_mid <- do.call("rbind", lapply(fires_sub$fire_id, function(f) {
  # f = fires_sub$fire_id[1]
  
  # get fwi values
  row <- which(fires_sub$fire_id == f)
  # print(row)
  
  dtemp <- fwi_data[fwi_data$fire_id == f, , drop = F]
  fwi_vals <- dtemp[1, grep("fwi", names(fwi_data))] %>% as.numeric
  # subset previous dates
  date_seq_av <- seq(dtemp$date_start, dtemp$date_end, 1)
  date_seq_prev <- seq(fires_sub$date[row] - 120, fires_sub$date[row], 1)
  
  if(length(date_seq_av) != length(fwi_vals)) stop("Date sequence does not match FWI values.")
  days_use <- which(date_seq_av %in% date_seq_prev)
  
  dd <- data.frame(fwi = fwi_vals, date = date_seq_av)
  dd <- dd[order(dd$date, decreasing = TRUE), ]
  dd_use <- dd[dd$date %in% date_seq_prev, ]
  
  m <- matrix(dd_use$fwi, 1)
  colnames(m) <- as.character(seq(0, -120, -1))
  
  return(m)
}))
rownames(fwi_previous_mid) <- fires_sub$fire_id

# Data for Stan ------------------------------------------------------
rows_fit <- fires_clim$days_focal < 20

# center fwi_previous based on ls = 20 or rho = 0.3 (unit = 15)
days_seq <- 0:(120)
fff <- fwi_previous_mid %*% normalize(exp(-(days_seq / 20) ^ 2))
fwi_center_expquad_mid <- mean(fff)
fwi_scale_expquad_mid <- sd(fff)

sdata <- list(
  y = log(fires$area_ha[rows_fit]),
  fwi_previous = fwi_previous_mid,
  N_days = 121,
  days_distance = 0:(120), # day zero for previous period is middate (fires$date)
  N = sum(rows_fit),
  
  fwi_center = fwi_center_expquad_mid,
  fwi_scale = fwi_scale_expquad_mid,
  
  use_mean = 1,
  unit = 15, # days unit for exponential model
  
  prior_sd_int = 20, 
  prior_mean_int = log(1000),
  prior_sd_b = 10,
  prior_sd_ls = 150,
  prior_sd_sigma = 8
)

# Compile both stan models
model_expquad_mid <- stan_model("temporal scale for climatic predictors/FWI model expquad middate.stan")

# Model 2: expquad cumulative from mid date --------------------------------

m2 <- sampling(
  model_expquad_mid, data = sdata, refresh = 100, 
  chains = 4, 
  cores = 4, 
  iter = 2000,
  pars = c("alpha", "intercept", "b", "ls", "sigma", "mu", "b_raw")
)
# una luz.
sm2 <- summary(m2)[[1]]

pairs(m2, pars = c("alpha", "b", "ls", "sigma"))

# legthscale
ls_hat <- as.matrix(m2, "ls") %>% as.numeric
npost <- length(ls_hat)
curve(exp(-(x / ls_hat[1]) ^ 2), from = 0, to = 120, col = rgb(0, 0, 0, 0.05),
      ylim = c(0, 1))
for(i in 2:npost) {
  ls <- ls_hat[i]
  # ls <- runif(1, 0, 500)
  curve(exp(-(x / ls) ^ 2), add = TRUE, col = rgb(0, 0, 0, 0.05))
}

# prior-posterior
plot(density(ls_hat, from = 0))
curve(dnorm(x, 0, sdata$prior_sd_ls) * 2, add = TRUE, col = 2)
# aprende un poco más que antes eh

# DHARMa
mu_hat <- as.matrix(m2, "mu") %>% t
sigma_hat <- as.matrix(m2, "sigma") %>% as.numeric
ysim <- matrix(rnorm(sdata$N * npost), sdata$N, npost) * sigma_hat + mu_hat

# en log
m2_res <- createDHARMa(ysim, sdata$y)
plot(m2_res)

# en exp
m2_res_exp <- createDHARMa(exp(ysim), exp(sdata$y))
plot(m2_res_exp, rank = F)

# se ve lo mismo que sin FWI acumulado.

# ajuste:

# en log
mu_var <- apply(mu_hat, 2, var)
r2log <- mu_var / (mu_var + sigma_hat ^ 2)
plot(density(r2log, from = 0, to = 1)) # re bajo
abline(v = mean(r2log), lty = 2, col = 2)

# en exp
mu_hat_exp <- exp(mu_hat + 0.5 * sigma_hat ^ 2)
mu_var_exp <- apply(mu_hat_exp, 2, var)
sigma2_exp <- (exp(sigma_hat ^ 2) - 1) * exp(2 * mu_hat + sigma_hat ^ 2)
r2exp <- mu_var_exp / (mu_var_exp + sigma2_exp)
plot(density(r2exp, from = 0, to = 1)) # re bajo
abline(v = mean(r2exp), lty = 2, col = 2)

# comparar con modelo del máximo, y también ver los pesos de los datos, a
# ver si alguno anda muy mal.
plot(sdata$y ~ rowMeans(mu_hat)); abline(0, 1)

# Ajusta un poquitito peor, pero todo muy parecido. De esta forma sí puede 
# estimar el lscale.


# Reflexiones: mucho ruido ------------------------------------------------

# Dado lo ruidosa que es la relación, quizás convenga usar promedios seamanales,
# quincenales o mensuales de FWI. Eso nos evita tener que simular diariamente
# luego. Entonces, nos basamos en la fecha central y calculamos la media
# de la semana, de dos semanas o de el mes.

# get mean aggregated over week, fortnite and month
fwi_agg <- do.call("rbind", lapply(fires$fire_id, function(f) {
  # f = fires$fire_id[51]
  
  # get fwi values
  row <- which(fires$fire_id == f)
  dtemp <- fwi_data[fwi_data$fire_id == f, , drop = F]
  fwi_vals <- dtemp[1, grep("fwi", names(fwi_data))] %>% as.numeric
  # subset focal dates
  date_seq_av <- seq(dtemp$date_start, dtemp$date_end, 1)
  middate <- fires$date[row]
  
  if(length(date_seq_av) != length(fwi_vals)) stop("Date sequence does not match FWI values.")
  
  dd <- data.frame(fwi = fwi_vals, 
                   date = date_seq_av,
                   week = week(date_seq_av),
                   month = month(date_seq_av),
                   year = year(date_seq_av))
  
  dd <- dd[order(dd$date), ]
  
  # if there is a year change in the period, restart weeks backwards from 
  # december 31st. This is because the 7-day counting may create a last week 
  # with too few days.
  if(max(diff(dd$year)) == 1) {
    restart <- dd$year == min(dd$year)  
    weeks_n <- ceiling(sum(restart) / 7)
    week_tmp <- rep(seq(0, -(weeks_n - 1)), each = 7)[1:sum(restart)]
    week_use <- rev(week_tmp)
    dd$week[restart] <- week_use
  }
  
  # view(dd)
  # create fortnight, centered at the new year
  odd_weeks <- (dd$week %% 2) == 0
  dd$fortnight <- dd$week
  dd$fortnight[odd_weeks & dd$week > 0] <- dd$week[odd_weeks & dd$week > 0] - 1
  dd$fortnight[!odd_weeks & dd$week < 0] <- dd$week[!odd_weeks & dd$week < 0] + 1

  dd_mid <- dd[dd$date == middate, , drop = F]
  vals <- c(
    "fwi_week" = max(dd$fwi[dd$week == dd_mid$week]),
    "fwi_fort" = max(dd$fwi[dd$fortnight == dd_mid$fortnight]),
    "fwi_month" = max(dd$fwi[dd$month == dd_mid$month])
  )
  
  return(vals)
})) %>% as.data.frame()

# the three summaries are highly correlated.

# merge and longanize
fires_agg <- cbind(fires, fwi_agg)

fires_agg_long <- pivot_longer(
  fires_agg, all_of(which(names(fires_agg) %in% c("fwi_week", "fwi_fort", "fwi_month"))),
  names_to = "fwi_agg", values_to = "fwi"
)

fires_not_uncertain <- fires_clim$fire_id[fires_clim$days_focal < 5] 
fires_agg_long <- fires_agg_long[fires_agg_long$fire_id %in% fires_not_uncertain, ]

# plots
ggplot(fires_agg_long, aes(x = fwi, y = area_ha)) +
  geom_smooth(method = "glm", method.args = list(family = Gamma(link = "log")),
              se = F) +
  geom_point() +
  facet_wrap(vars(fwi_agg), nrow = 2)

ggplot(fires_agg_long, aes(x = fwi, y = area_ha)) +
  geom_smooth(method = "lm") +
  geom_point() +
  facet_wrap(vars(fwi_agg), ncol = 2) + 
  scale_y_continuous(trans = "log")




# check 3 models at log scale with DHARMa

fires_agg_ok <- fires_agg[fires_agg$fire_id %in% fires_not_uncertain,]

# location-scale gaussian model at log scale

# week
m1 <- gam(list(log(area_ha) ~ fwi_week, ~ log(fwi_week)), data = fires_agg_ok,
          family = gaulss())
r1 <- createDHARMa(simulatedResponse = simulate(m1, nsim = 3000),
                   observedResponse = log(fires_agg_ok$area_ha),
                   fittedPredictedResponse = fitted(m1)[, 1])
# fort
m2 <- gam(list(log(area_ha) ~ fwi_fort, ~ log(fwi_fort)), data = fires_agg_ok,
          family = gaulss())
r2 <- createDHARMa(simulatedResponse = simulate(m2, nsim = 3000),
                   observedResponse = log(fires_agg_ok$area_ha),
                   fittedPredictedResponse = fitted(m2)[, 1])
# month
m3 <- gam(list(log(area_ha) ~ fwi_month, ~ log(fwi_month)), data = fires_agg_ok,
          family = gaulss())
r3 <- createDHARMa(simulatedResponse = simulate(m3, nsim = 3000),
                   observedResponse = log(fires_agg_ok$area_ha),
                   fittedPredictedResponse = fitted(m3)[, 1])


# now repeat the previous models, using mean, median and max at each fires' range

fires_clim_ok <- fires_clim[fires_clim$fire_id %in% fires_not_uncertain, ]

# mean
m4 <- gam(list(log(area_ha) ~ mean, ~ log(mean)), data = fires_clim_ok,
          family = gaulss())
r4 <- createDHARMa(simulatedResponse = simulate(m4, nsim = 3000),
                   observedResponse = log(fires_clim_ok$area_ha),
                   fittedPredictedResponse = fitted(m4)[, 1])

# median
m5 <- gam(list(log(area_ha) ~ median, ~ log(median)), data = fires_clim_ok,
          family = gaulss())
r5 <- createDHARMa(simulatedResponse = simulate(m5, nsim = 3000),
                   observedResponse = log(fires_clim_ok$area_ha),
                   fittedPredictedResponse = fitted(m5)[, 1])

# max
m6 <- gam(list(log(area_ha) ~ max, ~ log(max)), data = fires_clim_ok,
          family = gaulss())
r6 <- createDHARMa(simulatedResponse = simulate(m6, nsim = 3000),
                   observedResponse = log(fires_clim_ok$area_ha),
                   fittedPredictedResponse = fitted(m6)[, 1])




plot(r1, main = "week", rank = F)
plot(r2, main = "fort", rank = F)
plot(r3, main = "month", rank = F)

plot(r4, main = "fwi mean", rank = F)
plot(r5, main = "fwi median", rank = F)
plot(r6, main = "fwi max", rank = F)

# with good fires as <20 days uncertainty:

                     # mean in period // max in period
summary(m1)$dev.expl # 0.1752908      // 0.225506       | week
summary(m2)$dev.expl # 0.2242881      // 0.3039021      | fort
summary(m3)$dev.expl # 0.2053487      // 0.2801199      | month

                     # mean, median or max in the supposed period
summary(m4)$dev.expl # 0.1804293  | mean
summary(m5)$dev.expl # 0.1819342  | median
summary(m6)$dev.expl # 0.3100527  | max

AIC(m1, m2, m3, m4, m5, m6)

# Models are very similar. In all cases, predicting sigma with log(fwi) helps a lot.
# The best according to r2 is fortnight model, but according to AIC is month model.

# If we use the max fwin in period (week, fort, month), the r2 increase
# fortnight model seems to be the best.

# try using log(fwi as predictor)


# Fortnight max cumulative models -----------------------------------------

# In this case, the focal fortnight is going to have the maximum weight, 
# and it will decrease with the delay. 
# Compare 3 models: no weight, weighted exp and weighted expquad. 

# for all this, keep only fires with <= 20 days date uncertainty.
fires_sub <- fires[fires_clim$days_focal <= 20, ]

# get mean aggregated over week, fortnite and month
fwi_fort_list <- lapply(fires_sub$fire_id, function(f) {
  # f = fires$fire_id[1]
  
  # get fwi values
  row <- which(fires$fire_id == f)
  dtemp <- fwi_data[fwi_data$fire_id == f, , drop = F]
  fwi_vals <- dtemp[1, grep("fwi", names(fwi_data))] %>% as.numeric
  # subset focal dates
  date_seq_av <- seq(dtemp$date_start, dtemp$date_end, 1)
  middate <- fires$date[row]
  
  if(length(date_seq_av) != length(fwi_vals)) stop("Date sequence does not match FWI values.")
  
  dd <- data.frame(fwi = fwi_vals, 
                   date = date_seq_av,
                   week = week(date_seq_av),
                   year = year(date_seq_av))
  
  dd <- dd[order(dd$date), ]
  
  # if there is a year change in the period, restart weeks backwards from 
  # december 31st. This is because the 7-day counting may create a last week 
  # with too few days.
  if(max(diff(dd$year)) == 1) {
    restart <- dd$year == min(dd$year)  
    weeks_n <- ceiling(sum(restart) / 7)
    week_tmp <- rep(seq(0, -(weeks_n - 1)), each = 7)[1:sum(restart)]
    week_use <- rev(week_tmp)
    dd$week[restart] <- week_use
  }
  
  # view(dd)
  # create fortnight, centered at the new year
  odd_weeks <- (dd$week %% 2) == 0
  dd$fortnight <- dd$week
  dd$fortnight[odd_weeks & dd$week > 0] <- dd$week[odd_weeks & dd$week > 0] - 1
  dd$fortnight[!odd_weeks & dd$week < 0] <- dd$week[!odd_weeks & dd$week < 0] + 1
  
  dd_mid <- dd[dd$date == middate, , drop = F]
  
  # subset dd using only the focal and previous fortnights
  dd_sub <- dd[dd$fortnight <= dd_mid$fortnight, ]
  
  # aggregate fwi max by fortnight
  dd_agg <- aggregate(fwi ~ fortnight, dd_sub, max)
  dd_agg$fortnight <- -(nrow(dd_agg) - 1):0
  
  return(dd_agg)
})
min_forts <- lapply(fwi_fort_list, nrow) %>% unlist %>% min

fwi_fort <- do.call("rbind", lapply(fwi_fort_list, function(x) {
  # x <- fwi_fort_list[[1]]
  xx <- x[order(x$fortnight, decreasing = TRUE), ]
  fw <- xx$fwi[1:min_forts]
  mm <- matrix(fw, 1)
  colnames(mm) <- xx$fortnight[1:min_forts]
  return(mm)
}))

# Prior check fortnight ---------------------------------------------------

# expquad model
ls <- 1
curve(exp(-(x / ls) ^ 2), from = 0, to = 9)
for(i in 1:500) {
  ls <- abs(rnorm(1, 0, sd = 9))
  # ls <- runif(1, 0, 500)
  curve(exp(-(x / ls) ^ 2), add = TRUE, col = rgb(0, 0, 0, 0.05))
}
# ls ~ normal(0, 6) allows for almost no decrease in importance

# exp model
rho <- 0.3
unit <- 1
curve(rho ^ (x / unit), from = 0, to = 9, ylim = c(0, 1))
for(i in 1:500) {
  rho <- runif(1, 0, 1)
  curve(rho ^ (x / unit), add = TRUE, col = rgb(0, 0, 0, 0.05))
}
# rho ~ unif(0, 1) allows for no decrease in importance, with more likelihood of
#                  fast decreases, when unit = 1



# Data for Stan (fortnight) ------------------------------------------------s

# center fwi_fort based on ls = 1 or rho = 0.3 (unit = 1)
fseq <- 0:(ncol(fwi_fort)-1)


fwi_fort_center_exp <- mean(fwi_fort %*% normalize(0.3 ^ fseq))
fwi_fort_center_expquad <- mean(fwi_fort %*% normalize(exp(-(fseq / 1) ^ 2)))

fwi_fort_scale_exp <- sd(fwi_fort %*% normalize(0.3 ^ fseq))
fwi_fort_scale_expquad <- sd(fwi_fort %*% normalize(exp(-(fseq / 1) ^ 2)))

# the same in log, for sigma
fwi_fort_scale_exp_log <- sd(log(fwi_fort %*% normalize(0.3 ^ fseq)))
fwi_fort_scale_expquad_log <- sd(log(fwi_fort %*% normalize(exp(-(fseq / 1) ^ 2))))

sdata_exp <- list(
  y = log(fires_sub$area_ha),
  fwi_fort = fwi_fort,
  N_fort = min_forts,
  N = nrow(fires_sub),
  
  fwi_center = fwi_fort_center_exp,
  fwi_scale = fwi_fort_scale_exp,
  fwi_scale_log = fwi_fort_scale_exp_log, # not center because cant log a negative
  
  prior_mean_int_mu = log(1000),
  prior_mean_int_sigma = log(2), # sd(log(fires_sub$area_ha))
  
  prior_sd_int_mu = 20,
  prior_sd_int_sigma = 10, # sd(log(fires_sub$area_ha))
  
  prior_sd_b_mu = 10,
  prior_sd_b_sigma = 2,
  
  prior_sd_ls = 6,
  prior_sd_sigma = 8
)

sdata_expquad <- list(
  y = log(fires_sub$area_ha),
  fwi_fort = fwi_fort,
  N_fort = min_forts,
  N = nrow(fires_sub),
  
  fwi_center = fwi_fort_center_expquad,
  fwi_scale = fwi_fort_scale_expquad,
  fwi_scale_log = fwi_fort_scale_exp_log, # not center because cant log a negative
  
  prior_mean_int_mu = log(1000),
  prior_mean_int_sigma = log(2), # sd(log(fires_sub$area_ha))
  
  prior_sd_int_mu = 20,
  prior_sd_int_sigma = 10, # sd(log(fires_sub$area_ha))
  
  prior_sd_b_mu = 10,
  prior_sd_b_sigma = 2,
  
  prior_sd_ls = 6,
  prior_sd_sigma = 8
)

sdata_focal <- list(
  y = log(fires_sub$area_ha),
  fwi_fort = fwi_fort[, 1],
  N = nrow(fires_sub),
  
  prior_mean_int_mu = log(1000),
  prior_mean_int_sigma = log(2), # sd(log(fires_sub$area_ha))
  
  prior_sd_int_mu = 20,
  prior_sd_int_sigma = 10, # sd(log(fires_sub$area_ha))
  
  prior_sd_b_mu = 10,
  prior_sd_b_sigma = 2,
  
  prior_sd_ls = 6,
  prior_sd_sigma = 8
)

# Compile stan models
model_expquad <- stan_model("temporal scale for climatic predictors/FWI model fort expquad.stan") # m1
model_exp <- stan_model("temporal scale for climatic predictors/FWI model fort exp.stan")         # m2
model_focal <- stan_model("temporal scale for climatic predictors/FWI model fort focal.stan")     # m3

m1 <- sampling(
  model_expquad, sdata_expquad, refresh = 100, 
  cores = 4, chains = 4, iter = 2000,
  pars =  c("alpha_mu", "alpha_sigma", 
            "intercept_mu", "intercept_sigma", 
            "b_mu", "b_sigma", "ls", 
            "sigma", "mu")
)

m2 <- sampling(
  model_exp, sdata_exp, refresh = 100, 
  cores = 4, chains = 4, iter = 2000,
  pars =  c("alpha_mu", "alpha_sigma", 
            "intercept_mu", "intercept_sigma", 
            "b_mu", "b_sigma", "rho", 
            "sigma", "mu")
)

m3 <- sampling(
  model_focal, sdata_focal, refresh = 100, 
  cores = 4, chains = 4, iter = 2000,
  pars =  c("alpha_mu", "alpha_sigma", 
            "intercept_mu", "intercept_sigma", 
            "b_mu", "b_sigma", 
            "sigma", "mu")
)


# function to compute r2
r2_compute <- function(model, log = TRUE) {
  # model = m1
  mu_hat <- as.matrix(model, "mu")
  sigma_hat <- as.numeric(as.matrix(model, "sigma"))
  
  # dim(mu_hat)
  
  if(log) {
    mu_var <- apply(mu_hat, 1, var)
    r2log <- mu_var / (mu_var + sigma_hat ^ 2)
    return(r2log)
  } else {
    mu_hat_exp <- exp(mu_hat + 0.5 * sigma_hat ^ 2)
    # dim(mu_hat_exp)
    mu_var_exp <- apply(mu_hat_exp, 1, var)
    sigma2_exp <- (exp(sigma_hat ^ 2) - 1) * exp(2 * mu_hat + sigma_hat ^ 2)
    # dim(sigma2_exp)
    r2exp <- mu_var_exp / (mu_var_exp + rowMeans(sigma2_exp))
    return(r2exp)
  }
}

r2m1 <- r2_compute(m1)
r2m2 <- r2_compute(m2)
r2m3 <- r2_compute(m3)

plot(density(r2m1, from = 0, to = 1), ylim = c(0, 10)) # re bajo
abline(v = mean(r2m1), lty = 2, col = 1)

lines(density(r2m2, from = 0, to = 1), col = 2) # re bajo
abline(v = mean(r2m2), lty = 2, col = 2)

lines(density(r2m3, from = 0, to = 1), col = 3) # re bajo
abline(v = mean(r2m3), lty = 2, col = 3)

# dan re igual, y quizás exp y focal son mejores. por poco.
 


# En exp tienden a 0
r2m1 <- r2_compute(m1, F)
r2m2 <- r2_compute(m2, F)
r2m3 <- r2_compute(m3, F)

plot(density(r2m1, from = 0, to = 1), ylim = c(0, 300), xlim = c(0, 0.05)) # re bajo
abline(v = mean(r2m1), lty = 2, col = 1)

lines(density(r2m2, from = 0, to = 1), col = 2) # re bajo
abline(v = mean(r2m2), lty = 2, col = 2)

lines(density(r2m3, from = 0, to = 1), col = 3) # re bajo
abline(v = mean(r2m3), lty = 2, col = 3)




# Conclu ------------------------------------------------------------------

# Usar máximo del fortnight focal. Todo es muy ruidoso.
# Limpiar código.
