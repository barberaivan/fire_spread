###
### Had troubles to fit temp and pp models, so for simplicity we just kept
### the FWI models.
###


# Exploring which climatic variables explain better fire size, to be used as
# predictor of the steps parameter in the spread model.
# Temperature, precipitation and FWI, will be considered, based on ERA5 data.

# We need to define a temporal window
# relevant for fire spread, e.g., one week, one month, or 3 months
# before the ignition. However, a more continuos approach is to consider that the
# climatic effect on fire spread decreases as a function of the temporal lag.
# Hence, we can fit cumulative effects models.

# In fires with dates checked by records we define the date as the beginning
# date, as the ending is usually many days after the main spread events.
# In other fires which lasted too long we choose mid-dates in which most of
# the fire spread.

# The cumulative climate is based on that date, taken as date 0. Hence, these
# models compute a weighted average, with a weight decreasing with lag.
# For daily values we'll consider L = 120 days as the largest lag.

# log(mu[d]) = alpha + b * sum_{l = 0}^(-L) w[l] * climate[d-l]
# w[l] = wu[l] / sum(wu)
# wu[l] = exp(- (l / lengthscale) ^ 2)

# where d is the day, and l is the lag index in days (negative).
# alpha, b and lengthscale are parameters to estimate. The larger the scale,
# the more lasting is the climatic effect. Bear in mind, for interpretation, that
# FWI is already a cumulative variable. However, it presents considerable
# daily variation.

# However, defining periods independent of fire dates could be better for
# further landscape simulations.
# We could take the mean, median or max by week, fortnight, and month, and even
# consider cumulative effects of those values.

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(terra)
library(DHARMa)
library(lubridate)
library(rstan)
library(mgcv)
library(viridis)

# Functions ---------------------------------------------------------------

normalize <- function(x) x / sum(x)

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

r2_bayes <- function(mu, sigma) {
  var(mu) / (var(mu) + mean(sigma ^ 2))
}


# Constants ---------------------------------------------------------------

var_names <- c("fwi", "temp", "pp")
nvar <- length(var_names)

# periods for temporal aggregation
per_names <- c("day", "week", "fortnight")
nper <- length(per_names)

# model types
mod_names <- c("raw", "cumulative")
nmod <- length(mod_names)


# Data --------------------------------------------------------------------

fires_vec <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp")
fires <- as.data.frame(fires_vec)
# remove fires with uncertain year
uncertain <- which(fires$obs == "uncertain_year")
fires <- fires[-uncertain, ]

fwi_data <- read.csv("data/climatic_data_by_fire_fwi.csv")
temp_data <- read.csv("data/climatic_data_by_fire_temp.csv")
pp_data <- read.csv("data/climatic_data_by_fire_pp.csv")

fwi_data <- fwi_data[order(fwi_data$fire_id), ]
temp_data <- temp_data[order(temp_data$fire_id), ]
pp_data <- pp_data[order(pp_data$fire_id), ]
fires <- fires[order(fires$fire_id), ]

all(fires$fire_id == fwi_data$fire_id) # OK

# Turn dates into dates
dd <- c("date_start", "date_end", "date_l")
for(d in dd) {
  fwi_data[, d] <- as.Date(fwi_data[, d])
  temp_data[, d] <- as.Date(temp_data[, d])
  pp_data[, d] <- as.Date(pp_data[, d])
}

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

# date uncertainty
fires$datediff <- fires$date_u - fires$date_l
rows_fit <- fires$datediff < 20

# fires without much date uncertainty
fires_sub <- fires[rows_fit, ]

# compare middate modis with thomas dates
fires$middate_sat <- fires$date_l + floor((fires$date_u - fires$date_l) / 2)

delta_sat_th <- fires$date - fires$middate_sat
# hist(as.numeric(delta_sat_th[delta_sat_th < 100]), breaks = 30)
# hist(as.numeric(delta_sat_th[grep("modis", fires$datesrc)]), breaks = 30)
# plot(table(abs(as.numeric(delta_sat_th[grep("modis", fires$datesrc)]))))
# plot(table(abs(as.numeric(delta_sat_th[abs(delta_sat_th) < 100]))))


# Subset also climate datasets
fwi_sub <- fwi_data[rows_fit, ]
temp_sub <- temp_data[rows_fit, ]
pp_sub <- pp_data[rows_fit, ]


# For each variable, create matrices of lagged data, with different temporal
# aggregates

fires_sub$fwi_mat <- fwi_sub[, grep("fwi", colnames(fwi_sub))] %>% as.matrix
fires_sub$temp_mat <- temp_sub[, grep("temp", colnames(temp_sub))] %>% as.matrix
fires_sub$pp_mat <- pp_sub[, grep("pp", colnames(pp_sub))] %>% as.matrix

fires_sub$date_start <- fwi_sub$date_start
fires_sub$date_end <- fwi_sub$date_end

# a list with fires, containing a list with periods
clim_list <- lapply(1:length(fires_sub$fire_id), function(f) {
  # f = 1

  # get climate values
  dtemp <- fires_sub[f, , drop = F]

  clim_vals <- data.frame(
    fwi = as.numeric(dtemp$fwi_mat),
    temp = as.numeric(dtemp$temp_mat),
    pp = as.numeric(dtemp$pp_mat)
  )

  # subset focal dates
  date_seq_av <- seq(dtemp$date_start, dtemp$date_end, 1)
  middate <- fires_sub$date[f]

  if(length(date_seq_av) != nrow(clim_vals)) stop("Date sequence does not match FWI values.")

  dd <- data.frame(day = date_seq_av,
                   week = week(date_seq_av),
                   month = month(date_seq_av),
                   year = year(date_seq_av))
  dd <- cbind(dd, clim_vals)

  dd <- dd[order(dd$day), ]

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

  # create fortnight, centered at the new year
  odd_weeks <- (dd$week %% 2) == 0
  dd$fortnight <- dd$week
  dd$fortnight[odd_weeks & dd$week > 0] <- dd$week[odd_weeks & dd$week > 0] - 1
  dd$fortnight[!odd_weeks & dd$week < 0] <- dd$week[!odd_weeks & dd$week < 0] + 1

  dd_mid <- dd[dd$day == middate, , drop = F]

  # Subset data and aggregate by periods

  per_list <- vector("list", nper)
  names(per_list) <- per_names

  for(period in per_names) {
    # period = "day"
    # subset dd using only the focal and previous periods
    dd_sub <- dd[dd[, period] <= dd_mid[, period], ]
    dd_sub$period <- dd_sub[, period]

    # aggregate fwi by period
    if(period != "day") {
      dd_agg <- aggregate(cbind(fwi, temp, pp) ~ period, dd_sub, mean)
      dd_agg$period <- -(nrow(dd_agg) - 1):0
      dd_agg <- dd_agg[order(dd_agg$period, decreasing = T), ]

      per_list[[period]] <- dd_agg
    } else {
      dd_agg <- dd_sub[, c("period", "fwi", "temp", "pp")]
      dd_agg$period <- -(nrow(dd_agg) - 1):0
      dd_agg <- dd_agg[dd_agg$period >= -120, ]
      dd_agg <- dd_agg[order(dd_agg$period, decreasing = T), ]

      per_list[[period]] <- dd_agg
    }
  }

  return(per_list)
})

names(clim_list) <- fires_sub$fire_id

# get min number of dates across fire, by period
min_periods <- numeric(nper)
names(min_periods) <- per_names
for(p in 1:nper) {
  # p = 1
  min_periods[p] <- lapply(clim_list, function(x) nrow(x[[p]])) %>% unlist %>% min
}

# make list of arrays, one by period
clim_arrs <- vector("list", nper)
names(clim_arrs) <- c(per_names)

for(i in 1:nper) {
  # i = 1
  aa <- array(NA, dim = c(nrow(fires_sub),
                          min_periods[i],
                          nvar),
              dimnames = list(
                fire_id = fires_sub$fire_id,
                time = 0 : (-(min_periods[i] - 1)),
                var = var_names
              ))
  for(f in 1:nrow(fires_sub)) {
    xx <- clim_list[[f]][[per_names[i]]]
    xx <- xx[1:min_periods[i], ]
    aa[f, , ] <- t(as.matrix(xx[, var_names]))
  }

  clim_arrs[[i]] <- aa
}


# Fit bayesian cumulative models ------------------------------------------
# to get the point estimates of parameters. Then, the predictors will be
# computed in a cumulative way, based on the point estimates, and many
# lms will be fitted.

lengthscales <- matrix(NA, nvar, nper)
colnames(lengthscales) <- per_names
rownames(lengthscales) <- var_names

models <- vector("list", nper)
names(models) <- per_names
for(i in 1:nper) {
  models[[i]] <- vector("list", nvar)
  names(models[[i]]) <- var_names
}

fires_sub$area_ha_log <- log(fires_sub$area_ha)

model_expquad <- stan_model("temporal scale for climatic predictors/cumulative model.stan")

for(p in per_names) {
  for(v in var_names) {
    # p = "day"; v = "temp"

    clim_mat <- clim_arrs[[p]][, , v]
    # standardize
    clim_mat <- (clim_mat - mean(clim_mat[, 1])) / sd(clim_mat[, 1])

    times <- ncol(clim_mat)
    y <- as.numeric(scale(fires_sub$area_ha_log))

    stan_data <- list(
      y = y,
      clim_mat = clim_mat,
      N = nrow(clim_mat),
      N_times = times,

      prior_sd_alpha = 5,
      prior_mean_alpha = 0,
      prior_sd_beta = 5,
      prior_sd_sigma = 5,
      prior_sd_ls = times * 0.5
    )

    smodel <- sampling(
      model_expquad, data = stan_data, seed = 1235425,
      chains = 4, cores = 4
    )

    print(paste(p, v))
    (lengthscales[v, p] <- mean(as.matrix(smodel, "ls")))
    models[[p]][[v]] <- smodel
  }
}

saveRDS(models, "files/climatic_cumulative_models.rds")

# Function to fit cumulative-aggregate models -----------------------------

# Compile stan models
model_exp <- stan_model("temporal scale for climatic predictors/FWI model cumulative exp.stan")

cum_agg_fit <- function(period = c("day", "week", "fortnight", "month"),
                        light = TRUE) {

  ## test:
  # period <- "day"

  # get data for the previous and focal day
  fwi_list <- lapply(fires_sub$fire_id, function(f) {
    # f = fires_sub$fire_id[1]

    # get fwi values
    row <- which(fires$fire_id == f)
    dtemp <- fwi_data[fwi_data$fire_id == f, , drop = F]
    fwi_vals <- dtemp[1, grep("fwi", names(fwi_data))] %>% as.numeric
    # subset focal dates
    date_seq_av <- seq(dtemp$date_start, dtemp$date_end, 1)
    middate <- fires$date[row]

    if(length(date_seq_av) != length(fwi_vals)) stop("Date sequence does not match FWI values.")

    dd <- data.frame(fwi = fwi_vals,
                     day = date_seq_av,
                     week = week(date_seq_av),
                     month = month(date_seq_av),
                     year = year(date_seq_av))

    dd <- dd[order(dd$day), ]

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

    # create fortnight, centered at the new year
    odd_weeks <- (dd$week %% 2) == 0
    dd$fortnight <- dd$week
    dd$fortnight[odd_weeks & dd$week > 0] <- dd$week[odd_weeks & dd$week > 0] - 1
    dd$fortnight[!odd_weeks & dd$week < 0] <- dd$week[!odd_weeks & dd$week < 0] + 1

    dd_mid <- dd[dd$day == middate, , drop = F]

    # subset dd using only the focal and previous periods
    dd_sub <- dd[dd[, period] <= dd_mid[, period], ]
    dd_sub$period <- dd_sub[, period]

    # aggregate fwi by period
    if(period != "day") {
      dd_agg <- aggregate(fwi ~ period, dd_sub, mean)
      dd_agg$period <- -(nrow(dd_agg) - 1):0
      return(dd_agg)
    } else {
      dd_agg <- dd_sub[, c("period", "fwi")]
      dd_agg$period <- -(nrow(dd_agg) - 1):0
      dd_agg <- dd_agg[dd_agg$period >= -120, ]
      return(dd_agg)
    }

})
  min_periods <- lapply(fwi_list, nrow) %>% unlist %>% min

  fwi_periods <- do.call("rbind", lapply(fwi_list, function(x) {
    # x <- fwi_fort_list[[1]]
    xx <- x[order(x$period, decreasing = TRUE), ]
    fw <- xx$fwi[1:min_periods]
    mm <- matrix(fw, 1)
    colnames(mm) <- xx$period[1:min_periods]
    return(mm)
  }))

  # parameters to scale temporal correlation
  times <- min_periods
  time_scale_expquad <- times / 4 # the new time distance is t / time_scale_
  time_scale_exp <- times / 10
  prior_ls_sd <- times * 0.75
  timeseq <- 0:(ncol(fwi_periods)-1)

  # reference values to center predictors
  ls_ref <- time_scale_expquad;
  rho_ref <- 0.5

  # ## Prior check
  # # expquad model
  # curve(exp(-((x) / ls_ref) ^ 2), from = 0, to = times,
  #       ylim = c(0, 1))
  # for(i in 1:500) {
  #   ls <- abs(rnorm(1, 0, sd = prior_ls_sd))
  #   curve(exp(-(x / ls) ^ 2), add = TRUE, col = rgb(0, 0, 0, 0.05))
  # }
  # # exp model
  # curve(dlogitnorm(x, 0, 1.5), ylim = c(0, 4), n = 300)
  # curve(rho_ref ^ (x / time_scale_exp), from = 0, to = times, ylim = c(0, 1))
  # for(i in 1:500) {
  #   # rho <- runif(1, 0, 1)
  #   # rho <- rbeta(1, 0.8, 0.8)
  #   rho <- plogis(rnorm(1, 0, 1.5))
  #   curve(rho ^ (x / time_scale_exp), add = TRUE, col = rgb(0, 0, 0, 0.05))
  # }

  # center predictors based on reference values
  fwi_center_exp <- mean(fwi_periods %*% normalize(rho_ref ^ (timeseq / time_scale_exp)))
  fwi_center_expquad <- mean(fwi_periods %*% normalize(exp(-(timeseq / ls_ref) ^ 2)))

  fwi_scale_exp <- sd(fwi_periods %*% normalize(rho_ref ^ (timeseq / time_scale_exp)))
  fwi_scale_expquad <- sd(fwi_periods %*% normalize(exp(-(timeseq / ls_ref) ^ 2)))

  # the same in log, for sigma
  fwi_scale_exp_log <- sd(log(fwi_periods %*% normalize(rho_ref ^ (timeseq / time_scale_exp))))
  fwi_scale_expquad_log <- sd(log(fwi_periods %*% normalize(exp(-(timeseq / ls_ref) ^ 2))))

  # data for stan
  sdata_exp <- list(
    y = log(fires_sub$area_ha),
    fwi_mat = fwi_periods,
    N_times = ncol(fwi_periods),
    N = nrow(fwi_periods),

    fwi_center = fwi_center_exp,
    fwi_scale = fwi_scale_exp,
    fwi_scale_log = fwi_scale_exp_log, # not center because cant log a negative

    prior_mean_int_mu = log(1000),
    prior_mean_int_sigma = log(2), # sd(log(fires_sub$area_ha))

    prior_sd_int_mu = 20,
    prior_sd_int_sigma = 10, # sd(log(fires_sub$area_ha))

    prior_sd_b_mu = 10,
    prior_sd_b_sigma = 2,

    prior_sd_ls = prior_ls_sd,
    prior_sd_rho = 1.5, # logitnormal prior

    time_scale_expquad = time_scale_expquad,
    time_scale_exp = time_scale_exp
  )

  sdata_expquad <- sdata_exp
  sdata_expquad$fwi_center = fwi_center_expquad
  sdata_expquad$fwi_scale = fwi_scale_expquad
  sdata_expquad$fwi_scale_log = fwi_scale_expquad_log

  # Fit both models
  m_expquad <- sampling(
    model_expquad, sdata_expquad, seed = 123,
    cores = 4, chains = 4, iter = 2000,
    pars =  c("alpha_mu", "alpha_sigma",
              "intercept_mu", "intercept_sigma",
              "b_mu", "b_sigma", "ls",
              "sigma", "mu")
  )

  m_exp <- sampling(
    model_exp, sdata_exp, seed = 123,
    cores = 4, chains = 4, iter = 2000,
    pars =  c("alpha_mu", "alpha_sigma",
              "intercept_mu", "intercept_sigma",
              "b_mu", "b_sigma", "rho",
              "sigma", "mu")
  )

  # extract correlation parameters to compute the cumulative FWI
  map_iter_expquad <- which.max(as.numeric(as.matrix(m_expquad, "lp__")))
  ls_hat <- as.matrix(m_expquad, "ls")[map_iter_expquad]

  map_iter_exp <- which.max(as.numeric(as.matrix(m_exp, "lp__")))
  rho_hat <- as.matrix(m_exp, "rho")[map_iter_exp]

  # cumulative fwi
  w_expquad <- normalize(exp(-(timeseq / ls_hat) ^ 2))
  fwi_expquad <- fwi_periods %*% w_expquad

  w_exp <- normalize(rho_hat ^ (timeseq / time_scale_exp))
  fwi_exp <- fwi_periods %*% w_exp

  # return the focal value and the cumulative values
  dat <- data.frame(
    "focal" = fwi_periods[, 1],
    "expquad" = fwi_expquad,
    "exp" = fwi_exp
  )

  names(dat) <- paste(rep(period, 3), names(dat), sep = "_")

  if(light) {
    return(dat)
  } else {
    return(
      list(fwi_values = dat,
           m_exp = m_exp,
           m_expquad = m_expquad,
           sdata_exp = sdata_exp,
           sdata_expquad = sdata_expquad,
           ls_hat = ls_hat,
           rho_hat = rho_hat)
    )
  }
}


# Cumulative FWI from cumulative bayesian models --------------------------

day_fit <- cum_agg_fit("day", light = FALSE) # ok
week_fit <- cum_agg_fit("week") # ok
fortnight_fit <- cum_agg_fit("fortnight", light = FALSE) # ok
month_fit <- cum_agg_fit("month") # okok

fwi_cum <- cbind(day_fit, week_fit, fortnight_fit, month_fit)


# Fitting all models with mgcv gaulss ------------------------------------

models_r2 <- data.frame(
  predictor = colnames(fwi_cum),
  ll = NA,
  r2_bayes = NA, # this is deviance explained according to gam
  r2_freq = NA
)

for(c in 1:ncol(fwi_cum)) {
  # c = 1
  vv <- names(fwi_cum)[c]
  y <- log(fires_sub$area_ha)
  x <- fwi_cum[, c]
  xlog <- log(x)

  m <- gam(list(y ~ x, ~ xlog), family = gaulss())
  ff <- fitted(m)
  sm <- summary(m)

  models_r2[c, "ll"] <- logLik(m) %>% as.numeric
  models_r2[c, "r2_bayes"] <- r2_bayes(ff[, 1], 1 / ff[, 2])
  models_r2[c, "r2_freq"] <- sm$dev.expl
  # In this way I find the same results as before, with very high
  # r2 for the max model.
}

models_r2

# plots
fwi_plot <- fwi_cum
fwi_plot$area_ha <- fires_sub$area_ha
fpl <- pivot_longer(fwi_plot, 1:ncol(fwi_cum),
                    names_to = "predictor", values_to = "fwi")
fpl <- separate(fpl, "predictor", into = c("period", "model"), sep = "_",
                remove = F)
fpl$period <- factor(fpl$period, levels = c("day", "week", "fortnight", "month"))
fpl$model <- factor(fpl$model, levels = c("focal", "exp", "expquad"))

# data for plotting r2
r2plot <- models_r2
r2plot$fwi <- 58
r2plot$area_ha <- 30

r2plot$r2_freq <- paste(round(r2plot$r2_freq * 100, 2), "%")
r2plot$r2_bayes <- paste(round(r2plot$r2_bayes * 100, 2), "%")

r2plot <- separate(r2plot, "predictor", into = c("period", "model"), sep = "_",
                remove = F)
r2plot$period <- factor(r2plot$period, levels = c("day", "week", "fortnight", "month"))
r2plot$model <- factor(r2plot$model, levels = c("focal", "exp", "expquad"))

# plot
ggplot(fpl, aes(fwi, area_ha)) +
  geom_point(size = 1.5, alpha = 0.25) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.25,
              fill = "purple") +
  geom_text(data = r2plot, mapping = aes(label = r2_freq),
            size = 3, color = "blue") +
  geom_text(data = r2plot, mapping = aes(y = area_ha + 35, label = r2_bayes),
            size = 3, color = "darkred") +
  # facet_wrap(vars(predictor), scales = "free_x") +
  facet_grid(rows = vars(period), cols = vars(model)) +
  scale_y_continuous(trans = "log10") +
  theme(panel.grid.minor = element_blank()) +
  ylab("Fire size (ha)") +
  xlab("Fire weather index")
ggsave("temporal scale for climatic predictors/comparison.png",
       width = 18, height = 20, units = "cm")

GGally::ggpairs(fortnight_fit)
# usar fortnight expquad


# fortnight models  -------------------------------------------------------

fortnight_mod <- cum_agg_fit("fortnight", light = F) # ok

ls_samples <- as.matrix(fortnight_mod$m_expquad, "ls") %>% as.numeric()
curve(exp(- (x / fortnight_mod$ls_hat) ^ 2),
      to = fortnight_mod$sdata_expquad$N_times)
for(i in 1:length(ls_samples)) {
  curve(exp(- (x / ls_samples[i]) ^ 2), add = TRUE, col = rgb(0, 0, 0, 0.05))
}
curve(exp(- (x / fortnight_mod$ls_hat) ^ 2), add = TRUE, col = "blue", lwd = 2)

# prior-posterior
plot(density(ls_samples, from = 0), main = "Expquad prior-posterior")
curve(dnorm(x, 0, fortnight_mod$sdata_expquad$prior_sd_ls) * 2, add = TRUE, col = 2)
abline(v = fortnight_mod$ls_hat, lty = 2)
# usa hasta 3 meses atrás

# fortnight distance with 0.5 importance
dis <- 0.5
(fort_dis <- fortnight_mod$ls_hat * sqrt(log(1 / 0.5)))
dw <- data.frame(
  fortnight = 0:10,
  weight = exp(- (0:10 / fortnight_mod$ls_hat) ^ 2) %>% normalize
)
dw$weight_cum <- cumsum(dw$weight)
dw
# fortnight       weight weight_cum
#         0 3.561154e-01  0.3561154
#         1 3.072998e-01  0.6634153
#         2 1.974590e-01  0.8608743
#         3 9.447883e-02  0.9553531
#         4 3.366163e-02  0.9890147
#         5 8.930564e-03  0.9979453
#         6 1.764272e-03  0.9997096
#         7 2.595347e-04  0.9999691
#         8 2.842944e-05  0.9999975
#         9 2.318912e-06  0.9999999
#        10 1.408455e-07  1.0000000
fortnight_mod$ls_hat # 2.604387 fortnights


# daily models (the best) --------------------------------------------------

ls_samples <- as.matrix(day_fit$m_expquad, "ls") %>% as.numeric()
curve(exp(- (x / day_fit$ls_hat) ^ 2),
      to = day_fit$sdata_expquad$N_times)
for(i in 1:length(ls_samples)) {
  curve(exp(- (x / ls_samples[i]) ^ 2), add = TRUE, col = rgb(0, 0, 0, 0.05))
}
curve(exp(- (x / day_fit$ls_hat) ^ 2), add = TRUE, col = "blue", lwd = 2)

# prior-posterior
plot(density(ls_samples, from = 0), main = "Expquad prior-posterior")
curve(dnorm(x, 0, day_fit$sdata_expquad$prior_sd_ls) * 2, add = TRUE, col = 2)
abline(v = day_fit$ls_hat, lty = 2)
# usa hasta 20 días atrás

# day distance with 0.5 importance
dis <- 0.5
(day_dis <- day_fit$ls_hat * sqrt(log(1 / 0.5))) # 8.4 días

options(scipen = 999)
dw <- data.frame(
  day = 0:120,
  weight = exp(- (0:120 / day_fit$ls_hat) ^ 2) %>% normalize
)
dw$weight_cum <- cumsum(dw$weight)
dw

# en 15 días acumula el 95 % del weight

# export data with FWI cumulative values (using all fires!) --------------

# get data for the previous and focal day
# this repeats the beginning of the cum_fit
period <- "day"
fwi_list <- lapply(fires$fire_id, function(f) {
  # f = fires_sub$fire_id[1]

  # get fwi values
  row <- which(fires$fire_id == f)
  dtemp <- fwi_data[fwi_data$fire_id == f, , drop = F]
  fwi_vals <- dtemp[1, grep("fwi", names(fwi_data))] %>% as.numeric
  # subset focal dates
  date_seq_av <- seq(dtemp$date_start, dtemp$date_end, 1)
  middate <- fires$date[row]

  if(length(date_seq_av) != length(fwi_vals)) stop("Date sequence does not match FWI values.")

  dd <- data.frame(fwi = fwi_vals,
                   day = date_seq_av,
                   week = week(date_seq_av),
                   month = month(date_seq_av),
                   year = year(date_seq_av))

  dd <- dd[order(dd$day), ]

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

  # create fortnight, centered at the new year
  odd_weeks <- (dd$week %% 2) == 0
  dd$fortnight <- dd$week
  dd$fortnight[odd_weeks & dd$week > 0] <- dd$week[odd_weeks & dd$week > 0] - 1
  dd$fortnight[!odd_weeks & dd$week < 0] <- dd$week[!odd_weeks & dd$week < 0] + 1

  dd_mid <- dd[dd$day == middate, , drop = F]

  # subset dd using only the focal and previous periods
  dd_sub <- dd[dd[, period] <= dd_mid[, period], ]
  dd_sub$period <- dd_sub[, period]

  # aggregate fwi by period
  if(period != "day") {
    dd_agg <- aggregate(fwi ~ period, dd_sub, mean)
    dd_agg$period <- -(nrow(dd_agg) - 1):0
    return(dd_agg)
  } else {
    dd_agg <- dd_sub[, c("period", "fwi")]
    dd_agg$period <- -(nrow(dd_agg) - 1):0
    dd_agg <- dd_agg[dd_agg$period >= -120, ]
    return(dd_agg)
  }

})
min_periods <- lapply(fwi_list, nrow) %>% unlist %>% min

fwi_periods <- do.call("rbind", lapply(fwi_list, function(x) {
  # x <- fwi_fort_list[[1]]
  xx <- x[order(x$period, decreasing = TRUE), ]
  fw <- xx$fwi[1:min_periods]
  mm <- matrix(fw, 1)
  colnames(mm) <- xx$period[1:min_periods]
  return(mm)
}))

# extract correlation parameters to compute the cumulative FWI
map_iter_expquad <- which.max(as.numeric(as.matrix(day_fit$m_expquad, "lp__")))
ls_hat <- as.matrix(day_fit$m_expquad, "ls")[map_iter_expquad]

map_iter_exp <- which.max(as.numeric(as.matrix(day_fit$m_exp, "lp__")))
rho_hat <- as.matrix(day_fit$m_exp, "rho")[map_iter_exp]

# cumulative fwi
timeseq <- 0:120
w_expquad <- normalize(exp(-(timeseq / ls_hat) ^ 2))
fwi_expquad <- fwi_periods %*% w_expquad

w_exp <- normalize(rho_hat ^ (timeseq / day_fit$sdata_exp$time_scale_exp))
fwi_exp <- fwi_periods %*% w_exp

# return the focal value and the cumulative values
dat <- data.frame(
  "focal" = fwi_periods[, 1],
  "expquad" = fwi_expquad,
  "exp" = fwi_exp
)

names(dat) <- c("fwi_focal_day", "fwi_expquad_day", "fwi_exp_day")
export <- cbind(
  fires[, c("fire_id", "date", "date_l", "date_u", "date_lr", "date_ur", "obs",
            "datediff")],
  dat
)

write.csv(export, "data/climatic_data_by_fire_fwi-day-cumulative.csv",
          row.names = F)



# export data with FWI cumulative values (using all fires! -- FORTNIGHT) ------

# get data for the previous and focal fortnight
# this repeats the beginning of the cum_fit
period <- "fortnight"
fwi_list <- lapply(fires$fire_id, function(f) {
  # f = fires_sub$fire_id[1]

  # get fwi values
  row <- which(fires$fire_id == f)
  dtemp <- fwi_data[fwi_data$fire_id == f, , drop = F]
  fwi_vals <- dtemp[1, grep("fwi", names(fwi_data))] %>% as.numeric
  # subset focal dates
  date_seq_av <- seq(dtemp$date_start, dtemp$date_end, 1)
  middate <- fires$date[row]

  if(length(date_seq_av) != length(fwi_vals)) stop("Date sequence does not match FWI values.")

  dd <- data.frame(fwi = fwi_vals,
                   day = date_seq_av,
                   week = week(date_seq_av),
                   month = month(date_seq_av),
                   year = year(date_seq_av))

  dd <- dd[order(dd$day), ]

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

  # create fortnight, centered at the new year
  odd_weeks <- (dd$week %% 2) == 0
  dd$fortnight <- dd$week
  dd$fortnight[odd_weeks & dd$week > 0] <- dd$week[odd_weeks & dd$week > 0] - 1
  dd$fortnight[!odd_weeks & dd$week < 0] <- dd$week[!odd_weeks & dd$week < 0] + 1

  dd_mid <- dd[dd$day == middate, , drop = F]

  # subset dd using only the focal and previous periods
  dd_sub <- dd[dd[, period] <= dd_mid[, period], ]
  dd_sub$period <- dd_sub[, period]

  # aggregate fwi by period
  if(period != "day") {
    dd_agg <- aggregate(fwi ~ period, dd_sub, mean)
    dd_agg$period <- -(nrow(dd_agg) - 1):0
    return(dd_agg)
  } else {
    dd_agg <- dd_sub[, c("period", "fwi")]
    dd_agg$period <- -(nrow(dd_agg) - 1):0
    dd_agg <- dd_agg[dd_agg$period >= -120, ]
    return(dd_agg)
  }

})
min_periods <- lapply(fwi_list, nrow) %>% unlist %>% min

fwi_periods <- do.call("rbind", lapply(fwi_list, function(x) {
  # x <- fwi_fort_list[[1]]
  xx <- x[order(x$period, decreasing = TRUE), ]
  fw <- xx$fwi[1:min_periods]
  mm <- matrix(fw, 1)
  colnames(mm) <- xx$period[1:min_periods]
  return(mm)
}))

# extract correlation parameters to compute the cumulative FWI
map_iter_expquad <- which.max(as.numeric(as.matrix(fortnight_fit$m_expquad, "lp__")))
ls_hat <- as.matrix(fortnight_fit$m_expquad, "ls")[map_iter_expquad]

map_iter_exp <- which.max(as.numeric(as.matrix(fortnight_fit$m_exp, "lp__")))
rho_hat <- as.matrix(fortnight_fit$m_exp, "rho")[map_iter_exp]

# cumulative fwi
timeseq <- 0:(min_periods-1)
w_expquad <- normalize(exp(-(timeseq / ls_hat) ^ 2))
fwi_expquad <- fwi_periods %*% w_expquad

w_exp <- normalize(rho_hat ^ (timeseq / fortnight_fit$sdata_exp$time_scale_exp))
fwi_exp <- fwi_periods %*% w_exp

# return the focal value and the cumulative values
dat <- data.frame(
  "focal" = fwi_periods[, 1],
  "expquad" = fwi_expquad,
  "exp" = fwi_exp
)

names(dat) <- c("fwi_focal_fortnight", "fwi_expquad_fortnight", "fwi_exp_fortnight")
export <- cbind(
  fires[, c("fire_id", "date", "date_l", "date_u", "date_lr", "date_ur", "obs",
            "datediff", "area_ha")],
  dat
)

# check
# pairs(export[, c("fwi_focal_fortnight", "fwi_expquad_fortnight", "fwi_exp_fortnight")])

write.csv(export, "data/climatic_data_by_fire_fwi-fortnight-cumulative.csv",
          row.names = F)
