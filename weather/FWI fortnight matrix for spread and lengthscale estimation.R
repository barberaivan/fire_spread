# This code inherits from <FWI temporal scale.R>, where several temporal
# aggreation schemes were tested, choosing fortnights. In that case, the
# raw FWI was used, but now I use the FWI temporal anomalies. This was needed
# because the spread model made huge fires in the steppe, under extreme FWI raw
# values that were not present in the fires used for spread. Hence, we ignore
# the spatial variation in FWI.

# Here we only estimate the lenghtscale based on a regression of fire size and
# (cumulative) FWI. The FWI raster standardized and aggregated by fortnight
# was produced in <FWI standardize and aggregate by fortnight.R>.

# The FWI for each fire is interpolated at the ignition point, and at the polygon
# centroid if the ignition point is unknown. Getting ignition points data will
# take a considerable portion of the script.

# In the parent script, the FWI of each fire was obtained in GEE, creating the
# "data/climatic_data_by_fire_fwi.csv" file. Here I bring the data from the
# raster.

# The exports ending in FWIZ2 use the FWI at 24 km.
 
# Here I also export FWI matrix with lags, in case it is sensible to fit the
# lengthscale alongside spread parameters, in the second stage of fitting.
# The exports have the "matrix" word to indicate that the lagged FWI matrix
# is included.


# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(terra)
library(DHARMa)
library(lubridate)
library(rstan)
library(mgcv)
source(file.path("weather", "fortnight_functions.R"))

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

tseq <- 0:10 # previous fortnights considered (0 is the focal)
nt <- length(tseq)
lagnames <- stringr::str_pad(tseq, width = 2, side = "left", pad = "0")

# Data --------------------------------------------------------------------

# fwi_rast <- rast(file.path("data", "fwi_daily_1998-2022", "52km",
#                            "fwi_fortnights_19980101_20230630_standardized.tif"))
fwi_rast <- rast(file.path("data", "fwi_daily_1998-2022", "24km",
                           "fwi_fortnights_19970701_20230630_standardized.tif"))

fires_vec <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp")

igpoints <- vect(file.path("data", "ignition_points_checked.shp"))

# Tidy fires dates --------------------------------------------------------

# remove fires with uncertain year
uncertain <- which(fires_vec$obs == "uncertain_year")
fires_vec <- fires_vec[-uncertain, ]

# Cast dates
dd <- c("date", "date_l", "date_u")
fires_vec$date <- as.Date(fires_vec$date)
fires_vec$date_l <- as.Date(fires_vec$date_l)
fires_vec$date_u <- as.Date(fires_vec$date_u)

# Correct date_l and date_u for Landsat-based images, so they are included in the
# focal period. For landsat, date_l and date_u are the dates between which the
# fire occurred, but not included. For modis, they are inclusive. Turn all into
# inclusive. Do this only for date_diff larger than 2
ll <- (fires_vec$datesrc %in% c("landsat", "Mermoz")) &
  (fires_vec$date_u - fires_vec$date_l > 2)
fires_vec$date_l[ll] <- fires_vec$date_l[ll] + 1
fires_vec$date_u[ll] <- fires_vec$date_u[ll] - 1

# compare middate modis with thomas dates
fires_vec$middate_sat <- fires_vec$date_l + floor((fires_vec$date_u - fires_vec$date_l) / 2)
delta_sat_th <- fires_vec$date - fires_vec$middate_sat
# hist(as.numeric(delta_sat_th[delta_sat_th < 100]), breaks = 30)
# hist(as.numeric(delta_sat_th[grep("modis", fires_vec$datesrc)]), breaks = 30)
# plot(table(abs(as.numeric(delta_sat_th[grep("modis", fires_vec$datesrc)]))))
# plot(table(abs(as.numeric(delta_sat_th[abs(delta_sat_th) < 100]))))

# Get fortnight
fires_vec$fort <- date2fort(fires_vec$date)

# column to match names with ignition data
fires_vec$fire_id_match <- fires_vec$fire_id

# Get coordinates for each fire -------------------------------------------

# A few names differ between the polygons and the ignition points.
igpoints$fire_id_ig <- igpoints$Name
igpoints$fire_id_match <- igpoints$Name

# rename repeated names (fires that were divided) to match
igpoints$fire_id_match[grep("2015_47", igpoints$fire_id_ig)] <- "2015_47"
igpoints$fire_id_match[grep("2011_19", igpoints$fire_id_ig)] <- "2011_19"

# get date and fortnight
igdata <- as.data.frame(igpoints)
igdata <- left_join(
  igdata,
  as.data.frame(fires_vec)[, c("fire_id_match", "date", "fort")],
  by = "fire_id_match"
)
igpoints$date <- igdata$date
igpoints$fort <- igdata$fort

## write ignition points with date and fire_id_match
# writeVector(igpoints, file.path("data", "ignition_points_checked_with_date.shp"))
igpoints <- vect(file.path("data", "ignition_points_checked_with_date.shp"))
names(igpoints)[names(igpoints) == "fire_id_ma"] <- "fire_id_match"

# get FWI for the fires with ignition
ig_fwi <- terra::extract(fwi_rast, igpoints, method = "bilinear",
                         raw = T)[, -1] # remove point id
ig_fwi_mat <- matrix(NA, length(igpoints), nt)
colnames(ig_fwi_mat) <- lagnames

for(i in 1:length(igpoints)) {
  focfort <- igpoints$fort[i]
  if(is.na(focfort)) next
  fortseq <- seq(focfort, focfort-(nt-1), by = -1) |> as.character()
  ig_fwi_mat[i, ] <- ig_fwi[i, fortseq]
}

## get FWI for the fire polygons
poly_fwi <- terra::extract(
  fwi_rast,
  centroids(fires_vec),
  method = "bilinear", raw = T
)[, -1] # remove point id

poly_fwi_mat <- matrix(NA, length(fires_vec), nt)
colnames(poly_fwi_mat) <- tseq * (-1)

for(i in 1:length(fires_vec)) {
  print(i)
  # i = 1
  # If fire has ignition point, use its FWI
  fid <- fires_vec$fire_id_match[i]
  if(fid %in% igpoints$fire_id_match) { # fire_id_match was shortened
    rows <- which(igpoints$fire_id_match == fid)
    fwi_vals <- colMeans(ig_fwi_mat[rows, , drop = F])
    poly_fwi_mat[i, ] <- fwi_vals
  } else {
    focfort <- fires_vec$fort[i]
    if(is.na(focfort)) next
    fortseq <- seq(focfort, focfort-(nt-1), by = -1) |> as.character()
    poly_fwi_mat[i, ] <- poly_fwi[i, fortseq]
  }
}

# Export lagged FWI matrix ------------------------------------------------
# to fit the lengthscale in the MCMC for spread.

# full fire database
fires_data <- as.data.frame(fires_vec)
fwi_mat_full <- poly_fwi_mat
colnames(fwi_mat_full) <- paste("fwi", lagnames, sep = "_")
fires_data <- cbind(fires_data, fwi_mat_full)

write.csv(
  fires_data,
  file.path("data", "climatic_data_by_fire_fwi-fortnight-matrix_FWIZ.csv"),
  row.names = F
)

# ignition points
ig_fwi_df <- as.data.frame(ig_fwi_mat) 
colnames(ig_fwi_df) <- colnames(fwi_mat_full)
igpoints_export <- cbind(igpoints, ig_fwi_df)
names(igpoints_export)

writeVector(
  igpoints_export,
  file.path("data", "ignition_points_checked_with_date-fort-matrix-fwiz.shp"),
  overwrite = T
)

# Fit fire size model -----------------------------------------------------

# Compile stan models
model_expquad <- stan_model("weather/FWI model cumulative expquad_simpler.stan")

## Prior check lengthscale
prior_ls_sd <- 10 * 0.75 # 10 before is the furthest
# ls <- abs(rnorm(1, 0, sd = prior_ls_sd))
# curve(exp(-(x / ls) ^ 2), from = 0, to = nt,
#       ylim = c(0, 1))
# for(i in 1:500) {
#   ls <- abs(rnorm(1, 0, sd = prior_ls_sd))
#   curve(exp(-(x / ls) ^ 2), add = TRUE, col = rgb(0, 0, 0, 0.05))
# }

# data for stan
sdata <- list(
  y = log(fires_vec$area_ha),
  fwi_mat = poly_fwi_mat,
  N_times = nt,
  N = nrow(poly_fwi_mat),

  prior_mean_int_mu = mean(log(fires_vec$area_ha)),
  prior_mean_int_sigma = sd(log(fires_vec$area_ha)),

  prior_mean_a = mean(log(fires_vec$area_ha)),
  prior_sd_a = sd(log(fires_vec$area_ha)) * 5,
  prior_sd_b = 10,

  prior_mean_c = sd(log(fires_vec$area_ha)) |> log(),
  prior_sd_c = log(sd(log(fires_vec$area_ha))) * 10,
  prior_sd_d = 5,

  prior_sd_ls = prior_ls_sd
)

# sample
m_expquad <- sampling(
  model_expquad, sdata, seed = 123, cores = 4, chains = 4, iter = 2000,
  pars =  c("a", "b", "c", "d", "ls")
)
pairs(m_expquad)

# extract correlation parameters to compute the cumulative FWI
map_iter <- which.max(as.numeric(as.matrix(m_expquad, "lp__")))
ls_map <- as.matrix(m_expquad, "ls")[map_iter]
ls_mean <- as.matrix(m_expquad, "ls") |> mean()

# Use mean or map?
lsamples <- as.matrix(m_expquad, "ls") |> as.numeric()

curve(exp(-(x / ls_map) ^ 2), from = 0, to = nt,
      ylim = c(0, 1))
for(i in sample(1:length(lsamples), 200, replace = F)) {
  ls <- lsamples[i]
  curve(exp(-(x / ls) ^ 2), add = TRUE, col = rgb(0, 0, 0, 0.05))
}
curve(exp(-(x / ls_map) ^ 2), add = T, lwd = 2, col = "blue")
curve(exp(-(x / ls_mean) ^ 2), add = T, lwd = 2, col = "red")
# use map

# cumulative fwi
w_expquad <- normalize(exp(-(tseq / ls_map) ^ 2))
fwi_expquad <- poly_fwi_mat %*% w_expquad

data.frame(
  "fortnight" = tseq * (-1),
  "weight" = w_expquad,
  "weight_cum" = cumsum(w_expquad)
)
## Based on 52km FWI:

# fortnight          weight weight_cum
#         0 0.3500754639842  0.3500755
#        -1 0.3039057445114  0.6539812
#        -2 0.1988248291702  0.8528060
#        -3 0.0980295026115  0.9508355
#        -4 0.0364248232163  0.9872604
#        -5 0.0101998216669  0.9974602
#        -6 0.0021524951278  0.9996127
#        -7 0.0003423310152  0.9999550
#        -8 0.0000410303093  0.9999960
#        -9 0.0000037061057  0.9999997
#        -10 0.0000002522814 1.0000000
ls_map # 2.659056 fortnights (it decreased when variance was modelled)

## Based on 24km FWI:

# fortnight       weight weight_cum
# 1          0 3.333769e-01  0.3333769
# 2         -1 2.939973e-01  0.6273742
# 3         -2 2.016354e-01  0.8290096
# 4         -3 1.075489e-01  0.9365586
# 5         -4 4.461294e-02  0.9811715
# 6         -5 1.439234e-02  0.9955638
# 7         -6 3.610917e-03  0.9991747
# 8         -7 7.045621e-04  0.9998793
# 9         -8 1.069145e-04  0.9999862
# 10        -9 1.261740e-05  0.9999988
# 11       -10 1.158027e-06  1.0000000
ls_map # 2.820506 fortnights 

# Export data with FWI cumulative values (using all fires!) ---------------

igpoints$fwi_fort_expquad <- as.numeric(ig_fwi_mat %*% w_expquad)
igpoints$fwi_fort_focal <- as.numeric(ig_fwi_mat[, 1])

fires_data <- as.data.frame(fires_vec)
fires_data$fwi_fort_expquad <- as.numeric(poly_fwi_mat %*% w_expquad)
fires_data$fwi_fort_focal <- as.numeric(poly_fwi_mat[, 1])

write.csv(
  fires_data,
  file.path("data", "climatic_data_by_fire_fwi-fortnight-cumulative_FWIZ2.csv"),
  row.names = F
)

writeVector(
  igpoints,
  file.path("data", "ignition_points_checked_with_date-fort-fwiz2.shp")
)



# Compare FWIZ based on the 24 and 52 km ----------------------------------

f1 <- read.csv(
  file.path("data", "climatic_data_by_fire_fwi-fortnight-cumulative_FWIZ.csv")
)

f2 <- read.csv(
  file.path("data", "climatic_data_by_fire_fwi-fortnight-cumulative_FWIZ2.csv")
)

plot(f1$fwi_fort_expquad ~ f2$fwi_fort_expquad); abline(0, 1)
plot(f1$fwi_fort_focal ~ f2$fwi_fort_focal); abline(0, 1)