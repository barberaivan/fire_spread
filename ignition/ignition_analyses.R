# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(deeptime)
library(rstan)
library(lubridate)
library(terra)
library(bayesplot)

# Figure size settings ----------------------------------------------------

a4h <- 29.7
a4w <- 21.0
margins <- 2.5
fig_width_max <- a4w - margins * 2
fig_height_max <- a4h - margins * 2

# Functions ---------------------------------------------------------------

# Create a reference table with date and fortnight-date, to match observations.
# The fortnight-date (14 days) is labelled as ending-year_fortnight.
# Years run from July to June, named by the year of the ending month. Fortnights
# are centred at January first, so irregular-length ones occur at the first or
# last month, which correspond to the southern winter (less relevant for fire).

nfy <- floor(365/14)
year_low <- 1996; year_high <- 2023
years <- year_low:year_high

# Date table
dtable <- do.call("rbind", lapply(years, function(y) {
  # y = 1997
  textlow <- paste(y-1, "07", "01", sep = "-")
  textmid <- paste(y, "01", "01", sep = "-")
  texthigh <- paste(y, "06", "30", sep = "-")

  dlow <- as.Date(textlow, format = "%Y-%m-%d")
  dmid <- as.Date(textmid, format = "%Y-%m-%d")
  dhigh <- as.Date(texthigh, format = "%Y-%m-%d")

  # upper half of the year
  nd_upper <- as.numeric(dhigh - dmid) + 1
  nf_upper <- floor(nd_upper / 14)
  rem_upper <- nd_upper %% 14

  if(rem_upper >= 7) {
    fort_upper <- rep(1:(nf_upper + 1), c(rep(14, nf_upper), rem_upper))
  } else {
    lengths <- rep(14, nf_upper)
    lengths[nf_upper] <- lengths[nf_upper] + rem_upper
    fort_upper <- rep(1:nf_upper, lengths)
  }

  df_up <- data.frame(date = seq(dmid, dhigh, 1),
                      fort0 = fort_upper)

  # lower half of the year
  nd_lower <- as.numeric(dmid - dlow)
  nf_lower <- floor(nd_lower / 14)
  rem_lower <- nd_lower %% 14

  if(rem_lower >= 7) {
    fort_lower <- rep(-nf_lower:0, c(rem_lower, rep(14, nf_lower)))
  } else {
    lengths <- rep(14, nf_lower)
    lengths[1] <- lengths[1] + rem_lower
    fort_lower <- rep(-(nf_lower-1):0, lengths)
  }

  df_low <- data.frame(date = seq(dlow, dmid - 1, 1),
                       fort0 = fort_lower)

  # merge and tidy
  dd <- rbind(df_low, df_up)
  dd$fort_focal <- dd$fort0 - min(dd$fort0) + 1
  dd$year <- as.numeric(y)

  return(dd[, c("date", "fort_focal", "year")])
}))

# Get continuous fortnight identifier, from year_low to year_high.
year_ref <- dtable$year - min(dtable$year)
dtable$fort <- dtable$fort_focal + year_ref * nfy
# plot(dtable$fort)
# max(date_table$fort_cont) / length(years) # OK

# based on the table, turn date into (continuous) fortnight
date2fort <- function(d) {
  dd <- data.frame(date = d)
  dd <- left_join(dd, dtable, by = "date")
  if(anyNA(dd$fort)) {
    warning("Found dates outside the reference table, returning NA.")
  }
  return(dd$fort)
}

# Negative binomial parameterization like in Stan's neg_binomial_2:
# https://distribution-explorer.github.io/discrete/negative_binomial.html
rnegbin <- function(n = 10, mu = 10, phi = 5) {
  alpha <- phi
  beta <- alpha / mu
  lambdas <- rgamma(n, alpha, beta)
  return(rpois(n, lambdas))
}


# northing importance function
# x is the slope steepness, in degrees
northing_weight <- function(x, a = -5, b = 0.35) {
  plogis(a + b * x)
}
# northing_weighted <- cos(aspect * pi / 180) *
#                      northing_weight(slope)

normalize <- function(x) x / sum(x)

mean_ci <- function(x) {
  qq <- quantile(x, probs = c(0.025, 0.975), method = 8) %>% unname
  return(c("mean" = mean(x), "lower" = qq[1], "upper" = qq[2]))
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

# Compute prediction from ordinal model

ordinal_predict <- function(pdata, xname, bname) {
  eta <- as.matrix(scmod, bname) %*% pdata[, xname]
  cmf <- array(NA, dim = c(nrow(pdata), K-1, npost))
  for(k in 1:(K-1)) cmf[, k, ] <- plogis(ahat[, k] - eta) |> t()
  pmf <- array(NA, dim = c(nrow(pdata), K, npost))
  pmf[, 1, ] <- cmf[, 1, ]
  pmf[, K, ] <- 1 - cmf[, K-1, ]
  for(k in 2:(K-1)) pmf[, k, ] <- cmf[, k, ] - cmf[, k-1, ]

  out <- do.call("rbind", lapply(1:K, function(k) {
    summ <- apply(pmf[, k, ], 1, mean_ci) |> t() |> as.data.frame()
    summ$class <- k
    summ$class_name <- factor(class_names[k], levels = class_names)
    return(cbind(pdata, summ))
  }))

  return(out)
}

# Import vegetation class transforms --------------------------------------

dveg <- readxl::read_excel("/home/ivan/Insync/Mapa vegetación WWF - Lara et al. 1999/clases de vegetacion y equivalencias.xlsx",
                           sheet = "Sheet2")
# Make urban as forest, and get zero-indexing
dveg$cnum2[dveg$class1 == "Urban"] <- 1
dveg$class2[dveg$class1 == "Urban"] <- "Wet forest"
# urban is taken as forest so its burn probability changes markedly with NDVI.
n_veg_types <- V <- 5

# rename to use left_join
dveg$veg_focal <- dveg$cnum1
dveg$veg_num <- dveg$cnum2

# Aggregate daily FWI raster into fortnights ------------------------------

# # daily FWI stack
# fwi_day <- rast(file.path("data", "fwi_daily_1998-2022",
#                           "fwi_daily_19980101_20230630.tif"))
# fwi_day_dates <- time(fwi_day)
# fwi_day_dates <- as.Date(fwi_day_dates, format = "%Y-%m-%d")
# fwi_day_forts <- date2fort(fwi_day_dates)
# fwi_forts <- unique(fwi_day_forts)
# # aggregate raster by fortnight
# fwi_fort <- do.call("c", lapply(fwi_forts, function(f) {
#   print(f)
#   ids <- which(fwi_day_forts == f)
#   return(mean(fwi_day[[ids]]))
# })) # this takes long
# # Assign the first day to every fortnight, so the fortnight number does not
# # depend on the reference table.
# ddd <- left_join(data.frame(fort = fwi_forts), y = dtable, by = "fort")
# names(fwi_fort) <- ddd$date
# writeRaster(fwi_fort, file.path("data", "fwi_daily_1998-2022",
#                                 "fwi_fortnights_19980101_20230630.tif"),
#             overwrite = T)


# Data --------------------------------------------------------------------

# FWI
fwi_fort <- rast(file.path("data", "fwi_daily_1998-2022",
                           "fwi_fortnights_19980101_20230630.tif"))
fwi_fort <- project(fwi_fort, "EPSG:5343")

# Nahuel Huapi National Park
pnnh <- vect(file.path("data", "protected_areas", "apn_limites.shp"))
pnnh <- pnnh[pnnh$nombre == "Nahuel Huapi", ]
pnnh <- project(pnnh, "EPSG:5343")

# Ignition data is not in the fire_spread repo, it's not public
igdata_dir <- file.path("..", "ignition_data")

# Ignition points from PNNH (provided by Marcelo Bari)
ig_pnnh_data <- readxl::read_excel(file.path(igdata_dir, "Total_focos_NH_nov89-mar21.xlsx"))
ig_pnnh <- vect(ig_pnnh_data, geom = c("long", "lat"), crs = "EPSG:4326")
ig_pnnh <- project(ig_pnnh, "EPSG:5343")

# Ignitions points from Kitzberger, all caused by lightening
igl_data <- readxl::read_excel(
  file.path(igdata_dir, "base_ampliado_kitzberger_rayos.xlsx"),
  sheet = 2
)
igl_data <- igl_data[complete.cases(igl_data[, c("Lat", "Long", "Fecha")]), ]
igl0 <- vect(igl_data, geom = c("Long", "Lat"), crs = "EPSG:4326")
igl0 <- project(igl0, "EPSG:5343")


# parameters to compute flammability indices
fi_params <- readRDS(file.path("data", "NDVI_regional_data",
                               "flammability_indices.rds"))

# summary of predictors that make the flammability indices
# (used for predictions)
data_summ <- readRDS(file.path("data", "NDVI_regional_data",
                               "ndvi_elevation_summary.rds"))

# model to detrend ndvi, for the vfi
mdetrend <- readRDS(file.path("data", "NDVI_regional_data",
                              "ndvi_detrender_model.rds"))



# Remove redundant points -------------------------------------------------

igl <- igl0
igl$date <- as.Date(igl0$Fecha, format = "%Y-%m-%d")
igl$cause <- "lightning"
igl$area <- igl0$Area_est # ha
igl$id <- paste("kitz", 1:nrow(igl0), sep = "_")  # Thomas Kitzberger database
igl <- igl[, c("date", "cause", "area", "id")]

igpn <- ig_pnnh
igpn$date <- as.Date(ig_pnnh$Fecha, format = "%Y-%m-%d")
igpn$cause <- "human"
igpn$cause[ig_pnnh$Causa == "Rayo"] <- "lightning"
igpn$cause[ig_pnnh$Causa == "sin determinar"] <- "unknown"
igpn$area <- ig_pnnh$`Superficie (ha.)`
igpn$id <- paste("bari", 1:nrow(ig_pnnh), sep = "_") # Marcelo Bari database
igpn <- igpn[, c("date", "cause", "area", "id")]
# unique(igpn$cause)

# Loop over lightening ignitions in PNNH
ids_light_pn <- igpn$id[igpn$cause == "lightning"]

conflict_light <- do.call("rbind", lapply(ids_light_pn, function(id) {
  # filter by date, searching in a 5 days radius
  ## id = igpn$id[13]
  i <- which(igpn$id == id)
  dd <- igpn$date[i]
  dseq <- seq(dd-5, dd+5, 1)

  rows_kitz <- which(igl$date %in% dseq)
  ids_kitz <- igl$id[rows_kitz]
  if(length(rows_kitz) > 0) {
    crds_pn <- crds(igpn[i, ])
    crds_kitz <- crds(igl[rows_kitz, ])

    distances <- sqrt(
      (crds_pn[, "x"] - crds_kitz[, "x"]) ^ 2 +
      (crds_pn[, "y"] - crds_kitz[, "y"]) ^ 2
    )

    datediff <- abs(dd - igl$date[rows_kitz])

    out <- data.frame(id_kitz = ids_kitz,
                      distances = distances,
                      date_diff = datediff)

    # In all cases there is one fire at zero distance
    out <- out[which.min(out$distances), ]
    out$id_pn <- id

    return(out)
  } else {
    out <- data.frame(id_kitz = NA,
                      distances = 1e6,
                      date_diff = 100,
                      id_pn = id)
    return(out)
  }
}))

# Almost all lightning ignitions from PNNH are in Kitzberger database.
# Use the Kitzberger table. Only one fire is present at PNNH, so it will be added
# to Kitzberger table.

# for all fires in conflict_ligth with ditances == 0, set the date as in Bari
# database, and the area as an average of both sources.
for(i in 1:nrow(conflict_light)) {
  if(conflict_light$distances[i] < 10) {
    row_kitz <- which(igl$id == conflict_light$id_kitz[i])
    row_pn <- which(igpn$id == conflict_light$id_pn[i])
    igl$date[row_kitz] <- igpn$date[row_pn]
    aa <- mean(c(igl$area[row_kitz], igpn$area[row_pn]), na.rm = T)
    igl$area[row_kitz] <- aa
  }
}

# move to kitzberger table the unmatched Bari lightning ignitions
ids_move <- conflict_light$id_pn[is.na(conflict_light$id_kitz)]

igl <- rbind(igl, igpn[igpn$id %in% ids_move, ])
igpn <- igpn[!(igpn$id %in% ids_move), ]
# remove lightning-caused ignitions from pn database
igpn <- igpn[igpn$cause != "lightning", ]

### Try to match unknown cause ignitions in Bari database with lightning from
### Kitzberger.

ids_unk_pn <- igpn$id[igpn$cause == "unknown"]

conflict_unk <- do.call("rbind", lapply(ids_unk_pn, function(id) {
  # filter by date, searching in a 5 days radius
  ## id = igpn$id[13]
  i <- which(igpn$id == id)
  dd <- igpn$date[i]
  dseq <- seq(dd-5, dd+5, 1)

  rows_kitz <- which(igl$date %in% dseq)
  ids_kitz <- igl$id[rows_kitz]
  if(length(rows_kitz) > 0) {
    crds_pn <- crds(igpn[i, ])
    crds_kitz <- crds(igl[rows_kitz, ])

    distances <- sqrt(
      (crds_pn[, "x"] - crds_kitz[, "x"]) ^ 2 +
        (crds_pn[, "y"] - crds_kitz[, "y"]) ^ 2
    )

    datediff <- abs(dd - igl$date[rows_kitz])

    out <- data.frame(id_kitz = ids_kitz,
                      distances = distances,
                      date_diff = datediff)

    # In all cases there is one fire at zero distance
    out <- out[which.min(out$distances), ]
    out$id_pn <- id

    return(out)
  } else {
    out <- data.frame(id_kitz = NA,
                      distances = 1e6,
                      date_diff = 100,
                      id_pn = id)
    return(out)
  }
}))

## visualize the unknown
pnnh_widen <- buffer(pnnh, 500)
kitz_in <- relate(igl, pnnh_widen, relation = "within")[, 1]
bari_in <- relate(igpn, pnnh_widen, relation = "within")[, 1]

# plot(pnnh)
# plot(igl[kitz_in, ], add = T, col = 2, cex = 1)
# plot(igpn[igpn$cause == "human" & bari_in, ], add = T, col = "green", cex = 1)
# plot(igpn[igpn$cause == "unknown" & bari_in, ], add = T, col = "blue", cex = 0.8)
# # Most seem to be human-caused, but this is not so evident for some of them.


# Weight FWI pixels according to their intersection with PNNH -------------

fwi_fort <- crop(fwi_fort, pnnh, snap = "out")
template <- fwi_fort[[1]]
npix_fwi <- ncell(template)
values(template) <- 1:npix_fwi
pixels <- as.polygons(template, values = T, extent = F)
parts <- intersect(pixels, pnnh)
# plot(parts, col = 1:8)
parts_size <- expanse(parts)
pix_weights <- parts_size / sum(parts_size) # pixels by row!
study_area_size <- expanse(pnnh, "km")
# 7161.577 km2, that will be the unit for the ignition rate.


# Define temporal range ---------------------------------------------------

# range(igpn$date)
# range(igl$date)
time_human <- c(as.Date("1998-07-01"), as.Date("2021-06-30"))
time_light <- c(as.Date("1998-07-01"), as.Date("2023-06-30"))

# Tidy data for temporal model -------------------------------------------

# ig_count[fortnight(nfort), cause(3)]: number of ignitions by cause in columns
#   (human, lightning, unknown).
# fwi[nfort, npix]: standardized fwi in all pixels of the study area.
# weights[npix]: normalized weights or area

# The ignition rate will be at the spatial scale of the study area size

timefilter <- dtable$date >= min(c(time_human, time_light)) &
              dtable$date <= max(c(time_human, time_light))
time_data <- dtable[timefilter, ]
time_data <- time_data[!duplicated(time_data[, "fort", drop = F]), ]
# in which dates are both human and lightning ignitions?
time_data$cause_both <- 1
time_data$cause_both[time_data$date > max(time_human)] <- 0

# count of ignitions in the whole study area by cause
nfort <- nrow(time_data)
ig_counts <- matrix(0, nfort, 3)
rownames(ig_counts) <- time_data$fort
colnames(ig_counts) <- c("human", "lightning", "unknown")

# Apply a spatial and temporal filter to databases.
filter_pn <-
  igpn$date >= min(time_human) &
  igpn$date <= max(time_human) &
  bari_in
igpn_sub <- igpn[filter_pn, ]

filter_l <-
  igl$date >= min(time_light) &
  igl$date <= max(time_light) &
  kitz_in
igl_sub <- igl[filter_l, ]

# add fortnight to ignitions datasets and sort!
igpn_sub$fort <- date2fort(igpn_sub$date)
igl_sub$fort <- date2fort(igl_sub$date)
igpn_sub <- igpn_sub[order(igpn_sub$date), ]
igl_sub <- igl_sub[order(igl_sub$date), ]

# fill ig_counts
mm1 <- as.matrix(table(igpn_sub$fort, igpn_sub$cause))
mm2 <- as.matrix(table(igl_sub$fort))

ig_counts[rownames(mm1), "human"] <- mm1[, "human"]
ig_counts[rownames(mm1), "unknown"] <- mm1[, "unknown"]
ig_counts[rownames(mm2), "lightning"] <- mm2[, 1]

# FWI matrix (focal fortnight)
fwi_vals_all <- t(values(fwi_fort))
fwi_vals <- fwi_vals_all[as.character(time_data$fort), ]

# approximate mean and sd for standardization
fwi_vals_vec <- fwi_vals %*% pix_weights
fwi_mean <- mean(fwi_vals) # 6.340739
fwi_sd <- sd(fwi_vals)     # 8.907912

# focal values standardized and simplified
fwi_vals_all_z <- (fwi_vals_all - fwi_mean) / fwi_sd
fwi_vals_z <- (fwi_vals - fwi_mean) / fwi_sd
fwi_vals_vec_z <- (fwi_vals_vec - fwi_mean) / fwi_sd

# Create array with lagged values
nlags <- 11
lagtime <- 0:(nlags-1)
fwi_arr <- array(NA, dim = c(nfort, nlags, npix_fwi))
for(f in 1:nfort) {
  fort_focal <- time_data$fort[f]
  row_now <- which(rownames(fwi_vals_all_z) == as.character(fort_focal))
  rows <- row_now:(row_now-nlags+1)
  fwi_arr[f, , ] <- fwi_vals_all_z[rows, ]
}

# How to compute the likelihood for ignition counts?
condition_like <- rep(1, nfort) # evaluate both human and lightning-caused, known
condition_like[ig_counts[, 3] > 0] <- 2 # consider unknown cause
condition_like[time_data$cause_both == 0] <- 3 # known, but only lightning

# Try simple models to check whether modelling with overdispersion is necessary.

# library(mgcv)
# library(DHARMa)
# # likely human
# yy1 <- rowSums(ig_counts[, c("human", "unknown")])
# mm1 <- gam(yy1 ~ fwi_vals_vec_z, family = nb())
# plot(simulateResiduals(mm1)) # nb muy bien
# # likely lightning
# yy2 <- ig_counts[, "lightning"]
# mm2 <- gam(yy2 ~ fwi_vals_vec_z, family = nb())
# plot(simulateResiduals(mm2)) # nb perfecto

# Export points to sample spatial covariates -----------------------------

# merge all ignition data
igboth <- rbind(igpn_sub, igl_sub)
npoint <- nrow(igboth)
igboth$cause_num <- as.numeric(factor(igboth$cause,
                                      levels = c("human", "lightning", "unknown")))
nbycause <- as.numeric(table(igboth$cause_num))

# identify the corresponding row in the time_data for each ignition point, to
# match the expected rate of every ignition type.
igboth$row_rate <- NA
for(i in 1:npoint) {
  fortt <- igboth$fort[i]
  row <- which(time_data$fort == fortt)
  igboth$row_rate[i] <- row
}

# # write database to be uploaded to Earth engine, so we can sample the
# # spatial covariates.
# writeVector(project(igboth, "EPSG:4326"),
#             file.path(igdata_dir, "ignition_points_pnnh_bari-kitzberger.shp"),
#             overwrite = T)
# writeVector(project(pnnh, "EPSG:4326"),
#             file.path(igdata_dir, "pnnh.shp"),
#             overwrite = T) # used for the population sample


# Prepare data for spatial model ------------------------------------------

igspat <- vect(file.path(igdata_dir, "ignition_points_pnnh_bari-kitzberger_data.shp"))
igspat_flat <- project(igspat, "EPSG:5343") # match FWI projection
# will be used to resample FWI in escape model

igspat <- as.data.frame(igspat)

varnames <- c("veg_focal",
              "fwi",
              "elevation",
              "aspect",
              "slope",
              "dist_humans",
              "dist_roads")

varnames_ndvi <- sort(colnames(igspat)[grep("b_", colnames(igspat))])
colnames(igspat)[colnames(igspat) == "dist_human"] <- "dist_humans"

# merge with previous data
ig <- left_join(
  as.data.frame(igboth),
  igspat[, c("id", varnames, varnames_ndvi)],
  by = "id"
)
# nrow(ig) - nrow(igspat) # 3 missing points
ig <- ig[complete.cases(ig[, c(varnames, varnames_ndvi)]), ]

# sample of population of pixels
pop <- vect(file.path(igdata_dir, "population_points_pnnh_bari-kitzberger_data.shp"))
pop <- as.data.frame(pop)
colnames(pop)[colnames(pop) == "dist_human"] <- "dist_humans"
pop <- pop[, c(varnames, varnames_ndvi)]
pop <- pop[complete.cases(pop), ]

# # vegetation type of ignition points
# table(ig$veg_focal) # 5 urban
# ig$id[ig$veg_focal > 7] # to check in earth engine

# Correct vegetation type in a few points, by hand
ig$veg_focal[ig$id == "bari_5"] <- 6
ig$veg_focal[ig$id == "bari_35"] <- 5
ig$veg_focal[ig$id == "bari_37"] <- 10 # truely high andean

ig$veg_focal[ig$id == "bari_68"] <- 5
ig$veg_focal[ig$id == "bari_73"] <- 6
ig$veg_focal[ig$id == "bari_74"] <- 6

ig$veg_focal[ig$id == "bari_102"] <- 1
ig$veg_focal[ig$id == "bari_116"] <- 1
ig$veg_focal[ig$id == "bari_242"] <- 1

ig <- ig[ig$veg_focal < 8, ]
pop <- pop[pop$veg_focal < 8, ]

# recode vegetation type
ig <- left_join(ig, dveg[, c("veg_focal", "veg_num")], by = "veg_focal")
pop <- left_join(pop, dveg[, c("veg_focal", "veg_num")], by = "veg_focal")

# detrend NDVI
years_ndvi <- sapply(varnames_ndvi, function(x) strsplit(x, "_")[[1]][2]) |> as.numeric()
ny <- length(years_ndvi)

ig_ndvi <- as.matrix(ig[, varnames_ndvi])
ig_ndvi_dt <- matrix(NA, nrow(ig), ny)
ig_ndvi_dt[, ny] <- ig_ndvi[, ny]

pop_ndvi <- as.matrix(pop[, varnames_ndvi])
pop_ndvi_dt <- matrix(NA, nrow(pop), ny)
pop_ndvi_dt[, ny] <- pop_ndvi[, ny]

colnames(pop_ndvi_dt) <- colnames(ig_ndvi_dt) <- years_ndvi

# detrend NDVI, convert to 2022 equivalent
for(j in 1:(ny-1)) {
  # ignition points
  n22 <- ig_ndvi[, ny]
  nnow <- ig_ndvi[, j]
  dpred_ndvi <- data.frame(
    ndvi_dyn_logit = qlogis((nnow + 1) / 2),
    ndvi01_22 = (n22 + 1) / 2,
    year = years_ndvi[j]
  )
  dpred_ndvi$diff_logit <- mgcv::predict.gam(mdetrend, dpred_ndvi, se.fit = F)
  ig_ndvi_dt[, j] <- plogis(dpred_ndvi$ndvi_dyn_logit - dpred_ndvi$diff_logit) * 2 - 1

  # population points
  n22 <- pop_ndvi[, ny]
  nnow <- pop_ndvi[, j]
  dpred_ndvi <- data.frame(
    ndvi_dyn_logit = qlogis((nnow + 1) / 2),
    ndvi01_22 = (n22 + 1) / 2,
    year = years_ndvi[j]
  )
  dpred_ndvi$diff_logit <- mgcv::predict.gam(mdetrend, dpred_ndvi, se.fit = F)
  pop_ndvi_dt[, j] <- plogis(dpred_ndvi$ndvi_dyn_logit - dpred_ndvi$diff_logit) * 2 - 1
}
# plot(as.vector(ig_ndvi) ~ as.vector(ig_ndvi_dt))

# compute VFI
ig_vfi <- ig_ndvi   # placeholders
pop_vfi <- pop_ndvi

for(v in 1:5) {
  rows <- ig$veg_num == v
  ig_vfi[rows, ] <-
    fi_params$a[v] +
    fi_params$b[v] *
    (ig_ndvi_dt[rows, ] - fi_params$o[v]) ^ 2

  rows <- pop$veg_num == v
  pop_vfi[rows, ] <-
    fi_params$a[v] +
    fi_params$b[v] *
    (pop_ndvi_dt[rows, ] - fi_params$o[v]) ^ 2
}
# standardize
ig_vfi <- (ig_vfi - fi_params$vfi_mean) / fi_params$vfi_sd
pop_vfi <- (pop_vfi - fi_params$vfi_mean) / fi_params$vfi_sd

# Compute TFI (Topographic flammability index)

ig$northing <- cos(ig$aspect * pi / 180) * northing_weight(ig$slope)
ig$tfi <-
  fi_params$b_elev_ori * ig$elevation +
  fi_params$b_north_ori * ig$northing
ig$tfi <- (ig$tfi - fi_params$tfi_mean) / fi_params$tfi_sd

pop$northing <- cos(pop$aspect * pi / 180) * northing_weight(pop$slope)
pop$tfi <-
  fi_params$b_elev_ori * pop$elevation +
  fi_params$b_north_ori * pop$northing
pop$tfi <- (pop$tfi - fi_params$tfi_mean) / fi_params$tfi_sd

# Assign the correct VFI to ignitions, using the one computed for its previous
# fire-year.
ig <- left_join(ig, dtable[, c("date", "year")], by = "date")
ig$vfi <- NA
for(i in 1:nrow(ig)) {
  column <- ig$year[i] - 1 - 1997 # because ndvi years start at 1998
  ig$vfi[i] <- ig_vfi[i, column]
}

# subset npop pixels to use
npop <- 10000; set.seed(23542896)
ids_pop <- sample(1:nrow(pop), npop, replace = F)

# turn distances to km and standardize
ig$dist_humans_km <- ig$dist_humans / 1000
ig$dist_roads_km <- ig$dist_roads / 1000
pop$dist_humans_km <- pop$dist_humans / 1000
pop$dist_roads_km <- pop$dist_roads / 1000

dh_mean <- mean(pop$dist_humans_km); dh_sd <- sd(pop$dist_humans_km)
dr_mean <- mean(pop$dist_roads_km); dr_sd <- sd(pop$dist_roads_km)

ig$dhz <- (ig$dist_humans_km - dh_mean) / dh_sd
ig$drz <- (ig$dist_roads_km - dr_mean) / dr_sd
pop$dhz <- (pop$dist_humans_km - dh_mean) / dh_sd
pop$drz <- (pop$dist_roads_km - dr_mean) / dr_sd

# standardize spatial fwi
fwi_spat_mean <- mean(pop$fwi)
fwi_spat_sd <- sd(pop$fwi)
pop$fwi_z <- (pop$fwi - fwi_spat_mean) / fwi_spat_sd
ig$fwi_z <- (ig$fwi - fwi_spat_mean) / fwi_spat_sd

# Design matrices for ignition points
formula_fi <- ~ vfi + tfi - 1
formula_dist <- ~ dhz + drz - 1

X_ig_fi <- model.matrix(formula_fi, ig)
X_ig_dist <- model.matrix(formula_dist, ig)

X_pop_dist <- model.matrix(formula_dist, pop)
X_pop_fi <- lapply(1:(ny-1), function(y) {
  cbind(pop_vfi[, y], pop$tfi)
  # ig years are 1999:2022; ndvi years are 1998:2022
  # so the vfi is computed using the y-1 NDVI (1998:2021)
})
X_pop_fi <- abind::abind(X_pop_fi, along = 3)
dimnames(X_pop_fi) <- list("row" = 1:nrow(X_pop_fi),
                           "var" = c("vfi", "tfi"),
                           "year" = min(ig$year):max(ig$year))

# number of years in the study period
ny2 <- length(min(ig$year):max(ig$year))
ig$year_id <- ig$year - min(ig$year) + 1 # to identify which VFI to use
ig$fort_id <- ig$row_rate # better naming

# # Check vfi
# plot(density(X_pop_fi[, 1, 1], n = 2^10), xlim = c(-5, 2), ylim = c(0, 2))
# for(v in 2:ny) {
#   lines(density(X_pop_fi[, 1, v], n = 2^10), col = v)
# }

# Ignition model -------------------------------------------------

# compute average ignition rate to set a non-influential prior on the median
ig_transient <- ig_counts[, 1:2]
ig_transient[, 1] <- ig_transient[, 1] + ig_counts[, 3]

sdata <- list(
  ## ignition rate model
  nfort = nfort, npix = npix_fwi, nlag = nlags, nig = 2,
  ig_counts = ig_counts,
  condition = condition_like,
  fwi = aperm(fwi_arr, c(3, 1, 2)),
  pix_weights = pix_weights,
  nfort_both = sum(condition_like < 3),

  prior_a_sd = 5,
  prior_a_mean = log(mean(ig_transient)),
  prior_b_sd = 5,
  prior_phi_sd = 50,
  prior_ls_sd = nlags * 0.75,

  ## ignition location model
  npoint = nrow(ig), nland = npop,
  ny = max(ig$year_id), nfi = 2, ndist = 2,

  cause = ig$cause_num,
  fort_id = ig$fort_id,
  year_id = ig$year_id,

  X_ig_fi = X_ig_fi,
  # ig_fwi = ig$fwi_z,
  X_ig_dist = X_ig_dist,
  X_pop_fi = aperm(X_pop_fi[ids_pop, , ], c(3, 1, 2)),
  # pop_fwi = pop$fwi_z[ids_pop],
  X_pop_dist = X_pop_dist[ids_pop, ],

  prior_c_sd = 10
)

# smodel <- stan_model(file.path("ignition", "ignition_model3.stan"))
# igmod <- sampling(
#   smodel, data = sdata, seed = 85963, refresh = 20,
#   # cores = 1, chains = 1, iter = 5,
#   cores = 8, chains = 8, iter = 2000, warmup = 1000,
#   pars = c("a", "b", "phi", "ls",
#            "c_fi", "c_dist", "human_prop",
#            "c_vfi_h", "c_vfi_l",  # separate c parameters
#            "c_tfi_h", "c_tfi_l",
#            # "c_fwi_h", "c_fwi_l",
#            "c_dist_h", "c_dist_r")
# )
# saveRDS(igmod, file.path("files", "ignition", "ignition_model_samples.rds"))
# # 2382.66 / 60 = 40 min

# Model samples named *_samples2.rds includes FWI effect. It was not used because
# it estimates a huge effect of FWI on human ignitions (almost null with
# lightning). This is due to points at the eastern limit of the study area,
# where fwi is really high, and little population points fall there.

igmod <- readRDS(file.path("files", "ignition", "ignition_model_samples.rds"))
sig <- summary(igmod)[[1]]
min(sig[, "n_eff"]) # 1636.857
max(sig[, "Rhat"])  # 1.004732

pairs(igmod, pars = c("c_vfi_h", "c_vfi_l", "c_tfi_h", "c_tfi_l"))
pairs(igmod, pars = c("c_dist_r", "c_dist_h")) #"c_fwi_h", "c_fwi_l",

# glimpse at the posteriors
mcmc_dens(igmod, pars = c("a[1]", "a[2]"),
          facet_args = list(scales = "fixed", nrow = 2))
mcmc_dens(igmod, pars = c("b[1]", "b[2]"),
          facet_args = list(scales = "fixed", nrow = 2))
mcmc_dens(igmod, pars = c("phi[1]", "phi[2]"),
          facet_args = list(scales = "fixed", nrow = 2)) + xlim(0, 8)
mcmc_dens(igmod, pars = "ls") # biutiful, poca corr temporal

mcmc_dens(igmod, pars = c("c_vfi_h", "c_vfi_l", "c_tfi_h", "c_tfi_l"),
          facet_args = list(scales = "fixed", nrow = 2))
# mcmc_dens(igmod, pars = c("c_fwi_h", "c_fwi_l"),
#           facet_args = list(scales = "fixed", nrow = 2))
mcmc_dens(igmod, pars = c("c_dist_r", "c_dist_h"),
          facet_args = list(scales = "fixed", nrow = 2))

mcmc_dens(igmod, pars = "human_prop") + xlim(0, 1)
rows_compare <- igboth$date <= max(time_human)
tt <- table(igboth$cause[rows_compare])
pnorm <- tt / sum(tt)
humanp_obs <- 1 - pnorm[2]
humanp_obs2 <- pnorm[1] / sum(pnorm)
plot(density(as.matrix(igmod, "human_prop"), from = 0, to = 1, n = 2 ^ 10),
     main = NA, xlab = "human proportion", ylim = c(0, 18))
# lines(density(as.matrix(igmod2, "humanp"), from = 0, to = 1, n = 2 ^ 10),
#       col = "red")
abline(v = humanp_obs, lty = 2)
abline(v = humanp_obs2, lty = 2, col = 2)
# lines show the proportions if all unknown were lightning (red) or human

# Ignition model predictions ----------------------------------------------

# FWI weights ____________________________________________________________

ls_hat <- as.numeric(as.matrix(igmod, "ls"))
npost <- length(ls_hat)

tseq <- seq(0, nlags-1, by = 0.1)
np <- length(tseq)
ids_norm <- seq(1, np, by = 10)

wun <- matrix(NA, np, npost)
wn <- wun

for(i in 1:npost) {
  wun[, i] <- exp(-(tseq / ls_hat[i]) ^ 2)
  norm_const <- sum(wun[ids_norm, i])
  wn[, i] <- wun[, i] / norm_const # normalize only at discrete timepoints
}

plot(wun[, 1] ~ tseq, type = "l")
for(i in 1:100) lines(wun[, i+1] ~ tseq, col = rgb(0, 0, 0, 0.05))

plot(wn[, 1] ~ tseq, type = "l")
for(i in 1:100) lines(wn[, i+1] ~ tseq, col = rgb(0, 0, 0, 0.05))

# summarize
wsumm <- apply(wn, 1, mean_ci) |> t() |> as.data.frame()
wsumm$t <- tseq

ggplot(wsumm, aes(t, mean, ymin = lower, ymax = upper)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  geom_point(data = wsumm[ids_norm, ], size = 2, alpha= 0.8) +
  scale_x_continuous(breaks = seq(0, 10, 2),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1),
                     expand = c(0.01, 0.01)) +
  ylab("FWI weight") +
  xlab("Fortnights before fire") +
  nice_theme()

nn <- file.path("ignition", "figures", "ignition_fwi_weights.png")
ggsave(nn, width = 10, height = 9, units = "cm")


# Ig rate ~ FWI __________________________________________________________

# To create the FWIseq, use the posterior mean of lengthscale:
lsmean <- mean(ls_hat)
wpoint <- exp(-((0:10) / lsmean) ^ 2) |> normalize()

# Compute avg FWI for all the time sequence and get limits.

# reduce time:
fwi_summ0 <- fwi_arr[, 1, ]
for(j in 1:dim(fwi_arr)[3]) { # loop over pixels
  fwi_summ0[, j] <- fwi_arr[, , j] %*% wpoint # weighted average of previous
                                               # fortnights, matrix mult
}
# reduce pixels:
fwi_summ <- fwi_summ0 %*% pix_weights

np <- 200
fwiseq <- seq(min(fwi_summ), max(fwi_summ), length.out = np)
fwiseq_ori <- fwiseq * fwi_sd + fwi_mean

# alpha and beta for humans (h) and ligthning (l)
abh <- as.matrix(igmod, pars = c("a[1]", "b[1]")) |> t()
abl <- as.matrix(igmod, pars = c("a[2]", "b[2]")) |> t()
X <- cbind(rep(1, np), fwiseq)

lamh <- apply(exp(X %*% abh) / study_area_size, 1, mean_ci) |> t() |> as.data.frame()
laml <- apply(exp(X %*% abl) / study_area_size, 1, mean_ci) |> t() |> as.data.frame()
# scale to ignition rate / (fortnight, km2)

lamh$Cause <- "Human"
laml$Cause <- "Lightning"

lam <- rbind(lamh, laml)
lam$Cause <- factor(lam$Cause, levels = c("Human", "Lightning"))
lam$fwi <- rep(fwiseq_ori, 2)

ggplot(lam, aes(fwi, mean, ymin = lower, ymax = upper,
                color = Cause, fill = Cause, group = Cause)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  scale_color_viridis(discrete = T, end = 0.4) +
  scale_fill_viridis(discrete = T, end = 0.4) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(1e-6, 1e-6)) +
  ylab("Ignition rate (fires / (fortnight * km2))") +
  xlab("Fire Weather Index") +
  nice_theme()

nn <- file.path("ignition", "figures", "ignition_rate_fwi.png")
ggsave(nn, width = 12, height = 9, units = "cm")


# Ignition location functions ______________________________________________

# Ig probability (conditional logistic) as a function of

# 1) NDVI by veg type
# 2) elevation
# 3) northing
# 4) distance from roads and from settlements (different color)

# in 2 columns [human, lightning]

# (borrowing code from spread hierarchical model)

## NDVI _________________________

veg_levels <- data_summ$ndvi$vegetation
nveg <- length(veg_levels)

nseq <- 300

pd_ndvi <- do.call("rbind", lapply(1:nveg, function(v) {
  expand.grid(
    ndvi = seq(data_summ$ndvi$hdi_lower_95[v],
               data_summ$ndvi$hdi_upper_95[v],
               length.out = nseq),
    vegnum = v,
    vegetation = veg_levels[v]
  )
}))

pd_ndvi$vfi <- 0

# VFI
for(v in 1:nveg) {
  id_fill <- pd_ndvi$vegnum == v
  pd_ndvi$vfi[id_fill] <-
    fi_params$a[v] +
    fi_params$b[v] *
    (pd_ndvi$ndvi[id_fill] - fi_params$o[v]) ^ 2
}
# (standardize)
pd_ndvi$vfi <- (pd_ndvi$vfi - fi_params$vfi_mean) / fi_params$vfi_sd


mh <- matrix(NA, nrow(pd_ndvi), npost)
ml <- matrix(NA, nrow(pd_ndvi), npost)

# compute unnormalized probabilities
mh <- outer(pd_ndvi$vfi, as.numeric(as.matrix(igmod, "c_vfi_h"))) |> exp()
ml <- outer(pd_ndvi$vfi, as.numeric(as.matrix(igmod, "c_vfi_l"))) |> exp()

# normalize
for(i in 1:npost) {
  for(v in 1:nveg) {
    id_fill <- pd_ndvi$vegnum == v
    mh[id_fill, i] <- normalize(mh[id_fill, i])
    ml[id_fill, i] <- normalize(ml[id_fill, i])
  }
}

# summarize and plot
summh <- apply(mh, 1, mean_ci) |> t() |> as.data.frame()
summl <- apply(ml, 1, mean_ci) |> t() |> as.data.frame()
summh$Cause <- "Human"
summl$Cause <- "Lightning"

prob_ndvi <- cbind(
  rbind(pd_ndvi, pd_ndvi),
  rbind(summh, summl)
)

bgcol <- "#1a1a1aff"

cond_ndvi <-
  ggplot(prob_ndvi, aes(ndvi, mean, ymin = lower, ymax = upper,
                       color = vegetation, fill = vegetation)) +
  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +
  scale_color_viridis(discrete = T, begin = 0, end = 0.9, option = "D",
                      name = "Vegetation\ntype") +
  scale_fill_viridis(discrete = T, begin = 0, end = 0.9, option = "D",
                     name = "Vegetation\ntype") +
  facet_wrap(vars(Cause), ncol = 2, axes = "all",
             axis.labels = "margins") +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(1e-5, 1e-5),
                     labels = function(x) x * 100) +
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  xlab("NDVI") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(b = 5, unit = "mm"),
        strip.background = element_rect(fill = bgcol, color = bgcol),
        strip.text = element_text(color = "white"))
cond_ndvi

## Vegetation _________________________

pd_veg <- data_summ$ndvi[, 1:2]
colnames(pd_veg)[2] <- "ndvi"
pd_veg$vfi <- 0

# VFI
for(v in 1:nveg) {
  pd_veg$vfi[v] <-
    fi_params$a[v] +
    fi_params$b[v] *
    (pd_veg$ndvi[v] - fi_params$o[v]) ^ 2
}
# (standardize)
pd_veg$vfi <- (pd_veg$vfi - fi_params$vfi_mean) / fi_params$vfi_sd

mh <- matrix(NA, nrow(pd_veg), npost)
ml <- matrix(NA, nrow(pd_veg), npost)

# compute unnormalized probabilities
mh <- outer(pd_veg$vfi, as.numeric(as.matrix(igmod, "c_vfi_h"))) |> exp()
ml <- outer(pd_veg$vfi, as.numeric(as.matrix(igmod, "c_vfi_l"))) |> exp()

# normalize
for(i in 1:npost) {
  mh[, i] <- normalize(mh[, i])
  ml[, i] <- normalize(ml[, i])
}

# summarize and plot
summh <- apply(mh, 1, mean_ci) |> t() |> as.data.frame()
summl <- apply(ml, 1, mean_ci) |> t() |> as.data.frame()
summh$Cause <- "Human"
summl$Cause <- "Lightning"

prob_veg <- cbind(
  rbind(pd_veg, pd_veg),
  rbind(summh, summl)
)

cond_veg <-
  ggplot(prob_veg, aes(vegetation, mean, ymin = lower, ymax = upper,
                       fill = vegetation)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2,
           alpha = 0.7) +
  geom_linerange() +
  # scale_color_viridis(discrete = T, begin = 0, end = 0.9, option = "D",
  #                     name = "Vegetation\ntype") +
  scale_fill_viridis(discrete = T, begin = 0, end = 0.9, option = "D",
                     name = "Vegetation\ntype") +
  facet_wrap(vars(Cause), ncol = 2, axes = "all",
             axis.labels = "margins") +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 0),
                     labels = function(x) x * 100) +
  # scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  xlab("Vegetation type") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        panel.spacing.x = unit(4, "mm"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())
cond_veg


## Topography _________________________

pd_topo <- rbind(
  expand.grid(
    elevation = seq(data_summ$elevation["hdi_lower_95"],
                    data_summ$elevation["hdi_upper_95"],
                    length.out = nseq),
    northing = 0,
    varying_var = "elevation"
  ),
  expand.grid(
    elevation = data_summ$elevation["mean"],
    northing = seq(-1, 1, length.out = nseq),
    varying_var = "northing"
  )
)

# Compute tfi
pd_topo$tfi <-
  fi_params$b_elev_ori * pd_topo$elevation +
  fi_params$b_north_ori * pd_topo$northing
pd_topo$tfi <- (pd_topo$tfi - fi_params$tfi_mean) / fi_params$tfi_sd

# compute unnormalized probabilities
mh <- outer(pd_topo$tfi, as.numeric(as.matrix(igmod, "c_tfi_h"))) |> exp()
ml <- outer(pd_topo$tfi, as.numeric(as.matrix(igmod, "c_tfi_l"))) |> exp()

# normalize
for(i in 1:npost) {
  for(v in c("elevation", "northing")) {
    id_fill <- pd_topo$varying_var == v
    mh[id_fill, i] <- normalize(mh[id_fill, i])
    ml[id_fill, i] <- normalize(ml[id_fill, i])
  }
}

# summarize and plot
summh <- apply(mh, 1, mean_ci) |> t() |> as.data.frame()
summl <- apply(ml, 1, mean_ci) |> t() |> as.data.frame()
summh$Cause <- "Human"
summl$Cause <- "Lightning"

prob_topo <- cbind(
  rbind(pd_topo, pd_topo),
  rbind(summh, summl)
)

cond_elev <-
  ggplot(prob_topo[prob_topo$varying_var == "elevation", ],
         aes(elevation, mean, ymin = lower, ymax = upper)) +
  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +
  facet_wrap(vars(Cause), ncol = 2, axes = "all",
             axis.labels = "margins") +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(1e-5, 1e-5),
                     labels = function(x) x * 100) +
  xlab("Elevation (m a.s.l.)") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        strip.background = element_rect(fill = bgcol, color = bgcol),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(b = 5, unit = "mm"))
cond_elev

cond_north <-
  ggplot(prob_topo[prob_topo$varying_var == "northing", ],
         aes(northing, mean, ymin = lower, ymax = upper)) +
  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +
  facet_wrap(vars(Cause), ncol = 2, axes = "all",
             axis.labels = "margins") +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(1e-5, 1e-5),
                     labels = function(x) x * 100) +
  xlab("Northing") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        strip.text = element_blank(),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(b = 5, unit = "mm"))
cond_north


## FWI spatial _______________________

# pd_fwi <- data.frame(
#   fwi = seq(0, quantile(pop$fwi, prob = 0.98), length.out = nseq)
# )
# pd_fwi$fwiz <- (pd_fwi$fwi - fwi_spat_mean) / fwi_spat_sd
#
# # compute unnormalized probabilities
# mh <- outer(pd_fwi$fwiz, as.numeric(as.matrix(igmod, "c_fwi_h"))) |> exp()
# ml <- outer(pd_fwi$fwiz, as.numeric(as.matrix(igmod, "c_fwi_l"))) |> exp()
#
# # normalize
# for(i in 1:npost) {
#   mh[, i] <- normalize(mh[, i])
#   ml[, i] <- normalize(ml[, i])
# }
#
# # summarize and plot
# summh <- apply(mh, 1, mean_ci) |> t() |> as.data.frame()
# summl <- apply(ml, 1, mean_ci) |> t() |> as.data.frame()
# summh$Cause <- "Human"
# summl$Cause <- "Lightning"
#
# prob_fwi <- cbind(
#   rbind(pd_fwi, pd_fwi),
#   rbind(summh, summl)
# )
#
# cond_fwi <-
#   ggplot(prob_fwi,
#          aes(fwi, mean, ymin = lower, ymax = upper)) +
#   geom_ribbon(color = NA, alpha = 0.4) +
#   geom_line() +
#   facet_wrap(vars(Cause), ncol = 2, axes = "all",
#              axis.labels = "margins") +
#   expand_limits(y = 0) +
#   scale_y_continuous(expand = c(1e-5, 1e-5),
#                      labels = function(x) x * 100) +
#   scale_x_continuous(expand = c(0.01, 0.01)) +
#   xlab("Fire Weather Index") +
#   nice_theme() +
#   theme(axis.title.y = element_blank(),
#         strip.text = element_blank(),
#         panel.spacing.x = unit(4, "mm"),
#         plot.margin = margin(b = 5, unit = "mm"))
# cond_fwi

## Distance from roads _______________________

pd_roads <- data.frame(
  dist_roads = seq(0, quantile(pop$dist_roads_km, prob = 0.98), length.out = nseq)
)
pd_roads$drz <- (pd_roads$dist_roads - dr_mean) / dr_sd

# compute unnormalized probabilities
mh <- outer(pd_roads$drz, as.numeric(as.matrix(igmod, "c_dist_r"))) |> exp()

# normalize
for(i in 1:npost) {
  mh[, i] <- normalize(mh[, i])
}

# summarize and plot
summh <- apply(mh, 1, mean_ci) |> t() |> as.data.frame()
summl <- summh
summl[,] <- 1/nseq
summh$Cause <- "Human"
summl$Cause <- "Lightning"

prob_roads <- cbind(
  rbind(pd_roads, pd_roads),
  rbind(summh, summl)
)

cond_roads <-
  ggplot(prob_roads,
         aes(dist_roads, mean, ymin = lower, ymax = upper)) +
  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +
  facet_wrap(vars(Cause), ncol = 2, axes = "all",
             axis.labels = "margins") +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(1e-5, 1e-5),
                     labels = function(x) x * 100) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  xlab("Distance from roads (km)") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        strip.text = element_blank(),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(b = 5, unit = "mm"))
cond_roads

## Distance from humans _______________________

pd_humans <- data.frame(
  dist_humans = seq(0, quantile(pop$dist_humans_km, prob = 0.98),
                    length.out = nseq)
)
pd_humans$dhz <- (pd_humans$dist_humans - dh_mean) / dh_sd

# compute unnormalized probabilities
mh <- outer(pd_humans$dhz, as.numeric(as.matrix(igmod, "c_dist_h"))) |> exp()

# normalize
for(i in 1:npost) {
  mh[, i] <- normalize(mh[, i])
}

# summarize and plot
summh <- apply(mh, 1, mean_ci) |> t() |> as.data.frame()
summl <- summh
summl[,] <- 1/nseq
summh$Cause <- "Human"
summl$Cause <- "Lightning"

prob_humans <- cbind(
  rbind(pd_humans, pd_humans),
  rbind(summh, summl)
)

cond_humans <-
  ggplot(prob_humans,
         aes(dist_humans, mean, ymin = lower, ymax = upper)) +
  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +
  facet_wrap(vars(Cause), ncol = 2, axes = "all",
             axis.labels = "margins") +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(1e-5, 1e-5),
                     labels = function(x) x * 100) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  xlab("Distance from human settlements (km)") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        strip.text = element_blank(),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(b = 0, unit = "mm"))
cond_humans

## Merge all plots ______________

ytit1 <- grid::textGrob("Relative ignition probability (%)",
                        rot = 90, hjust = 0.25,
                        gp = grid::gpar(fontsize = 12, fontface = 'plain'))

condplot_veg <- ggarrange2(cond_ndvi, cond_veg,
                           ncol = 1, left = ytit1, padding = unit(8, "mm"))
ggsave("ignition/figures/ignition_conditional_prob_veg.png",
       width = fig_width_max, height = fig_height_max * 0.55, units = "cm",
       plot = condplot_veg)

ytit2 <- grid::textGrob("Relative ignition probability (%)",
                        rot = 90, hjust = 0.5,
                        gp = grid::gpar(fontsize = 12, fontface = 'plain'))

condplot_env <- ggarrange2(cond_elev, cond_north, #cond_fwi,
                           cond_roads, cond_humans,
                           ncol = 1, left = ytit2, padding = unit(8, "mm"))

ggsave("ignition/figures/ignition_conditional_prob_environment.png",
       width = fig_width_max * 0.8, height = fig_height_max * 0.8, units = "cm",
       plot = condplot_env)


# Prepare data for escape model -------------------------------------------

# Fire escape is defined as an ignition burning more than one 30 * 30 m pixel
# (> 0.09 ha). I will model this as a Bernoulli process, depending on
# FWI (with cumulative effect),
# vfi, tfi, distance to roads and distance to human settlements.
# FWI will be interpolated to each point to increase precision.

ig$escape <- numeric(nrow(ig))
for(i in 1:nrow(ig)) {
  if(!is.na(ig$area[i])) ig$escape[i] <- as.numeric(ig$area[i] > 0.09)
}
mean(ig$escape) # 0.23

# make equivalent spatial data to interpolate FWI to points
ig2 <- ig[order(ig$id), ]
igspat2 <- igspat_flat[order(igspat$id), ]
igspat2 <- igspat2[igspat2$id %in% ig2$id, ]
all.equal(ig2$id, igspat2$id)

# Get FWI time-series for each ignition (lagged values)
# (extract() does not interpolate!)
fwi_points <- matrix(NA, nrow(ig2), nlags)
colnames(fwi_points) <- 0:(-nlags+1)
for(i in 1:nrow(ig2)) {
  # i = 1
  ff <- ig2$fort[i]
  ffseq <- as.character(ff:(ff-(nlags-1)))
  ext <- extract(fwi_fort[[ffseq]], igspat2[i, ])[, -1] |> as.matrix()
  fwi_points[i, ] <- ext
}
fwi_points_z <- (fwi_points - fwi_mean) / fwi_sd


# Escape model (> 0.09 ha) ------------------------------------------------

# stan data
sdata_esc <- list(
  n = nrow(ig2),
  nlag = nlags,
  y = ig2$escape,
  fwi_mat_z = fwi_points_z,
  vfi = ig2$vfi,
  tfi = ig2$tfi,
  drz = ig2$drz,
  dhz = ig2$dhz,

  prior_a_sd = 10,
  prior_b_sd = 10,
  prior_ls_sd = nlags * 0.75
)

smodel_esc <- stan_model(file.path("ignition", "escape_model.stan"))
escmod <- sampling(
  smodel_esc, data = sdata_esc, seed = 1596142, refresh = 200,
  # cores = 1, chains = 1, iter = 5
  cores = 8, chains = 8, iter = 2000, warmup = 1000
)
saveRDS(escmod, file.path("files", "ignition", "escape_model_samples.rds"))

escmod <- readRDS(file.path("files", "ignition", "escape_model_samples.rds"))
sesc <- summary(escmod)[[1]]
min(sesc[, "n_eff"], na.rm = T) # 2730.538
max(sesc[, "Rhat"], na.rm = T)  # 1.004778

pairs(escmod, pars = c("a", "ls", "b_fwi", "b_vfi", "b_tfi", "b_drz", "b_dhz"))

# glimpse at the posteriors
mcmc_dens(escmod, pars = c("a", "ls", "b_fwi", "b_vfi", "b_tfi",
                           "b_drz", "b_dhz"),
          facet_args = list(scales = "free", ncol = 3))


# Escape model (fire size) -----------------------------------------------

ig2$area_stan <- ig2$area
ig2$area_stan[is.na(ig2$area)] <- 1e6

# stan data
sdata_size <- list(
  n = nrow(ig2),
  nlag = nlags,

  y = log(ig2$area_stan),
  cutoff = log(0.09), # one pixel size

  n_notna = sum(!is.na(ig2$area)),
  n_na = sum(is.na(ig2$area)),
  ids_na = which(is.na(ig2$area)),
  ids_notna = which(!is.na(ig2$area)),

  fwi_mat_z = fwi_points_z,
  vfi = ig2$vfi,
  tfi = ig2$tfi,
  drz = ig2$drz,
  dhz = ig2$dhz,

  prior_a_mean = mean(log(na.omit(ig2$area_stan))),
  prior_a_sd = 20,
  prior_b_sd = 5,
  prior_sigma_sd = 10,
  prior_g_sd = 10,
  prior_ls_sd = nlags * 0.75
)

smodel_size <- stan_model(file.path("ignition", "size_model.stan"))
sizemod <- sampling(
  smodel_size, data = sdata_size, seed = 1596142, refresh = 200,
  control = list(adapt_delta = 0.98),
  # cores = 1, chains = 1, iter = 5
  cores = 8, chains = 8, iter = 2000, warmup = 1000
)
saveRDS(sizemod, file.path("files", "ignition", "size_model_samples.rds"))

sizemod <- readRDS(file.path("files", "ignition", "size_model_samples.rds"))
ssize <- summary(sizemod)[[1]]
min(ssize[, "n_eff"], na.rm = T) # 2049.183
max(ssize[, "Rhat"], na.rm = T)  # 1.004151

pairs(sizemod, pars = c("a", "b_fwi", "b_vfi", "b_tfi", "b_drz", "b_dhz",
                        "ls", "sigma", "g"))

# glimpse at the posteriors
mcmc_dens(sizemod, pars = c("a", "b_fwi", "b_vfi", "b_tfi", "b_drz", "b_dhz",
                            "ls", "sigma", "g"),
          facet_args = list(scales = "free", ncol = 3))


# Fire Escape model (size class) -----------------------------------------

ig2$area_impute <- ig2$area
ig2$area_impute[is.na(ig2$area)] <- 0.045 # NA is less than 1pix

area_cuts <- c(0.09, 10, 100, max(ig2$area, na.rm = T) * 10)
K = length(area_cuts)

ig2$sizeclass <- sapply(area_cuts, function(cut) {
  as.numeric(ig2$area_impute > cut)
}) |> rowSums() + 1
table(ig2$sizeclass)

# stan data
sdata_sizeclass <- list(
  n = nrow(ig2),
  nlag = nlags,
  K = max(ig2$sizeclass),

  y = ig2$sizeclass,

  fwi_mat_z = fwi_points_z,
  vfi = ig2$vfi,
  tfi = ig2$tfi,
  drz = ig2$drz,
  dhz = ig2$dhz,

  prior_a_sd = 3,
  prior_b_sd = 10,
  prior_ls_sd = nlags * 0.75
)

smodel_sizeclass <- stan_model(file.path("ignition", "sizeclass_model.stan"))
scmod <- sampling(
  smodel_sizeclass, data = sdata_sizeclass, seed = 1596142, refresh = 200,
  # cores = 1, chains = 1, iter = 5,
  cores = 8, chains = 8, iter = 2000, warmup = 1000,
  pars = c("a", "b_fwi", "b_vfi", "b_tfi", "b_drz", "b_dhz", "ls", "pmf")
)
saveRDS(scmod, file.path("files", "ignition", "sizeclass_model_samples.rds"))

scmod <- readRDS(file.path("files", "ignition", "sizeclass_model_samples.rds"))
sscmod <- summary(scmod)[[1]]
min(sscmod[, "n_eff"], na.rm = T) # 2650.989
max(sscmod[, "Rhat"], na.rm = T)  # 1.002762

pairs(scmod, pars = c("a", "b_fwi", "b_vfi", "b_tfi", "b_drz", "b_dhz",
                             "ls"))

# glimpse at the posteriors
mcmc_dens(scmod, pars = c("a[1]", "a[2]", "a[3]",
                                 "b_fwi", "b_vfi", "b_tfi",
                                 "b_drz", "b_dhz", "ls"),
          facet_args = list(scales = "free", ncol = 3))

# Escape model predictions: vegetation -----------------------------------

# en función del NDVI: un panel por veg type.
# En función de la veg: distribución de clase por veg type.
# es un cbind de 2 facet_wraps, con legend abajo. Aunque el graf de abajo podría
# resolverse con position_dodge.

# Las demás covariables tendrían un panel cada una, con una curva para la prob
# de cada clase.

# Usar viridis C para las fire size classes.

ahat <- as.matrix(scmod, "a")
npost <- nrow(ahat)
class_names <- c("(0, 0.09]", "(0.09, 10]", "(10, 100]", "(100, ...)")
bgcol <- "#1a1a1aff"

## NDVI _________________________

veg_levels <- data_summ$ndvi$vegetation
nveg <- length(veg_levels)

nseq <- 300


pd_ndvi <- do.call("rbind", lapply(1:nveg, function(v) {
  expand.grid(
    ndvi = seq(data_summ$ndvi$hdi_lower_95[v],
               data_summ$ndvi$hdi_upper_95[v],
               length.out = nseq),
    vegnum = v,
    vegetation = veg_levels[v]
  )
}))

pd_ndvi$vfi <- 0

# VFI
for(v in 1:nveg) {
  id_fill <- pd_ndvi$vegnum == v
  pd_ndvi$vfi[id_fill] <-
    fi_params$a[v] +
    fi_params$b[v] *
    (pd_ndvi$ndvi[id_fill] - fi_params$o[v]) ^ 2
}
# (standardize)
pd_ndvi$vfi <- (pd_ndvi$vfi - fi_params$vfi_mean) / fi_params$vfi_sd

prob_ndvi <- ordinal_predict(pd_ndvi, "vfi", "b_vfi")

esc_ndvi <-
  ggplot(prob_ndvi, aes(ndvi, mean, ymin = lower, ymax = upper,
                        color = class_name, fill = class_name)) +
  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +
  scale_color_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                      name = "Fire size\nclass (ha)") +
  scale_fill_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                     name = "Fire size\nclass (ha)") +
  facet_wrap(vars(vegetation), ncol = 1, axes = "all",
             axis.labels = "margins", strip.position = "right") +
  expand_limits(y = c(0, 1)) +
  scale_y_continuous(expand = c(1e-5, 1e-5),
                     labels = function(x) x * 100) +
  # scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  xlab("NDVI") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(l = 5, unit = "mm"),
        strip.background = element_rect(fill = bgcol, color = bgcol),
        strip.text = element_text(color = "white"))
esc_ndvi


## Vegetation _________________________

pd_veg <- data_summ$ndvi[, 1:2]
colnames(pd_veg)[2] <- "ndvi"
pd_veg$vfi <- 0

# VFI
for(v in 1:nveg) {
  pd_veg$vfi[v] <-
    fi_params$a[v] +
    fi_params$b[v] *
    (pd_veg$ndvi[v] - fi_params$o[v]) ^ 2
}
# (standardize)
pd_veg$vfi <- (pd_veg$vfi - fi_params$vfi_mean) / fi_params$vfi_sd

prob_veg <- ordinal_predict(pd_veg, "vfi", "b_vfi")

esc_veg <-
  ggplot(prob_veg, aes(class_name, mean, ymin = lower, ymax = upper,
                        fill = class_name)) +
  geom_bar(stat = "identity", linewidth = 0.2,
           alpha = 0.7) +
  geom_linerange() +
  scale_color_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                      name = "Fire size\nclass (ha)") +
  scale_fill_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                     name = "Fire size\nclass (ha)") +
  facet_wrap(vars(vegetation), ncol = 1, axes = "all",
             axis.labels = "margins", strip.position = "right") +
  expand_limits(y = c(0, 1)) +
  scale_y_continuous(expand = c(1e-5, 1e-5),
                     labels = function(x) x * 100) +
  # scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 1)) +
  xlab("Fire size class (ha)") +
  ylab("Fire size class probability") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(b = 5, unit = "mm"),
        strip.background = element_blank(),
        strip.text = element_blank())
esc_veg


ggarrange2(esc_veg, esc_ndvi, ncol = 2)


## Topography _________________________

pd_topo <- rbind(
  expand.grid(
    elevation = seq(data_summ$elevation["hdi_lower_95"],
                    data_summ$elevation["hdi_upper_95"],
                    length.out = nseq),
    northing = 0,
    varying_var = "elevation"
  ),
  expand.grid(
    elevation = data_summ$elevation["mean"],
    northing = seq(-1, 1, length.out = nseq),
    varying_var = "northing"
  )
)

# Compute tfi
pd_topo$tfi <-
  fi_params$b_elev_ori * pd_topo$elevation +
  fi_params$b_north_ori * pd_topo$northing
pd_topo$tfi <- (pd_topo$tfi - fi_params$tfi_mean) / fi_params$tfi_sd

prob_topo <- ordinal_predict(pd_topo, "tfi", "b_tfi")

esc_elev <-
  ggplot(prob_topo[prob_topo$varying_var == "elevation", ],
         aes(elevation, mean, ymin = lower, ymax = upper,
             color = class_name, fill = class_name)) +
  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +
  expand_limits(y = 0:1) +
  scale_color_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                      name = "Fire size\nclass (ha)") +
  scale_fill_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                     name = "Fire size\nclass (ha)") +
  scale_y_continuous(expand = c(0.001, 0.001),
                     labels = function(x) x * 100) +
  xlab("Elevation (m a.s.l.)") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = bgcol, color = bgcol),
        strip.text = element_text(color = "white"),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(b = 5, unit = "mm"))
esc_elev

esc_north <-
  ggplot(prob_topo[prob_topo$varying_var == "northing", ],
         aes(northing, mean, ymin = lower, ymax = upper,
             color = class_name, fill = class_name)) +
  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +
  expand_limits(y = 0:1) +
  scale_color_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                      name = "Fire size\nclass (ha)") +
  scale_fill_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                     name = "Fire size\nclass (ha)") +
  scale_y_continuous(expand = c(0.001, 0.001),
                     labels = function(x) x * 100) +
  expand_limits(y = 0) +
  xlab("Northing") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(b = 5, l = 5, unit = "mm"))
esc_north


## Distance from roads _______________________

pd_roads <- data.frame(
  dist_roads = seq(0, quantile(pop$dist_roads_km, prob = 0.98), length.out = nseq)
)
pd_roads$drz <- (pd_roads$dist_roads - dr_mean) / dr_sd

prob_roads <- ordinal_predict(pd_roads, "drz", "b_drz")

esc_roads <-
  ggplot(prob_roads,
         aes(dist_roads, mean, ymin = lower, ymax = upper,
             color = class_name, fill = class_name)) +
  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +
  scale_color_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                      name = "Fire size\nclass (ha)") +
  scale_fill_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                     name = "Fire size\nclass (ha)") +
  scale_y_continuous(expand = c(0.001, 0.001),
                     labels = function(x) x * 100) +
  expand_limits(y = 0:1) +
  xlab("Distance from roads (km)") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(b = 5, unit = "mm"))
esc_roads

## Distance from humans _______________________

pd_humans <- data.frame(
  dist_humans = seq(0, quantile(pop$dist_humans_km, prob = 0.98), length.out = nseq)
)
pd_humans$dhz <- (pd_humans$dist_humans - dr_mean) / dr_sd

prob_humans <- ordinal_predict(pd_humans, "dhz", "b_dhz")

esc_humans <-
  ggplot(prob_humans,
         aes(dist_humans, mean, ymin = lower, ymax = upper,
             color = class_name, fill = class_name)) +
  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +
  scale_color_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                      name = "Fire size\nclass (ha)") +
  scale_fill_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                     name = "Fire size\nclass (ha)") +
  scale_y_continuous(expand = c(0.001, 0.001),
                     labels = function(x) x * 100) +
  expand_limits(y = 0:1) +
  xlab("Distance from humans (km)") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(b = 5, l = 5, unit = "mm"))
esc_humans

## FWI _________________________________________________

lshat <- as.matrix(scmod, "ls")
lsmean <- mean(lshat)
tseq <- 0:(nlags-1)
wfwi <- exp(-(tseq / lsmean) ^ 2) |> normalize()
fwi_points_z_dense <- as.numeric(fwi_points_z %*% wfwi)

pd_fwi <- data.frame(
  fwiz = seq(quantile(fwi_points_z_dense, prob = 0.01),
             quantile(fwi_points_z_dense, prob = 0.99), length.out = nseq)
)
pd_fwi$fwi <- pd_fwi$fwiz * fwi_sd + fwi_mean

prob_fwi <- ordinal_predict(pd_fwi, "fwiz", "b_fwi")

esc_fwi0 <-
  ggplot(prob_fwi,
         aes(fwi, mean, ymin = lower, ymax = upper,
             color = class_name, fill = class_name)) +
  geom_ribbon(color = NA, alpha = 0.4) +
  geom_line() +
  scale_color_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                      name = "Fire size\nclass (ha)") +
  scale_fill_viridis(discrete = T, begin = 0, end = 0.8, option = "B",
                     name = "Fire size\nclass (ha)") +
  scale_y_continuous(expand = c(0.001, 0.001),
                     labels = function(x) x * 100) +
  expand_limits(y = 0:1) +
  xlab("Fire Weather Index") +
  nice_theme() +
  theme(axis.title.y = element_blank(),
        strip.text = element_blank(),
        panel.spacing.x = unit(4, "mm"),
        plot.margin = margin(b = 5, unit = "mm"))
esc_fwi0

leg <- ggpubr::get_legend(esc_fwi0) |> ggpubr::as_ggplot()
# leg <- leg + theme(plot.margin = margin(1, 1, 1, 1, unit = "cm"))

esc_fwi <- esc_fwi0 + theme(legend.position = "none")

## Merge all plots ______________

ytit1 <- grid::textGrob("Fire size class probability (%)",
                        rot = 90, hjust = 0.5,
                        gp = grid::gpar(fontsize = 12, fontface = 'plain'))

escplot_veg <- ggarrange2(esc_veg,
                          esc_ndvi +
                            theme(plot.margin = margin(l = 5, unit = "mm")),
                          ncol = 2, left = ytit1, padding = unit(8, "mm"))

ggsave("ignition/figures/escape_veg.png",
       width = fig_width_max, height = fig_height_max * 0.8, units = "cm",
       plot = escplot_veg)

##
ytit1 <- grid::textGrob("Fire size class probability (%)",
                        rot = 90, hjust = 0.35,
                        gp = grid::gpar(fontsize = 12, fontface = 'plain'))

escplot_env <- ggarrange2(esc_elev, esc_north,
                          esc_roads, esc_humans,
                          esc_fwi, leg,
                          ncol = 2, left = ytit1, padding = unit(8, "mm"))

ggsave("ignition/figures/escape_environment.png",
       width = fig_width_max * 0.8, height = fig_height_max * 0.7, units = "cm",
       plot = escplot_env)


# FWI weights ____________________________________________________________

ls_hat <- as.numeric(as.matrix(scmod, "ls"))

tseq <- seq(0, nlags-1, by = 0.05)
np <- length(tseq)
ids_norm <- seq(1, np, by = 10)

wun <- matrix(NA, np, npost)
wn <- wun

for(i in 1:npost) {
  wun[, i] <- exp(-(tseq / ls_hat[i]) ^ 2)
  norm_const <- sum(wun[ids_norm, i])
  wn[, i] <- wun[, i] / norm_const # normalize only at discrete timepoints
}

plot(wun[, 1] ~ tseq, type = "l")
for(i in 1:100) lines(wun[, i+1] ~ tseq, col = rgb(0, 0, 0, 0.05))

plot(wn[, 1] ~ tseq, type = "l")
for(i in 1:100) lines(wn[, i+1] ~ tseq, col = rgb(0, 0, 0, 0.05))

# summarize
wsumm <- apply(wn, 1, mean_ci) |> t() |> as.data.frame()
wsumm$t <- tseq

ggplot(wsumm, aes(t, mean, ymin = lower, ymax = upper)) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  geom_point(data = wsumm[ids_norm, ], size = 2, alpha= 0.8) +
  scale_x_continuous(breaks = seq(0, 10, 2),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1),
                     expand = c(0.01, 0.01)) +
  ylab("Fortnight weight") +
  xlab("Fortnights before fire") +
  nice_theme()

nn <- file.path("ignition", "figures", "escape_fwi_weights.png")
ggsave(nn, width = 10, height = 9, units = "cm")



# TAREAS ------------------------------------------------------------------

# Mapear prob de ignition conditional y prob de escape en todo el parque.
# Mapear prob de propagación.
# Hacer algunos mapas de fuegos (spread).
# Simular fire size en función del FWI para el parque (spread).

# Poisson mixture ------------------------------------------------------

# Having 2 poisson, producing a number of ignitions from each source.
# When I take one ignition, which is the probability of having come from
# either component? It's the normalized lambda: lambda_focal / sum(lambdas)
it <- 5000
sample_prop <- numeric(it)
fac <- 1
l1 <- 100 * fac; l2 <- 70 * fac;
phi <- 1
prop_pop <- l1 / (l1 + l2)

for(i in 1:it) {
  n1 <- rpois(1, l1)
  n2 <- rpois(1, l2)
  sample_prop[i] <- n1 / (n1 + n2)
}
plot(density(sample_prop, from = 0, to = 1, n = 2 ^ 10), main = NA, xlab = NA)
abline(v = prop_pop, lty = 2)
# As ignitions are assumed to occur independently within each fortnight over
# the study area, it's sampling with replacement.

# And the negbin?
it <- 50000
sample_prop <- numeric(it)
fac <- 1000
l1 <- exp(-1) * fac; l2 <- exp(-3) * fac;
phi <- 1e6
prop_pop <- l1 / (l1 + l2)
for(i in 1:it) {
  n1 <- rnegbin(1, l1, phi)
  n2 <- rnegbin(1, l2, phi)
  sample_prop[i] <- n1 / (n1 + n2)
}
# remove zero total samples, which make NA
na_count <- sum(is.na(sample_prop))
tit <- paste("N =", it, "// notNA =", it-na_count,
             "// propNA =", round(na_count / it, 6))
sample_prop <- sample_prop[!is.na(sample_prop)]
plot(density(sample_prop, from = 0, to = 1, n = 2 ^ 10),
     main = tit, xlab = NA)
abline(v = prop_pop, lty = 2)
abline(v = mean(sample_prop), lty = 3, col = 2)
# Yes, on average, the proportion is the quotient of means. However, with small
# means, the distribution widens, and the mean is more variable, as a lower number
# of cases is used to compute the proportions.
