# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(rstan)
library(lubridate)
library(terra)
library(bayesplot)

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

# Later, change this: add 10 days before so it can be cumulative
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

# Ignition rate model (only time) ----------------------------------------

# no cumulative FWI effect, treating unknown as human.
smodel <- stan_model("ignition_model.stan")

# How to compute the likelihood?
condition_like <- rep(1, nfort)           # evaluate both human and lightning-caused, known
condition_like[ig_counts[, 3] > 0] <- 2   # consider unknown cause
condition_like[time_data$cause_both == 0] <- 3 # known, but only lightning

# compute average ignition rate to set a non-influential prior on the median
ig_transient <- ig_counts[, 1:2]
ig_transient[, 1] <- ig_transient[, 1] + ig_counts[, 3]

sdata <- list(
  nfort = nfort, npix = npix_fwi, nlag = nlags,
  ig_counts = ig_counts,
  condition = condition_like,
  fwi = aperm(fwi_arr, c(3, 1, 2)),
  pix_weights = pix_weights,
  nfort_both = sum(condition_like < 3),

  prior_a_sd = 5,
  prior_a_mu = log(mean(ig_transient)),
  prior_b_sd = 5,
  prior_phi_sd = 50,
  prior_ls_sd = nlags * 0.75
)

igmod <- sampling(
  smodel, data = sdata, seed = 85963, refresh = 100,
  # cores = 1, chains = 1, iter = 10,
  cores = 6, chains = 6, iter = 2000, warmup = 1000,
  pars = c("a", "b", "phi", "ls", "humanp")
)

# glimpse at the posteriors
mcmc_dens(igmod, pars = c("a[1]", "a[2]"),
          facet_args = list(scales = "fixed", nrow = 2))
mcmc_dens(igmod, pars = c("b[1]", "b[2]"),
          facet_args = list(scales = "fixed", nrow = 2))
mcmc_dens(igmod, pars = c("phi[1]", "phi[2]"),
          facet_args = list(scales = "fixed", nrow = 2)) + xlim(0, 8)
mcmc_dens(igmod, pars = "ls") # biutiful, poca corr temporal

# much higher variation in lightning, with lower mean and higher FWI eff

lshat <- as.matrix(igmod)[, "ls"]
curve(exp(-((x / mean(lshat))^2)), to = 5)
for(i in sample(1:length(lshat), 200, replace = F)) {
  curve(exp(-((x / lshat[i])^2)), add = T, col = rgb(0, 0, 0, 0.05))
}


# Tidy data for spatial analysis ------------------------------------------

# spatial predictor, sample from the population of non-burned pixels
nland <- 5000
xland <- rnorm(nland)

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

# write database to be uploaded to Earth engine, so we can sample the
# spatial covariates.
writeVector(project(igboth, "EPSG:4326"),
            file.path(igdata_dir, "ignition_points_pnnh_bari-kitzberger.shp"),
            overwrite = T)
writeVector(project(pnnh, "EPSG:4326"),
            file.path(igdata_dir, "pnnh.shp"),
            overwrite = T) # used for the population sample

# Spatio-temporal model (simulated predictor) ----------------------------

smodel2 <- stan_model("ignition_model2.stan")

# How to compute the likelihood for ignition counts?
condition_like <- rep(1, nfort)           # evaluate both human and lightning-caused, known
condition_like[ig_counts[, 3] > 0] <- 2   # consider unknown cause
condition_like[time_data$cause_both == 0] <- 3 # known, but only lightning

# compute average ignition rate to set a non-influential prior on the median
ig_transient <- ig_counts[, 1:2]
ig_transient[, 1] <- ig_transient[, 1] + ig_counts[, 3]

# spatial predictor, sample from the population of non-burned pixels
nland <- 5000
xland <- rnorm(nland)

# merge all ignition data
igboth <- rbind(igpn_sub, igl_sub)
npoint <- nrow(igboth)
igboth$cause_num <- as.numeric(factor(igboth$cause,
                               levels = c("human", "lightning", "unknown")))
nbycause <- as.numeric(table(igboth$cause_num))
# simulate the only predictor for now
igboth$x <- rnorm(npoint)
igboth$x[igboth$cause_num != 2] <- rnorm(sum(nbycause[c(1, 3)]), 3) # unk treated as human
# boxplot(x ~ cause, igboth)

# identify the corresponding row in the time_data for each ignition point, to
# match the expected rate of every ignition type.
igboth$row_rate <- NA
for(i in 1:npoint) {
  fortt <- igboth$fort[i]
  row <- which(time_data$fort == fortt)
  igboth$row_rate[i] <- row
}

# write database to be uploadede to Earth engine




sdata2 <- list(
  ## ignition rate model
  nfort = nfort, npix = npix_fwi, nlag = nlags,
  ig_counts = ig_counts,
  condition = condition_like,
  fwi = aperm(fwi_arr, c(3, 1, 2)),
  pix_weights = pix_weights,
  nfort_both = sum(condition_like < 3),

  prior_a_sd = 5,
  prior_a_mu = log(mean(ig_transient)),
  prior_b_sd = 5,
  prior_phi_sd = 50,
  prior_ls_sd = nlags * 0.75,

  ## ignition location model
  npoint = npoint, nland = nland,
  xland = xland, xpoints = igboth$x,
  cause = igboth$cause_num,
  row_rate = igboth$row_rate,

  prior_c_sd = 10
)

igmod2 <- sampling(
  smodel2, data = sdata2, seed = 85963, refresh = 100,
  # cores = 1, chains = 1, iter = 10,
  cores = 6, chains = 6, iter = 2000, warmup = 1000,
  pars = c("a", "b", "phi", "ls",
           "c", "humanp")
)
# 1.5 min

# glimpse at the posteriors
mcmc_dens(igmod2, pars = c("a[1]", "a[2]"),
          facet_args = list(scales = "fixed", nrow = 2))
mcmc_dens(igmod2, pars = c("b[1]", "b[2]"),
          facet_args = list(scales = "fixed", nrow = 2))
mcmc_dens(igmod2, pars = c("phi[1]", "phi[2]"),
          facet_args = list(scales = "fixed", nrow = 2)) + xlim(0, 8)
mcmc_dens(igmod2, pars = "ls") # biutiful, poca corr temporal

mcmc_dens(igmod2, pars = c("c[1]", "c[2]"),
          facet_args = list(scales = "fixed", nrow = 2))

mcmc_dens(igmod, pars = "humanp") + xlim(0, 1)
mcmc_dens(igmod2, pars = "humanp") + xlim(0, 1)


# observed proportion if all unknown were human-cuased, as in the simulation
rows_compare <- igboth$date <= max(time_human)
tt <- table(igboth$cause[rows_compare])
humanp_obs <- tt / sum(tt)
humanp_obs <- 1 - humanp_obs[2]

plot(density(as.matrix(igmod, "humanp"), from = 0, to = 1, n = 2 ^ 10),
     main = NA, xlab = "human proportion", ylim = c(0, 18))
lines(density(as.matrix(igmod2, "humanp"), from = 0, to = 1, n = 2 ^ 10),
      col = "red")
abline(v = humanp_obs, lty = 2)
# perfect. The difference between human and lighting proportions increases
# with their difference in x.

# Thinking of 2-poisson mixture.  ----------------------------------------

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
