# This script prepares the landscapes where fires are going to be simulated.

# The first step is to extract the elevation layer from all images to use as
# input in windninja. These files will be saved in data/focal fires data/wind ninja files,
# as *fire_id*.tif.

# Once the wind layers are created, the landscape arrays will be made as follows:

# vegetation flammability index (vfi)
# topographic flammability index (tfi)
# elevation (m) (to compute directional slope effect)
# wind direction (degrees, from where the wind comes)
# wind speed (km / h, but probably scaled)

# The FWI will be held constant within each fire, taking its value at the
# centroid. This variable will model parameters at the fire-level, so it's
# included in the landscapes

library(terra)
library(tidyverse)
library(lubridate)
source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for land_cube
source(file.path("flammability indices",
                 "flammability_indices_functions.R"))
# functions to compute flammability indices

# average wind direction
# cicular mean for aspect
# from
# https://stackoverflow.com/questions/32404222/circular-mean-in-r
mean_circular_raw <- function (x) { # takes angle in radians
  sinr <- sum(sin(x))
  cosr <- sum(cos(x))
  circmean <- atan2(sinr, cosr)
  return(circmean)
}
mean_circular_deg <- function (x) { # takes angle in degrees
  # x <- c(180, 160, 170)
  conv <- 2 * pi / 360 # degrees to radians factor
  mm_raw <- mean_circular_raw(conv * x) / conv
  mm <- (mm_raw + 360) %% 360
  return(mm)
}


# get paths
gee_dir <- file.path("data", "focal fires data", "raw data from GEE")
windninja_dir <- "/home/ivan/windninja_cli_fire_spread_files" # path with no spaces


# Import vegetation class transforms --------------------------------------

dveg <- readxl::read_excel("/home/ivan/Insync/Mapa vegetación WWF - Lara et al. 1999/clases de vegetacion y equivalencias.xlsx",
                           sheet = "Sheet2")

# Make urban as forest, and get zero-indexing
dveg$cnum2[dveg$class1 == "Urban"] <- 1
dveg$class2[dveg$class1 == "Urban"] <- "Wet forest"
dveg$cnum3 <- dveg$cnum2
dveg$cnum3[dveg$cnum3 < 6] <- dveg$cnum3[dveg$cnum3 < 6] - 1
dveg$cnum3[dveg$cnum2 == 6] <- 99
# urban is taken as forest so its burn probability changes markedly with NDVI.

n_veg_types <- V <- 5



# Fires stuff -------------------------------------------------------------

# get fire_ids
fnames <- list.files(gee_dir)

fire_ids <- sapply(fnames, function(x) {
  # x <- "fire_data_raw_CoVentana.tif"
  id_raw_1 <- strsplit(x, "_|[.]")[[1]]
  id_raw_2 <- id_raw_1[-c(1:3, length(id_raw_1))]
  return(paste(id_raw_2, collapse = "_"))
}) %>% unname
fire_ids <- fire_ids[order(fire_ids)]

n_fires <- length(fire_ids)

# Load image stacks --------------------------------------------------------

raw_imgs <- vector("list", n_fires)
for(i in 1:n_fires) {
  raw_imgs[[i]] <- rast(file.path(gee_dir, fnames[i]))
}

# Export elevation layers for windninja ------------------------------------

for(i in 1:length(fnames)) {
  # i = 1
  print(i)
  r <- raw_imgs[[i]][["elevation"]]
  v <- values(r)
  if(anyNA(v)) {
    mval <- mean(r)
    r <- subst(r, NA, mval)
  }
  id_raw <- strsplit(fnames[i], "_")[[1]][-(1:3)]
  id <- paste(id_raw, collapse = "_")
  writeRaster(r, file.path(windninja_dir, id))
}

# Create and import wind layers -------------------------------------------

# take wind data for windninja from this file:
wind_data_allfires <- read.csv("data/climatic_data_by_fire_FWI-wind_corrected.csv")

# merge with fires with known ignition point:
points_raw <- vect(file.path("data", "ignition_points_checked.shp"))
ids_ig <- points_raw$Name %>% unique

# check these are the same as downloaded data
fires_gee_ids <- sapply(1:length(fnames), function(i) {
  # i = 1
  id_raw_1 <- strsplit(fnames[i], "[.]")[[1]][-2]
  id_raw_2 <- strsplit(id_raw_1, "_")[[1]][-(1:3)]
  id <- paste(id_raw_2, collapse = "_")
  return(id)
})
all(fires_gee_ids %in% ids_ig) # OK

# detect fire ids not present in wind data
ids_lack <- ids_ig[which(!(ids_ig %in% wind_data_allfires$fire_id))]
ids_add <- c("2015_47", "2011_19")
ids_ig <- c(ids_ig, ids_add)

# create wind table to use for windninja
wind_table <- wind_data_allfires[wind_data_allfires$fire_id %in% ids_ig, ]
wind_table <- wind_table[order(wind_table$fire_id), c("fire_id", "direction_use")]
wind_table <- wind_table[, c("fire_id", "direction_use")]

# repeat the separated fires
rows_sep <- which(wind_table$fire_id %in% ids_add)
wind_table <- rbind(wind_table, wind_table[rows_sep, ])
wind_table <- wind_table[order(wind_table$fire_id), ]
wind_table$fire_id_sep <- wind_table$fire_id
wind_table$fire_id_sep[wind_table$fire_id == "2015_47"] <- c("2015_47N", "2015_47S")
wind_table$fire_id_sep[wind_table$fire_id == "2011_19"] <- c("2011_19E", "2011_19W")
wind_table <- wind_table[, c("fire_id", "fire_id_sep", "direction_use")]

wind_table$elev_file <- paste(wind_table$fire_id_sep, ".tif", sep = "")
# write.csv(wind_table, file.path(
#   "data", "focal fires data", "wind ninja files", "zzz_wind_direction_table.csv"
# ))

# Run WindNinja using a fixed windspeed.

# Wind speed avg (from TerraClimate?)
# 4.78628962675556 m/s  # seems reasonable, way higher than ERA5 values
# 4.78628962675556 * 3.6
# use 4 m/s = 14.4 km / h

# mean wind direction: circular_mean(wind_table$)
mean_circular_deg(wind_table$direction_use) # 293

for(i in 1:n_fires) {
  print(wind_table$fire_id_sep[i])

  wind_direction <- wind_table$direction_use[i]
  elev_path <- file.path(windninja_dir, wind_table$elev_file[i])

  # create config file
  cat(file = file.path(windninja_dir, "spread_config_file.cfg"),
      sep = "\n",
# constants
"num_threads                = 15
initialization_method      = domainAverageInitialization
input_speed                = 4.0
input_speed_units          = mps
output_speed_units         = mps
input_wind_height          = 10.0
units_input_wind_height    = m
output_wind_height         = 10.0
units_output_wind_height   = m
vegetation                 = trees
mesh_resolution            = 90.0
units_mesh_resolution      = m
output_buffer_clipping     = 0.0
write_ascii_output         = true
ascii_out_resolution       = 30.0
units_ascii_out_resolution = m
momentum_flag              = false",
# variables: DEM and wind direction
paste("elevation_file             =", elev_path),
paste("input_direction            =", wind_direction)
  )

  # Run WindNinja
  system("WindNinja_cli --config_file=/home/ivan/windninja_cli_fire_spread_files/spread_config_file.cfg")
}

# delete configuration file
unlink("/home/ivan/windninja_cli_fire_spread_files/spread_config_file.cfg")

# Import wind speed and direction ------------------------------------------

wind_dir_files <- list.files(windninja_dir, pattern = "*ang.asc")
wind_vel_files <- list.files(windninja_dir, pattern = "*vel.asc")
# view(cbind(wind_dir_files, wind_vel_files, fire_ids)) # ok

# compute windspeed scale to standardize later
windspeed_vals <- vector("list", n_fires)

wind_imgs <- vector("list", n_fires)
for(i in 1:n_fires) {
  # i = 1
  dir <- rast(file.path(windninja_dir, wind_dir_files[i]))
  vel <- rast(file.path(windninja_dir, wind_vel_files[i]))
  wind <- c(dir, vel)
  names(wind) <- c("direction", "speed")
  wind_imgs[[i]] <- wind

  windspeed_vals[[i]] <- values(vel)
}

wind_sd <- sd(unlist(windspeed_vals)) # = 1.464333 # to standardize

# Import ignition points --------------------------------------------------

# ignition points for 2014_1 and 2008_5 (Ñorquinco and Lolog) were edited
# because these fires had a change of wind. The point to be used is not the
# real ignition point, but is coherent with the dominant wind direction.
# The original point id end in "_original", while the fire id was assigned to
# the edited points. These points are not included in the kml file.
points_raw <- vect(file.path("data", "ignition_points_checked.shp"))
points <- project(points_raw, raw_imgs[[1]])

# Import FWI data ---------------------------------------------------------

fwi_raw <- read.csv(file.path("data", "climatic_data_by_fire_fwi-fortnight-cumulative.csv"))
fwi_raw <- rename(fwi_raw, "fire_id_match" = "fire_id")

fwi_left <- data.frame(fire_id = fire_ids,
                       fire_id_match = fire_ids)

# rename repeated names (fires that were divided) to match
fwi_left$fire_id_match[grep("2015_47", fwi_left$fire_id)] <- "2015_47"
fwi_left$fire_id_match[grep("2011_19", fwi_left$fire_id)] <- "2011_19"

fwi_data <- left_join(fwi_left, fwi_raw, by = "fire_id_match")


# Import NDVI parameters estimates and other stuff ------------------------

fi_params <- readRDS(file.path("data", "flammability indices",
                               "flammability_indices.rds"))

# model to detrend ndvi
mdetrend <- readRDS(file.path("data", "flammability indices",
                              "ndvi_detrender_model.rds"))


# Northing importance function (slope weighting) --------------------------

# x is the slope steepness, in degrees
northing_weight <- function(x, a = -5, b = 0.35) {
  plogis(a + b * x)
}

# northing_weighted <- cos(aspect * pi / 180) *
#                      northing_weight(slope)

# Make landscapes ---------------------------------------------------------

# Landscape layers names in FireSpread package
# enum land_names {
#   veg,    # {0: forest, 1: shrubland, 2: grassland, 99: non-burnable}
#   ndvi    # scale(pi[v] * (ndvi - optim[v]) ^ 2, center = F)
#   north,  # scale(slope-weighted northing, center = F)
#   elev,   # scale(elevation)
#   wdir,
#   wspeed, # scale(wspeed, center = F)
# };

# Note that elevation is centred and scaled, but the other terms are only scaled.
# The slope term is not going to be scaled because it would require to
# compute more things during the simulation. However, its sd will be used
# to scale its beta to be compared with other betas.

# (using the same order here)
n_fi <- 2          # number of flammability indices
n_nd <- 1 + n_fi   # number of non-directional terms
n_layers <- n_nd + 3
land_names <- c("veg", "vfi", "tfi", "elev", "wdir", "wspeed")

# list with landscapes
lands <- vector(mode = "list", length = n_fires)
names(lands) <- fire_ids

# every fire will be a list
elem_names <- c("landscape", "ig_rowcol",
                "burned_layer", "burned_ids",
                "counts_veg", "counts_veg_available",
                "landscape_img",
                "fire_id", "fire_id_spread")
# landscape includes vegetation, but the spread function uses this separately.
# landscape_img is saved for plotting purposes, and contains the original
# layers.
# fire_id1 is the id for fire spread, which has some fires separated;
# fire_id

for(i in 1:n_fires) {
  # i = 40
  print(i)

  lands[[i]] <- vector(mode = "list", length = length(elem_names))
  names(lands[[i]]) <- elem_names

  lands[[i]]$fire_id_spread <- fire_ids[i]
  lands[[i]]$fire_id <- wind_table$fire_id[wind_table$fire_id_sep == fire_ids[i]]

  # create landscape image
  landscape_img <- raw_imgs[[i]][[c("veg", "ndvi_prev", "elevation", "slope",
                                    "aspect", "burned")]]
  names(landscape_img) <- c("veg", "ndvi_prev_dt", "elevation", "slope",
                            "aspect", "burned")

  ## Vegetation type
  veg_img <- raw_imgs[[i]]$veg
  veg_img <- subst(veg_img, dveg$cnum1, dveg$cnum3) # make NaN non-burnable
  veg_img <- subst(veg_img, NaN, 99) # make NaN non-burnable

  ## Vegetation flammability index
  ndvi_img <- raw_imgs[[i]]$ndvi_prev
  vfi_img <- raw_imgs[[i]]$ndvi_prev

  # get year to detrend NDVI
  fire_date <- fwi_data$date[fwi_data$fire_id == fire_ids[i]]
  fire_date <- as.Date(fire_date, format = "%Y-%m-%d")
  fire_month <- month(fire_date)
  fire_year <- ifelse(fire_month >= 7, year(fire_date) + 1, year(fire_date))

  # detrend NDVI, convert to 2022 equivalent
  if((fire_year - 1) < 2022) {
    ndvi_22_img <- raw_imgs[[i]]$ndvi_22

    dpred_ndvi <- data.frame(
      ndvi_dyn_logit = as.numeric(qlogis((values(ndvi_img) + 1) / 2)),
      ndvi01_22 = as.numeric((values(ndvi_22_img) + 1) / 2),
      year = fire_year - 1 # because NDVI is from the previous year
    )
    # anyNA(dpred_ndvi)
    dpred_ndvi$diff_logit <- mgcv::predict.gam(mdetrend, dpred_ndvi, se.fit = F)
    dpred_ndvi$ndvi_dt <- plogis(dpred_ndvi$ndvi_dyn_logit - dpred_ndvi$diff_logit) * 2 - 1
    values(ndvi_img) <- as.numeric(dpred_ndvi$ndvi_dt)
  }
  # check detrend:
  # plot(values(ndvi_img) ~ values(raw_imgs[[i]]$ndvi_prev))
  # summary(values(ndvi_img))
  # summary(values(raw_imgs[[i]]$ndvi_prev))

  # Compute VFI
  ndvi_dt <- values(ndvi_img)

  # put detrended ndvi in the landscape_img
  values(landscape_img$ndvi_prev_dt) <- ndvi_dt

  vfi <- ndvi_dt # just initialize
  veg_values <- values(veg_img)

  for(v in 1:5) {
    # compute quadratic term
    # v = 1
    id_fill <- veg_values == (v-1)

    vfi[id_fill] <-
      fi_params$a[v] +
      fi_params$b[v] *
      (ndvi_dt[id_fill] - fi_params$o[v]) ^ 2
  }

  # standardize vfi
  vfi_z <- (vfi - fi_params$vfi_mean) / fi_params$vfi_sd
  vfi_z[veg_values == 99] <- 0 # non burnable makes nothing
  values(vfi_img) <- as.numeric(vfi_z)
  names(vfi_img) <- "vfi"

  names(ndvi_img) <- "ndvi"
  # plot(hist(values(ndvi_img)))
  # plot(hist(ndvi_term))
  # plot(hist(ndvi_dt))
  # # check ndvi term:
  # plot(ndvi_term[veg_values != 99] ~
  #      ndvi_dt[veg_values != 99])

  ## Topographic flammability index
  vtopo <- values(landscape_img[[c("elevation", "slope", "aspect")]])
  northing <- cos(vtopo[, "aspect"] * pi / 180) *
              northing_weight(vtopo[, "slope"])

  tfi <-
    fi_params$b_elev_ori * vtopo[, "elevation"] +
    fi_params$b_north_ori * northing
  tfi_z <- (tfi - fi_params$tfi_mean) / fi_params$tfi_sd
  tfi_z[veg_values == 99] <- 0 # non burnable makes nothing
  tfi_img <- vfi_img
  values(tfi_img) <- tfi_z
  names(tfi_img) <- "tfi"

  ## Wind
  # project wind direction to match extent
  wind_local <- project(wind_imgs[[i]],
                        raw_imgs[[i]],
                        method = "cubicspline")
  names(wind_local) <- c("wdir", "wspeed")

  # turn wind direction to radians
  wind_local$wdir <- wind_local$wdir * (pi / 180)
  # scale windspeed
  wind_local$wspeed <- wind_local$wspeed / wind_sd

  # merge all data in a single raster
  rall <- c(veg_img, vfi_img, tfi_img, landscape_img$elevation, wind_local)

  # identify NA in the predictors to make those pixels unburnable
  vv <- values(rall)

  # identify burnable cells with NA values
  na_cells <- which(apply(vv, 1, anyNA))

  # make the na non-burnable
  vv[na_cells, "veg"] <- 99

  # edit values in raster
  values(rall) <- vv

  # replace NA with -9999 to avoid problemas with C++
  rall <- subst(rall, NA, -9999)
  rall <- subst(rall, NaN, -9999)

  # caution if there are a lot of NA in burnable area
  na_prop <- length(na_cells) / ncell(rall)
  if(na_prop > 0.02) warning(paste("Many NA in landscape, fire: ", fire_ids[i], ", i: ", i, sep = ""))

  # make array for c++
  land_arr <- land_cube(rall)
  # str(land_arr)

  # subset in landscape and other stuff
  lands[[i]]$landscape <- land_arr
  lands[[i]]$landscape_img <- landscape_img # save as terra raster with original
  # values of predictors, to make maps

  ## Fire-related elements
  burned_img <- landscape_img$burned
  lands[[i]]$burned_layer <- land_cube(burned_img)[, , 1]
  lands[[i]]$burned_layer[is.na(lands[[i]]$burned_layer)] <- 0

  # compute burned_ids (with 0-indexing!)
  burned_cells <- which(values(burned_img) == 1)
  burned_rowcols <- rowColFromCell(raw_imgs[[i]], burned_cells)
  colnames(burned_rowcols) <- c("row", "col")
  lands[[i]]$burned_ids <- t(burned_rowcols) - 1 # for 0-indexing!

  # pixels by vegetation type in available and burned
  veg_vec <- as.vector(land_arr[, , "veg"])                       # - 1 # 0 is non burnable now.
  burned_vec <- as.vector(lands[[i]]$burned_layer)
  counts_veg <- integer(n_veg_types)
  counts_veg_available <- integer(n_veg_types)
  for(k in 1:n_veg_types) {
    available <- as.numeric(veg_vec == (k-1))
    counts_veg[k] <- sum(available * burned_vec)
    counts_veg_available[k] <- sum(available)
  }

  lands[[i]]$counts_veg <- counts_veg
  lands[[i]]$counts_veg_available <- counts_veg_available

  # get ignition points
  p_local <- points[points$Name == fire_ids[i]]

  # get coordinates and row_col
  cc <- crds(p_local)
  ig_rowcol <- rbind(rowFromY(raw_imgs[[i]], cc[, "y"]),
                     colFromX(raw_imgs[[i]], cc[, "x"]))
  row.names(ig_rowcol) <- c("row", "col")

  if(anyNA(ig_rowcol)) {
    stop(paste("Ignition point out of range,", "fire_id", fire_ids[i]))
  }

  lands[[i]]$ig_rowcol <- ig_rowcol - 1 # 0-indexing!!!
}

# Warning messages:
# 1: Many NA in landscape, fire: 2015_47S, i: 43
# 2: Many NA in landscape, fire: 2015_50, i: 44

# # Checks:
# f <- 44
# v <- 6
# l <- land_names[v]
# xxx <- lands[[f]]$landscape[, , l] %>% as.numeric
# sum(xxx != -9999) / length(xxx)
# hist(xxx[xxx != -9999])
# range(xxx[xxx != -9999])
# summary(xxx[xxx != -9999])

# Check all ignition points fall in burned and burnable cells,
# and all burned cells are burnable.
ccc <- numeric(n_fires)
# bbb <- numeric(n_fires)
for(i in 1:n_fires) {
  # i = 1
  ig <- lands[[i]]$ig_rowcol + 1 # because of 0-indexing

  point_checks <- sapply(1:ncol(ig), function(c) {
    # c = 1
    cbind(as.numeric(lands[[i]]$landscape[ig[1, c], ig[2, c], "veg"] < 99), # 1 is non-burnable
          lands[[i]]$burned_layer[ig[1, c], ig[2, c]])
  }) %>% colSums %>% unique()

  ccc[i] <- point_checks

  if(point_checks != 2) {
    stop(paste("Ignition point problems,", "fire_id:", fire_ids[i], "i:", i))
  }
}
all(ccc == 2) # OK

object.size(lands) / 1e6 # 2627.7 Mb

# sapply(lands, class) %>% unique # OK

# Save landscapes ---------------------------------------------------------

# saveRDS(lands, file.path("data", "landscapes_ig-known_non-steppe.rds"))
# lands <- readRDS(file.path("data", "landscapes_ig-known_non-steppe.rds"))

# save each landscape separately
for(i in 1:length(lands)) {
  print(i)
  saveRDS(lands[[i]], file.path("data", "focal fires data",
                                "landscapes",
                                paste(lands[[i]]$fire_id_spread, ".rds", sep = "")))
}


# Evaluate relative abundance of veg types in landscapes ------------------

dveg

# table with counts of veg types by fire
nv <- 7
veg_counts <- matrix(0, n_fires, nv)
colnames(veg_counts) <- dveg$class1[1:nv]

for(f in 1:n_fires) {
  # f = 2
  veg_img <- raw_imgs[[f]]$veg
  vals <- values(veg_img) %>% as.vector
  tt <- table(vals)
  names(tt)

  class_present <- as.numeric(names(tt))
  class_eval <- class_present[class_present %in% (1:nv)]
  ttsub <- tt[as.character(class_eval)]

  for(c in 1:length(class_eval)) {
    veg_counts[f, class_eval[c]] <- ttsub[c]
  }
}

veg_props <- apply(veg_counts, 1, function(x) x / sum(x)) %>% t

par(mfrow = c(2, 4))
for(v in 1:nv) {
  hist(veg_props[, v], breaks = seq(0, 1, by = 0.1),
       main = colnames(veg_props)[v],
       xlim = c(0, 1), xlab = "Relative abundance\nin landscape")
}
par(mfrow = c(1, 1))

# > 5 % in how many fires?
apply(veg_props, 2, function(x) sum(x < 0.05)) %>% t %>% t

# Condense araucaria and cypres as dry forest, and plantation as shrubland
veg_props_sub <- cbind(
  veg_props[, c("Wet forest", "Subalpine forest")],
  "Dry forest" = rowSums(veg_props[, c("Araucaria forest", "Cipres forest")]),
  "Shrubland" = rowSums(veg_props[, c("Shrubland", "Plantation")]),
  veg_props[, "Grassland"]
)

# Number of veg types with more than 5 % cover by landscape
veg_num <- apply(veg_props_sub, 1, function(x) sum(x >= 0.05)) %>% t %>% t
table(veg_num)

# Check landscape size of those with five
fire_size <- rowSums(veg_counts)
fire_size_rel <- fire_size / max(fire_size)
barplot(fire_size_rel)

fire_size_rel[veg_num == 5]
# most of them are small;



# PNNH landscape ----------------------------------------------------------

# we use the smaller landscape, which besides being smaller than the large
# buffer, has NDVI from year 2021... That avoids the low values around the
# steffen-martin fire occurred in 2022.

# load image
pnnh_rast <- rast(file.path("data", "pnnh_images",
                            "pnnh_data_spread_buffered_30m.tif"))

# export elevation for WindNinja
r <- pnnh_rast[["elevation"]]
v <- values(r)
if(anyNA(v)) {
  mval <- mean(v, na.rm = T)
  r <- subst(r, NA, mval)
}
writeRaster(r, file.path("data", "pnnh_images",
                         "pnnh_data_spread_elevation_30m.tif"))


## Make wind layers

# dominant wind direction (°)
elev_path <- "/home/ivan/Insync/Fire spread modelling/fire_spread/data/pnnh_images/pnnh_data_spread_elevation_30m.tif"
# create config file
cat(file = file.path(windninja_dir, "spread_config_file_pnnh.cfg"),
    sep = "\n",
# constants
"num_threads               = 6
initialization_method      = domainAverageInitialization
input_speed                = 4.0
input_speed_units          = mps
output_speed_units         = mps
input_wind_height          = 10.0
units_input_wind_height    = m
output_wind_height         = 10.0
units_output_wind_height   = m
vegetation                 = trees
mesh_resolution            = 120.0
units_mesh_resolution      = m
output_buffer_clipping     = 0.0
write_ascii_output         = true
ascii_out_resolution       = 30.0
units_ascii_out_resolution = m
momentum_flag              = false
input_direction            = 293",
paste("elevation_file      =", elev_path)
)
# The RAM was not enought to use mesh_resolution = 90.0

# Run WindNinja
# system("WindNinja_cli --config_file=/home/ivan/windninja_cli_fire_spread_files/spread_config_file_pnnh.cfg")

## Load wind images and project

wdir_rast <- rast(file.path("data", "pnnh_images",
                            "pnnh_data_spread_elevation_30m_293_4_30m_ang.asc"))
wspeed_rast <- rast(file.path("data", "pnnh_images",
                              "pnnh_data_spread_elevation_30m_293_4_30m_vel.asc"))

wind_rast <- c(wdir_rast, wspeed_rast)
names(wind_rast) <- c("direction", "speed")


# Landscape layers names in FireSpread package
# enum land_names {
#   veg,    # {0: forest, 1: shrubland, 2: grassland, 99: non-burnable}
#   ndvi    # scale(pi[v] * (ndvi - optim[v]) ^ 2, center = F)
#   north,  # scale(slope-weighted northing, center = F)
#   elev,   # scale(elevation)
#   wdir,
#   wspeed, # scale(wspeed, center = F)
# };

# Note that elevation is centred and scaled, but the other terms are only scaled.
# The slope term is not going to be scaled because it would require to
# compute more things during the simulation. However, its sd will be used
# to scale its beta to be compared with other betas.

# (using the same order here)
n_fi <- 2          # number of flammability indices
n_nd <- 1 + n_fi   # number of non-directional terms
n_layers <- n_nd + 3
land_names <- c("veg", "vfi", "tfi", "elev", "wdir", "wspeed")

# landscape includes vegetation, but the spread function uses this separately.

## Vegetation type
veg_img <- pnnh_rast$veg
veg_img <- subst(veg_img, dveg$cnum1, dveg$cnum3) # make NaN non-burnable
veg_img <- subst(veg_img, NaN, 99) # make NaN non-burnable

## Vegetation flammability index
veg5 <- subst(pnnh_rast[["veg"]], dveg$cnum1, dveg$cnum2) # cnum2 has 1:5
vfi_img <- pnnh_rast$veg # placeholder
values(vfi_img) <- vfi_calc(values(veg5), values(pnnh_rast$ndvi))
names(vfi_img) <- "vfi"

## Topographic flammability index
vtopo <- values(pnnh_rast[[c("elevation", "slope", "aspect")]])
tfi_img <- vfi_img # placeholder
values(tfi_img) <- tfi_calc(vtopo[, "elevation"], vtopo[, "aspect"],
                            vtopo[, "slope"])
names(tfi_img) <- "tfi"

## Wind
# project wind direction to match extent
wind_local <- project(wind_rast,
                      pnnh_rast,
                      method = "cubicspline")
names(wind_local) <- c("wdir", "wspeed")

# turn wind direction to radians
wind_local$wdir <- wind_local$wdir * (pi / 180)
# scale windspeed
wind_sd <- 1.464333
wind_local$wspeed <- wind_local$wspeed / wind_sd

# merge all data in a single raster
rall <- c(veg_img, vfi_img, tfi_img, pnnh_rast$elevation, wind_local)

# identify NA in the predictors to make those pixels unburnable
vv <- values(rall)

# identify burnable cells with NA values
na_cells <- which(apply(vv, 1, anyNA))

# make the na non-burnable
vv[na_cells, "veg"] <- 99

# edit values in raster
values(rall) <- vv

# replace NA with -9999 to avoid problemas with C++
rall <- subst(rall, NA, -9999)
rall <- subst(rall, NaN, -9999)

# caution if there are a lot of NA in burnable area
(na_prop <- length(na_cells) / ncell(rall)) * 100

# make array for c++
land_arr <- land_cube(rall)
# str(land_arr)

format(object.size(land_arr), units = "Mb")
saveRDS(land_arr, file.path("data", "pnnh_images", "pnnh_spread_landscape.rds"))

# check
apply(land_arr, 3, function(x) summary(as.vector(x)))
# OK


# PNNH landscape with unburnable urban ------------------------------------

# load image
pnnh_rast <- rast(file.path("data", "pnnh_images",
                            "pnnh_data_spread_buffered_30m.tif"))

## Load wind images and project

wdir_rast <- rast(file.path("data", "pnnh_images",
                            "pnnh_data_spread_elevation_30m_293_4_30m_ang.asc"))
wspeed_rast <- rast(file.path("data", "pnnh_images",
                              "pnnh_data_spread_elevation_30m_293_4_30m_vel.asc"))

wind_rast <- c(wdir_rast, wspeed_rast)
names(wind_rast) <- c("direction", "speed")


# Landscape layers names in FireSpread package
# enum land_names {
#   veg,    # {0: forest, 1: shrubland, 2: grassland, 99: non-burnable}
#   ndvi    # scale(pi[v] * (ndvi - optim[v]) ^ 2, center = F)
#   north,  # scale(slope-weighted northing, center = F)
#   elev,   # scale(elevation)
#   wdir,
#   wspeed, # scale(wspeed, center = F)
# };

# Note that elevation is centred and scaled, but the other terms are only scaled.
# The slope term is not going to be scaled because it would require to
# compute more things during the simulation. However, its sd will be used
# to scale its beta to be compared with other betas.

# (using the same order here)
n_fi <- 2          # number of flammability indices
n_nd <- 1 + n_fi   # number of non-directional terms
n_layers <- n_nd + 3
land_names <- c("veg", "vfi", "tfi", "elev", "wdir", "wspeed")

# landscape includes vegetation, but the spread function uses this separately.

## Vegetation type

# turn urban into unburnable ______________________________________
dveg$cnum4 <- dveg$cnum3
dveg$cnum4[dveg$class1 == "Urban"] <- 99

# the cnum5 is used for vfi
dveg$cnum5 <- dveg$cnum2
dveg$cnum5[dveg$class1 == "Urban"] <- 6

veg_img <- pnnh_rast$veg
veg_img <- subst(veg_img, dveg$cnum1, dveg$cnum4) 
veg_img <- subst(veg_img, NaN, 99) # make NaN non-burnable

## Vegetation flammability index
veg5 <- subst(pnnh_rast[["veg"]], dveg$cnum1, dveg$cnum5) # cnum2=5 has 1:5
vfi_img <- pnnh_rast$veg # placeholder
values(vfi_img) <- vfi_calc(values(veg5), values(pnnh_rast$ndvi))
names(vfi_img) <- "vfi"

## Topographic flammability index
vtopo <- values(pnnh_rast[[c("elevation", "slope", "aspect")]])
tfi_img <- vfi_img # placeholder
values(tfi_img) <- tfi_calc(vtopo[, "elevation"], vtopo[, "aspect"],
                            vtopo[, "slope"])
names(tfi_img) <- "tfi"

## Wind
# project wind direction to match extent
wind_local <- project(wind_rast,
                      pnnh_rast,
                      method = "cubicspline")
names(wind_local) <- c("wdir", "wspeed")

# turn wind direction to radians
wind_local$wdir <- wind_local$wdir * (pi / 180)
# scale windspeed
wind_sd <- 1.464333
wind_local$wspeed <- wind_local$wspeed / wind_sd

# merge all data in a single raster
rall <- c(veg_img, vfi_img, tfi_img, pnnh_rast$elevation, wind_local)

# identify NA in the predictors to make those pixels unburnable
vv <- values(rall)

# identify burnable cells with NA values
na_cells <- which(apply(vv, 1, anyNA))

# make the na non-burnable
vv[na_cells, "veg"] <- 99

# edit values in raster
values(rall) <- vv

# replace NA with -9999 to avoid problemas with C++
rall <- subst(rall, NA, -9999)
rall <- subst(rall, NaN, -9999)

# caution if there are a lot of NA in burnable area
(na_prop <- length(na_cells) / ncell(rall)) * 100

# make array for c++
land_arr <- land_cube(rall)
# str(land_arr)

format(object.size(land_arr), units = "Mb")
saveRDS(land_arr, file.path("data", "pnnh_images", "pnnh_spread_landscape_urban-nonburnable.rds"))

# check
apply(land_arr, 3, function(x) summary(as.vector(x)))
# OK