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

# get paths
gee_dir <- file.path("data", "focal fires data", "raw data from GEE")
windninja_dir <- "/home/ivan/windninja_cli_fire_spread_files" # path with no spaces

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

wind_sd <- mean(unlist(lapply(windspeed_vals, sd))) # ~ 1.14 # to standardize

# Import ignition points --------------------------------------------------

# ignition points for 2014_1 and 2008_5 (Ñorquinco and Lolog) were edited
# because these fires had a change of wind. The point to be used is not the
# real ignition point, but is coherent with the dominant wind direction.
# The original point id end in "_original", while the fire id was assigned to
# the edited points. These points are not included in the kml file.
points_raw <- vect(file.path("data", "ignition_points_checked.shp"))
points <- project(points_raw, raw_imgs[[1]])

# Import FWI data ---------------------------------------------------------

fwi_raw <- read.csv(file.path("data", "climatic_data_by_fire_fwi-day-cumulative.csv"))
fwi_raw <- rename(fwi_raw, "fire_id_match" = "fire_id")

fwi_left <- data.frame(fire_id = fire_ids,
                       fire_id_match = fire_ids)

# rename repeated names (fires that were divided) to match
fwi_left$fire_id_match[grep("2015_47", fwi_left$fire_id)] <- "2015_47"
fwi_left$fire_id_match[grep("2011_19", fwi_left$fire_id)] <- "2011_19"

fwi_data <- left_join(fwi_left, fwi_raw, by = "fire_id_match")


# Import vegetation class transforms --------------------------------------

dveg <- readxl::read_excel("/home/ivan/Insync/Mapa vegetación WWF - Lara et al. 1999/clases de vegetacion y equivalencias.xlsx",
                           sheet = "Sheet2")

# recode vegetation type in 5 classes

n_veg_types <- 5
veg_names <- c("wet", "subalpine", "dry", "shrubland", "grassland")

dveg$class5 <- plyr::revalue(dveg$class1,
                             replace = c(
                               "Araucaria forest" = "Dry forest",
                               "Cipres forest" = "Dry forest",
                               "Plantation" = "Shrubland",
                               "Urban" = "Grassland",
                               "Water" = "Non burnable",
                               "High andean" = "Non burnable",
                               "Ice and snow" = "Non burnable"
                             ))

dveg$class5 <- factor(dveg$class5, levels = c(
  "Wet forest", "Subalpine forest", "Dry forest", "Shrubland", "Grassland",
  "Non burnable"
))

dveg$cnum5 <- as.numeric(dveg$class5) - 1
dveg$cnum5[dveg$class5 == "Non burnable"] <- 99


# Make landscapes ---------------------------------------------------------

# The landscape array will have the vegetation layer, and then terraion:
#   elev = elevation, wdir = wind direction, wspeed = windspeed

# (using the same order here)
land_names <- c("veg", "elev", "wdir", "wspeed")

# list with landscapes
lands <- vector(mode = "list", length = n_fires)
names(lands) <- fire_ids

# every fire will be a list
elem_names <- c("landscape", "ig_rowcol",
                "burned_layer", "burned_ids",
                "landscape_img",
                "fire_id", "fire_id_spread",
                "cells_by_veg")
# landscape includes vegetation, but the spread function uses this separately.
# landscape_img is saved for plotting purposes, and contains the original
# layers.
# fire_id1 is the id for fire spread, which has some fires separated;
# fire_id

for(i in 1:n_fires) {
  # i = 1
  print(i)

  lands[[i]] <- vector(mode = "list", length = length(elem_names))
  names(lands[[i]]) <- elem_names

  lands[[i]]$fire_id_spread <- fire_ids[i]
  lands[[i]]$fire_id <- wind_table$fire_id[wind_table$fire_id_sep == fire_ids[i]]

  # create landscape image
  landscape_img <- raw_imgs[[i]][[c("veg", "elevation", "burned")]]

  ## Vegetation type
  veg_img <- raw_imgs[[i]]$veg
  veg_img <- subst(veg_img, NaN, 11) # make NaN Ice and Snow (non-burnable)
  # duplicate to reclassify
  veg_img_reclass <- veg_img
  values(veg_img_reclass) <- NA
  unique_vegs_raw <- unique(veg_img)[, 1]

  # cnum1 is the number label in the image.
  for(v in unique_vegs_raw) {
    veg_img_reclass[veg_img == v] <- dveg$cnum5[dveg$cnum1 == v]
  }

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

  ## elevation (standardized)
  elev_img <- raw_imgs[[i]]$elevation
  names(elev_img) <- "elev"

  # merge all data in a single raster
  rall <- c(veg_img_reclass, elev_img, wind_local)

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

  # burned image
  burned_img <- raw_imgs[[i]]$burned
  lands[[i]]$burned_layer <- land_cube(burned_img)[, , 1]
  lands[[i]]$burned_layer[is.na(lands[[i]]$burned_layer)] <- 0

  # compute burned_ids (with 0-indexing!)
  burned_cells <- which(values(burned_img) == 1)
  burned_rowcols <- rowColFromCell(raw_imgs[[i]], burned_cells)
  colnames(burned_rowcols) <- c("row", "col")
  lands[[i]]$burned_ids <- t(burned_rowcols) - 1 # for 0-indexing!

  # cells by vegetation type in available and burned
  veg_vec <- as.vector(land_arr[, , "veg"])
  cells_by_veg <- integer(n_veg_types)
  names(cells_by_veg) <- veg_names
  for(k in 1:n_veg_types) {
    cells_by_veg[k] <- sum(veg_vec == (k-1)) # zero-indexing
  }

  lands[[i]]$cells_by_veg <- cells_by_veg

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
# f <- 43
# v <- 4
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

object.size(lands) / 1e6 # 1884.2 Mb

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

# dveg
#
# # table with counts of veg types by fire
# nv <- 7
# veg_counts <- matrix(0, n_fires, nv)
# colnames(veg_counts) <- dveg$class1[1:nv]
#
# for(f in 1:n_fires) {
#   # f = 2
#   veg_img <- raw_imgs[[f]]$veg
#   vals <- values(veg_img) %>% as.vector
#   tt <- table(vals)
#   names(tt)
#
#   class_present <- as.numeric(names(tt))
#   class_eval <- class_present[class_present %in% (1:nv)]
#   ttsub <- tt[as.character(class_eval)]
#
#   for(c in 1:length(class_eval)) {
#     veg_counts[f, class_eval[c]] <- ttsub[c]
#   }
# }
#
# veg_props <- apply(veg_counts, 1, function(x) x / sum(x)) %>% t
#
# par(mfrow = c(2, 4))
# for(v in 1:nv) {
#   hist(veg_props[, v], breaks = seq(0, 1, by = 0.1),
#        main = colnames(veg_props)[v],
#        xlim = c(0, 1), xlab = "Relative abundance\nin landscape")
# }
# par(mfrow = c(1, 1))
#
# # > 5 % in how many fires?
# apply(veg_props, 2, function(x) sum(x < 0.05)) %>% t %>% t
#
# # Condense araucaria and cypres as dry forest, and plantation as shrubland
# veg_props_sub <- cbind(
#   veg_props[, c("Wet forest", "Subalpine forest")],
#   "Dry forest" = rowSums(veg_props[, c("Araucaria forest", "Cipres forest")]),
#   "Shrubland" = rowSums(veg_props[, c("Shrubland", "Plantation")]),
#   veg_props[, "Grassland"]
# )

veg_counts <- do.call("rbind", lapply(lands, function(x) x[["cells_by_veg"]]))
veg_props <- apply(veg_counts, 1, function(x) x / sum(x))  %>% t

# Number of veg types with more than thres % cover by landscape
thres <- 0.05
veg_num <- apply(veg_props, 1, function(x) sum(x >= 0.05)) %>% t %>% t
table(veg_num)

# Check landscape size of those with five
fire_size <- rowSums(veg_counts)
fire_size_rel <- fire_size / max(fire_size)
barplot(fire_size_rel)

fire_size_rel[veg_num == 5]
# most of them are small; Try GAMs with those, to see if the reduction is needed.
# In case of needing reduction, merge dry with shrubland

# 1999_25j           2002_33           2012_57           2012_58
# 0.01244363        0.03855622        0.08417732        0.03914146

# 2015_40           2015_42 2021_2146405150_W          2021_865
# 0.05576505        0.09326619        0.28873183        0.52477964

# CoVentana
# 0.00453907

plot(lands[["1999_25j"]]$landscape_img[["veg"]])
# Este es el smallest. Los otros son medio grandes como para hacer pruebitas

fire_size_rel[veg_num == 4]
# 1999_27j_N        1999_27j_S           1999_28           2002_36
# 0.0159860951      0.1628276201      0.0092614464      0.0933731752
# 2002_7           2004_23   2008_1379405717            2008_5
# 0.0052959697      0.0052064508      0.1922772424      0.0972374906
# 2009_2007583421           2009_26           2009_28           2011_18
# 0.0858063770      0.0506278089      0.0116577566      0.0172049482
# 2013_30            2014_1           2015_11           2015_46
# 0.0017321633      0.0918380310      0.0063268591      0.0162512881
# 2015_50           2015_53          2016_47j 2021_2146405150_E
# 1.0000000000      0.0341996823      0.0421235210      0.2582066896
# 2021_936 2022_2125136700_r
# 0.0003867937      0.3371249039


fire_size_rel[veg_num == 2]
# 2005_9     2006_20E      2013_12      2013_32      2015_41     2021_911
# 0.0004063275 0.0018928315 0.0003232768 0.0004305183 0.0017972326 0.0003874406
# CoCampana
# 0.0003076239