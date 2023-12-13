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
# centroid. This variable will model a fire-level random intercept, so it's not
# considered by the fire-spread functions.

library(terra)
library(tidyverse)
source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for land_cube

# get paths
gee_dir <- file.path("data", "focal fires data", "raw data from GEE")
windninja_dir <- file.path("data", "focal fires data", "wind ninja files")

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

# Trim a few landscapes with NA at their edges ----------------------------

# This happens for fires that were on the edge of the topography layer,
# the northernmost, easternmost and westernmost.

# # Cholila
# fire_ori <- rast(file.path(gee_dir, "fire_data_raw_2015_50.tif")) 
# fire_trim <- trim(fire_ori[["elev"]], value = NA) # no hace nada
# velev <- values(fire_ori[["elev"]])
# na_elev_ids <- which(is.na(velev))
# length(na_elev_ids) / ncell(fire_ori)
# cols_na <- colFromCell(fire_ori, na_elev_ids) %>% unique
# fire_trim <- fire_ori[, -cols_na, drop = F]
# fire_trim
# writeRaster(fire_trim, file.path(gee_dir, "fire_data_raw_2015_50-trimmed.tif"))
# 
# # San Ramón
# fire_ori <- rast(file.path(gee_dir, "fire_data_raw_1999_2140469994_r.tif"))
# fire_ori <- rast("/home/ivan/Insync/Fire spread modelling/fire_spread/data/focal fires data/fire_data_raw_1999_2140469994_r-original.tif")
# velev <- values(fire_ori[["elev"]])
# na_elev_ids <- which(is.na(velev))
# length(na_elev_ids) / ncell(fire_ori)
# cols_na <- colFromCell(fire_ori, na_elev_ids) %>% unique
# cols_na <- cols_na[order(cols_na)]
# diff(cols_na) %>% unique # ok, they are consecutive
# max(cols_na) == ncol(fire_ori) # ok, its in the latest end
# fire_trim <- fire_ori[, -cols_na, drop = F]
# fire_trim
# writeRaster(fire_trim, file.path(gee_dir, "fire_data_raw_1999_2140469994_r-trimmed.tif"))
# # they were renamed manually, so the original says -original, and the
# # trimmed says nothing.
 
# # save elevation separately.
# fire_all <- rast(file.path(gee_dir, "fire_data_raw_2015_50.tif")) 
# fire_elev <- fire_all[["elev"]]
# writeRaster(fire_elev, file.path(windninja_dir, "2015_50.tif"), overwrite = TRUE)
 
# fire_all <- rast(file.path(gee_dir, "fire_data_raw_1999_2140469994_r.tif")) 
# fire_elev <- fire_all[["elev"]]
# writeRaster(fire_elev, file.path(windninja_dir, "1999_2140469994_r.tif"), overwrite = TRUE)

# # ñorquinco, remove northern rows
# fire_ori <- rast(file.path(gee_dir, "fire_data_raw_2014_1.tif"))
# velev <- values(fire_ori[["elev"]])
# na_elev_ids <- which(is.na(velev))
# length(na_elev_ids) / ncell(fire_ori)
# rows_na <- rowFromCell(fire_ori, na_elev_ids) %>% unique
# (rows_na <- rows_na[order(rows_na)])
# diff(rows_na) %>% unique # ok, they are consecutive
# fire_trim <- fire_ori[-rows_na, , drop = F]
# fire_trim
# writeRaster(fire_trim[["elev"]], file.path(windninja_dir, "2014_1.tif"),
#             overwrite = TRUE)

# Export elevations -------------------------------------------------------

fnames <- list.files(gee_dir)

for(i in 1:length(fnames)) {
  print(i)
  r <- rast(file.path(gee_dir, fnames[i]))[["elev"]]
  id_raw <- strsplit(fnames[i], "_")[[1]][-(1:3)]
  id <- paste(id_raw, collapse = "_")
  writeRaster(r, file.path(windninja_dir, id))
}

# Inputs for wind ninja (avg over study area):

# averages used before:

# // Wind direction avg:
# // 284.0693132595486°

# // Wind speed avg:
# // 4.78628962675556 m/s  ## this is much higher than the values I got now. it's weird.

# Note that the .tif files will be the elevation layer in elevation_dir,
# while the others will be outputs from windninja.

# Details for wind ninja:
# mesh resolution = 90 m,
# result scale = 30 m,
# speed unit = m/s,
# mass conservation solver
# speed input = varying by fire
# windir input = varying by fire

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
wind_table <- wind_table[order(wind_table$fire_id), c("fire_id", "direction_use",
                                                      "speed_mean")]
wind_table$speed_mean_ms <- wind_table$speed_mean / 3.6
wind_table <- wind_table[, c("fire_id", "direction_use", "speed_mean_ms",
                             "speed_mean")]
# write.csv(wind_table, file.path(
#   "data", "focal fires data", "wind ninja files", "zzz_wind_table.csv"
# ))


# Import wind speed and direction -----------------------------------------

wind_files_raw <- list.files(windninja_dir)
wind_dir_files <- wind_files_raw[grep("_ang.asc", wind_files_raw)]
wind_vel_files <- wind_files_raw[grep("_vel.asc", wind_files_raw)]
# view(cbind(wind_dir_files, wind_vel_files, fire_ids)) # ok

wind_imgs <- vector("list", n_fires)
for(i in 1:n_fires) {
  # i = 1
  dir <- rast(file.path(windninja_dir, wind_dir_files[i]))
  vel <- rast(file.path(windninja_dir, wind_vel_files[i]))
  wind <- c(dir, vel)
  names(wind) <- c("direction", "speed")
  wind_imgs[[i]] <- wind
}


# Import elevation, vfi, tfi ----------------------------------------------

vegtopo_imgs <- vector("list", n_fires)
for(i in 1:n_fires) {
  vegtopo_imgs[[i]] <- rast(file.path(gee_dir, fnames[i]))
}


# Import ignition points --------------------------------------------------

# ignition points for 2014_1 and 2008_5 (Ñorquinco and Lolog) were edited
# because these fires had a change of wind. The point to be used is not the 
# real ignition point, but is coherent with the dominant wind direction. 
# The original point id end in "_original", while the fire id was assigned to
# the edited points. These points are not included in the kml file.
points_raw <- vect(file.path("data", "ignition_points_checked.shp"))
points <- project(points_raw, vegtopo_imgs[[1]])

# Import FWI data ---------------------------------------------------------

fwi_raw <- read.csv(file.path("data", "climatic_data_by_fire_fwi-day-cumulative.csv"))
fwi_raw <- rename(fwi_raw, "fire_id_match" = "fire_id")

fwi_left <- data.frame(fire_id = fire_ids,
                       fire_id_match = fire_ids)

# rename repeted names (fires that were divided) to match
fwi_left$fire_id_match[grep("2015_47", fwi_left$fire_id)] <- "2015_47"
fwi_left$fire_id_match[grep("2011_19", fwi_left$fire_id)] <- "2011_19"

fwi_data <- left_join(fwi_left, fwi_raw, by = "fire_id_match")

# Make landscapes ---------------------------------------------------------

# // Landscape layers names in FireSpread package
# enum land_names {
#   vfi,
#   tfi,
#   elev,
#   wdir,
#   wspeed
# };

# (using the same order here)
land_names <- c("vfi", "tfi", "elev", "wdir", "wspeed")

# list with landscapes
lands <- vector(mode = "list", length = n_fires)
names(lands) <- fire_ids

# every fire will be a list
elem_names <- c("landscape", "burnable", "vegetation",  "ig_rowcol",
                "burned_layer", "burned_ids",
                "counts_veg", "counts_veg_available",
                "fwi", "fire_id")
# landscape, burnable and ig_rowcol are used to simulate the fire, while the remaining
# elements are used to compare with simulated ones (these are the same outputs
# from the fire simulation function).

for(i in 1:n_fires) {
  print(i)
  # i = 44

  lands[[i]] <- vector(mode = "list", length = length(elem_names))
  names(lands[[i]]) <- elem_names
  lands[[i]]$fire_id <- fire_ids[i]

  # vegetation data
  #                                 class code
  # 1                        Non burnable    1
  # 2                    Subalpine forest    2
  # 3                          Wet forest    3
  # 4                          Dry forest    4
  # 5                           Shrubland    5
  # 6                           Grassland    6 (turned into shrubland)
  # 7 Anthropogenic prairie and shrubland    7 (turned into shrubland)
  # 8                          Plantation    8 (turned into dry forest)
  
  # old transformation:
  
  # veg_vec[veg_codes %in% c(5, 7)] <- 0 # shrubland, to be explicit
  # veg_vec[veg_codes == 2] <- 1            # subalpine
  # veg_vec[veg_codes == 3] <- 2            # wet
  # veg_vec[veg_codes %in% c(4, 8)] <- 3    # dry
  # 
  # veg_vec[veg_codes == 1] <- 99 # non-burnable

  # project wind direction to match extent
  wind_local <- project(wind_imgs[[i]],
                        vegtopo_imgs[[i]],
                        method = "cubicspline")
  names(wind_local) <- c("wdir", "wspeed")

  # turn wind direction to radians
  wind_local$wdir <- wind_local$wdir * pi / 180
  
  # merge all data in a single raster
  rall <- c(vegtopo_imgs[[i]], wind_local)
  
  # make non-burnable the nan vegetation
  # rall$burnable[is.nan(rall$veg)] <- 0 
  
  # identify NA in the predictors to make those pixels unburnable
  vv <- values(rall)

  # make burnable the burned
  burned_cells <- which(vv[, "burned"] == 1)
  vv[burned_cells, "burnable"] <- 1
  
  # identify burnable cells with NA values
  na_cells <- which(apply(vv[, land_names], 1, anyNA) & 
                    vv[, "burnable"] == 1)
  
  # make the na non-burnable
  vv[na_cells, "burnable"] <- 0
  
  # edit values in raster
  values(rall) <- vv
  
  # caution if there are a lot of NA in burnable area
  na_prop <- length(na_cells) / ncell(rall)
  if(na_prop > 0.02) warning(paste("Many NA in landscape, fire: ", fire_ids[i], ", i: ", i, sep = ""))
  
  # make array for c++
  land_arr <- land_cube(rall)
  
  # subset in landscape and other stuff
  lands[[i]]$landscape <- land_arr[, , land_names]
  lands[[i]]$burnable <- land_arr[, , "burnable"]
  lands[[i]]$vegetation <- land_arr[, , "veg"] # raw codes!!!
  lands[[i]]$burned_layer <- land_arr[, , "burned"]
  
  # compute burned_ids (with 0-indexing!)
  burned_cells <- which(values(rall$burned) == 1)
  # length(burned_cells) == sum(v[, "burned"])
  burned_rowcols <- rowColFromCell(vegtopo_imgs[[i]], burned_cells)
  colnames(burned_rowcols) <- c("row", "col")
  lands[[i]]$burned_ids <- t(burned_rowcols) - 1 # for 0-indexing!

  # DO NOT COUNT VEG_VEC FOR NOW
  
  # # pixels by vegetation type in available and burned
  # veg_vec <- as.vector(land_arr[, , "veg"]) - 1 # 0 is non burnable now. 
  # n_veg_types <- length(unique(veg_vec[veg_vec > 0]))
  # counts_veg <- integer(n_veg_types)
  # counts_veg_available <- integer(n_veg_types)
  # for(k in 1:n_veg_types) {
  #   available <- as.numeric(veg_vec == (k-1))
  #   counts_veg[k] <- sum(available * v[, "burned"])
  #   counts_veg_available[k] <- sum(available)
  # }
  # 
  # lands[[i]]$counts_veg <- counts_veg
  # lands[[i]]$counts_veg_available <- counts_veg_available

  # get ignition points
  p_local <- points[points$Name == fire_ids[i]]
  
  # get coordinates and row_col
  cc <- crds(p_local)
  ig_rowcol <- rbind(rowFromY(vegtopo_imgs[[i]], cc[, "y"]),
                     colFromX(vegtopo_imgs[[i]], cc[, "x"]))
  row.names(ig_rowcol) <- c("row", "col")

  if(anyNA(ig_rowcol)) {
    stop(paste("Ignition point out of range,", "fire_id", fire_ids[i]))
  }

  lands[[i]]$ig_rowcol <- ig_rowcol - 1 # 0-indexing!!!

  # get FWI values 
  lands[[i]]$fwi <- fwi_data[fwi_data$fire_id == fire_ids[i], 
                             grep("fwi_", names(fwi_data))]
}

# Warning messages:
# 1: Many NA in landscape, fire: 2015_16, i: 37 
# 2: Many NA in landscape, fire: 2015_47N, i: 42 
# 3: Many NA in landscape, fire: 2015_47S, i: 43 
# 4: Many NA in landscape, fire: 2015_50, i: 44 



# Check all ignition points fall in burned and burnable cells,
# and all burned cells are burnable.
ccc <- numeric(n_fires)
# bbb <- numeric(n_fires)
for(i in 1:n_fires) {
  # i = 43
  ig <- lands[[i]]$ig_rowcol + 1 # because of 0-indexing

  point_checks <- sapply(1:ncol(ig), function(c) {
    # c = 1
    cbind(lands[[i]]$vegetation[ig[1, c], ig[2, c]] > 1, # 1 is non-burnable
          lands[[i]]$burned_layer[ig[1, c], ig[2, c]])
  }) %>% colSums %>% unique()

  ccc[i] <- point_checks

  if(point_checks != 2) {
    stop(paste("Ignition point problems,", "fire_id:", fire_ids[i], "i:", i))
  }
}
all(ccc == 2) # OK?
# all(bbb == 1)  # problem at fire i = 43

# problem at one fire, i = 15, id = 2005_6:
# row   50
# col   48
plot(vegtopo_imgs[[15]][[c("veg", "burned", "burnable")]])
# from earth, it's evident it was burnable, all burned shrubland.
bid_temp <- lands[["2005_6"]]$burned_ids + 1
for(c in 1:ncol(bid_temp)) {
  lands[["2005_6"]]$vegetation[bid_temp[1, c], bid_temp[2, c]] <- 5 # make shrubland
}




object.size(lands) / 1e6 # 1472.3 Mb (before it was 2634.7 Mb, when fwi was not
# constant and veg was a matrix).

object.size(lands) / 1e6 # 2988.9 Mb more info is provided now.

# sapply(lands, class) %>% unique # OK

# Save landscapes ---------------------------------------------------------

# saveRDS(lands, file.path("data", "landscapes_ig-known_non-steppe.rds"))
# lands <- readRDS(file.path("data", "landscapes_ig-known_non-steppe.rds"))

# save each landscape separately
for(i in 1:length(lands)) {
  print(i)
  saveRDS(lands[[i]], file.path("data", "focal fires data",
                                "landscapes_ig-known",
                                paste(lands[[i]]$fire_id, ".rds", sep = "")))
}
