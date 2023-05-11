# This script prepares the landscapes where fires are going to be simulated.

# The first step is to extract the elevation layer from all images to use as
# input in windninja. These files will be saved as wind_*fire id*.tif.

# Once the wind layers are created, the landscape arrays will be made as follows:

# the vegetation_type will be an integer matrix, with the following values:

#   CLASS                        VALUE

#´  shrubland                      0
#´  subalpine forest               1
#´  wet forest                     2
#´  dry forest                     3

#´  non-burnable                   99

## IMPORTANT: the vegetation layer used here does not distinguish araucaria vs
## cypress dry forests. To separate them, the download from GEE should be 
## updated.


# The following variables will be stored in a 3D array, with the following
# layers:

#´   aspect            [-1, 1] slope-weighted northing
#´   wind direction    [-1, 1] windward vs leeward
#´   elevation,        (standardized; mean = 1163.3 m, sd = 399.5)

# (slope will be computed in the simulation, as it's a directional variable.)

# The FWI, as an anomaly standardized at the pixel level, will be held constant
# within each fire, taking its value at the ignition point. This variable will
# model a fire-level random intercept, so it's not considered by the fire-spread
# functions. 


library(terra)
library(tidyverse)

# get paths
gee_dir <- file.path("data", "focal fires data", "raw data from GEE")
windninja_dir <- file.path("data", "focal fires data", "wind ninja files")

# elevation constants (to standardize variable)
elev_mean <- 1163.3
elev_sd <- 399.5

# Number of vegetation types
n_veg_types <- 4 # (shrubland, subalpine, wet, dry)

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

# // Wind direction avg:
# // 284.0693132595486°
#
# // Wind speed avg:
# // 4.78628962675556 m/s

# Note that the .tif files will be the elevation layer in elevation_dir,
# while the others will be outputs from windninja.

# Details for wind ninja:
# mesh resolution = 90 m,
# result scale = 30 m,
# speed unit = m/s,
# speed input = 4.79 m/s,
# windir input = 284°
# mass conservation solver


# Import raw data, wind and fwi -------------------------------------------

# import fire shapes to get year (for fwi matching)
f <- vect(file.path("data", "patagonian_fires_spread.shp"))
# _spread has the fires edited so every burned separate patch counts as a
# separate fire, no matter whether they correspond to the same event or not.

# fwi image
fwi <- rast(file.path("data", "fwi_anomalies.tif"))
crs(fwi) <- "EPSG:4326"

fnames <- list.files(gee_dir)
fire_ids <- fnames
n_fires <- length(fire_ids)

# list with raw images
raw_imgs <- vector(mode = "list", length = length(fnames))
for(i in 1:length(fnames)) {
  raw_imgs[[i]] <- rast(file.path(gee_dir, fnames[i]))
  id_raw <- strsplit(fnames[i], c("_|[.]"))[[1]]
  remove <- c(1:3, length(id_raw))
  id <- paste(id_raw[-remove], collapse = c("_"))
  fire_ids[i] <- id
}
names(raw_imgs) <- fire_ids

# list with wind direction rasters

wind_files_raw <- list.files(windninja_dir)
wind_files <- wind_files_raw[grep("_ang.asc", wind_files_raw)]
wind_ids <- sapply(wind_files, function(x) {
  id_raw <- strsplit(x, c("_|[.]"))[[1]]
  remove <- length(id_raw) : (length(id_raw) - 4)
  id <- paste(id_raw[-remove], collapse = c("_"))
  return(id)
}) %>% unname()
# all.equal(wind_ids, fire_ids)

wind_imgs <- raw_imgs
for(i in 1:length(wind_imgs)) {
  fii <- file.path(windninja_dir, wind_files[i])
  wind_imgs[[i]] <- rast(fii)
}


# Import ignition points ---------------------------------------------------------

points_raw <- vect(file.path("data", "ignition_points_checked.shp"))
points <- project(points_raw, raw_imgs[[1]])

# Make landscapes ---------------------------------------------------------

# get years
fyears <- sapply(fire_ids, function(x) {
  f$year[f$fire_id == x]
})

# list with landscapes
lands <- vector(mode = "list", length = n_fires)
names(lands) <- fire_ids

for(i in 1:n_fires) {
  print(i)
  # i = 1

  # every fire will be a list
  elem_names <- c("terrain", "vegetation",  "ig_rowcol", 
                  "burned_layer", "burned_ids", 
                  "counts_veg", "counts_veg_available",
                  "fwi", "year", "fire_id")
  # terrain and ig_rowcol are used to simulate the fire, while the remaining
  # elements are used to compare with simulated ones (these are the same outputs
  # from the fire simulation function).

  lands[[i]] <- vector(mode = "list", length = length(elem_names))
  names(lands[[i]]) <- elem_names
  lands[[i]]$year <- fyears[i] %>% unname
  lands[[i]]$fire_id <- fire_ids[i]
  
  # get raw image values
  v <- values(raw_imgs[[i]])
  
  # vegetation
  veg_codes <- v[, "veg"]
  # replace nan with 1 (non burnable)
  veg_codes[is.nan(veg_codes)] <- 1

  # make veg matrix
  veg_vec <- integer(length(veg_codes)) 
  
  # vegetation data
  #                                 class code
  # 1                        Non burnable    1
  # 2                    Subalpine forest    2
  # 3                          Wet forest    3
  # 4                          Dry forest    4
  # 5                           Shrubland    5
  # 6                           Grassland    6 (turned into shrubland)
  # 7 Anthropogenic prairie and shrubland    7 (turned into shrubland)
  # 8                          Plantation    8 (turned into dry forest B)
  
  veg_vec[veg_codes %in% c(5, 6, 7)] <- 0 # shrubland, to be explicit
  veg_vec[veg_codes == 2] <- 1            # subalpine
  veg_vec[veg_codes == 3] <- 2            # wet
  veg_vec[veg_codes %in% c(4, 8)] <- 3    # dry
  
  veg_vec[veg_codes == 1] <- 99 # non-burnable
  
  # project wind direction to match extent
  wind_local <- project(wind_imgs[[i]],
                        raw_imgs[[i]],
                        method = "cubicspline")
  names(wind_local) <- "windir"
  
  # compute slope-weighted northing
  # (see code <northing importance function.R>)
  northing <- cos(v[, "aspect"] * pi / 180) * plogis(-5 + 0.35 * v[, "slope"])
  
  # get landscape in matrix form
  land_long <- cbind(
    veg = veg_vec,
    northing = northing, 
    windir = values(wind_local) * pi / 180, # in radians
    elev = (v[, "elev"] - elev_mean) / elev_sd, # standardized
    burned = v[, "burned"]
  )

  # get landscape in array form
  land_arr <- array(
    NA,
    dim = c(nrow(raw_imgs[[i]]),
            ncol(raw_imgs[[i]]),
            ncol(land_long)),
    dimnames = list(rows = NULL, cols = NULL,
                    layers = colnames(land_long))
  )

  for(j in 1:ncol(land_long)) {
    land_arr[, , j] <- matrix(land_long[, j],
                              nrow = nrow(land_arr),
                              ncol = ncol(land_arr),
                              byrow = TRUE)
  } # terra gives raster values by row

  # Fill terrain, vegetation matrix and burned_layer
  lands[[i]]$terrain <- land_arr[, , c("northing", "windir", "elev")] 
    # Order matters. Wind and elevation must be the last layers.
  lands[[i]]$vegetation <- land_arr[, , "veg"]
  lands[[i]]$burned_layer <- land_arr[, , "burned"]

  # compute burned_ids (with 0-indexing!)
  burned_cells <- which(v[, "burned"] == 1)
  # length(burned_cells) == sum(v[, "burned"])
  burned_rowcols <- rowColFromCell(raw_imgs[[i]], burned_cells)
  colnames(burned_rowcols) <- c("row", "col")
  lands[[i]]$burned_ids <- t(burned_rowcols) - 1 # for 0-indexing!

  # pixels by vegetation type in available and burned
  counts_veg <- integer(n_veg_types)
  counts_veg_available <- integer(n_veg_types)
  for(k in 1:n_veg_types) {
    available <- as.numeric(veg_vec == (k-1))
    counts_veg[k] <- sum(available * v[, "burned"]) 
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
  
  # get FWI values (mean, median and sd are used for additional analyses apart
  # from model fittng)
  fwi_local <- project(fwi[[as.character(fyears[i])]],
                       raw_imgs[[i]],
                       method = "cubicspline")
  names(fwi_local) <- "fwi"
  
  fwi_vals <- values(fwi_local)
  ig_cell <- cellFromRowCol(fwi_local, ig_rowcol["row", ], ig_rowcol["col", ])
  fwi_data <- c("ig" = mean(fwi_vals[ig_cell]),
                "mean" = mean(fwi_vals),
                "median" = median(fwi_vals),
                "sd" = sd(fwi_vals),
                "min" = min(fwi_vals),
                "max" = max(fwi_vals))
  lands[[i]]$fwi <- fwi_data
}

# Check all ignition points fall in burned and burnable cells,
# and all burned cells are burnable.
ccc <- numeric(n_fires)
bbb <- numeric(n_fires)
for(i in 1:n_fires) {
  # i = 14
  ig <- lands[[i]]$ig_rowcol + 1 # because of 0-indexing
  
  point_checks <- sapply(1:ncol(ig), function(c) {
    # c = 1
    cbind(lands[[i]]$vegetation[ig[1, c], ig[2, c]] < 99, # 99 is non-burnable
          lands[[i]]$burned_layer[ig[1, c], ig[2, c]])
  }) %>% colSums %>% unique()
  
  ccc[i] <- point_checks
  
  if(point_checks != 2) {
    stop(paste("Ignition point problems,", "fire_id:", fire_ids[i], "i:", i))
  }
  
  # check burned is burnable
  matches <- sum(lands[[i]]$burned_layer * (lands[[i]]$vegetation < 99))
  bbb[i] <- sum(lands[[i]]$counts_veg) == matches
}
all(ccc == 2) # OK?
all(bbb == 1)

# problem at one fire, i = 14, id = 2005_6:
# row   50
# col   48
plot(raw_imgs[["2005_6"]][[c("veg", "burned")]])
# from earth, it's evident it was burnable, all burned shrubland.
bid_temp <- lands[["2005_6"]]$burned_ids + 1
for(c in 1:ncol(bid_temp)) {
  lands[["2005_6"]]$vegetation[bid_temp[1, c], bid_temp[2, c]] <- 0
}
# correct the counts_veg from this fix.
for(k in 1:n_veg_types) {
  available <- as.numeric(lands[["2005_6"]]$vegetation == (k-1))
  lands[["2005_6"]]$counts_veg[k] <- sum(available * as.numeric(lands[["2005_6"]]$burned_layer)) 
  lands[["2005_6"]]$counts_veg_available[k] <- sum(available)
}

object.size(lands) / 1e6 # 1472.3 Mb (before it was 2634.7 Mb, when fwi was not 
                         # constant and veg was a matrix.
# sapply(lands, class) %>% unique # OK

# Save landscapes list ----------------------------------------------------

saveRDS(lands, file.path("data", "landscapes_ig-known_non-steppe.rds"))
# lands <- readRDS(file.path("data", "landscapes_ig-known_non-steppe.rds"))

# save each landscape separately
for(i in 1:length(lands)) {
  saveRDS(lands[[i]], file.path("data", "focal fires data",
                                "landscapes_ig-known_non-steppe",
                                paste(lands[[i]]$fire_id, ".rds", sep = "")))
}


# Export csv files to import them in C++ ----------------------------------

## WARNING, the following code follows an outdated data structure.

# File names are composed as <[FIRE_ID]-[OBJECT_NAME].csv>, and were produced
# in the <landscapes_preparation.R> script. All files have column names. 
# There is one csv by fire and by object type, which are the following:

# landscape.csv: landscape data, including the burnable layer (not the burned).
#   Each landscape layer (matrix) is vectorized by column. The csv has a row by 
#   pixel and a column by variable (landscape layer). 
#   NOTE: This is not how the terra R package converts the matrices into vectors.
#   They fill the vector by row, not by column.
# burned_ids.csv: ids of the burned cells in the observed fire, as {row, column} 
#   vectors. Every row is a burned pixel, first column is the row id, the second
#   is the column id. Row and column indexes of each pixel start at 0.
# ignition_points.csv: similar to burned_ids, but only containing the pixels were
#   the observed fire begun (usually just one).
# metadata.csv: a column by variable, as follows
#   size_rows = landscape size in rows,
#   size_cols = landscape size in columns,
#   counts_veg_shrubland = number of burned pixels in shrubland vegetation,
#   counts_veg_subalpine = number of burned pixels in subalpine vegetation,
#   counts_veg_wet = number of burned pixels in wet forest vegetation,
#   counts_veg_dry = number of burned pixels in dry forest vegetation.

# In addition, all fire lists will be saved separately, only named by the fire id.
# (useful for model fitting)

fires <- readRDS(file.path("..", "fire_spread_data", "landscapes_ig-known_non-steppe.rds"))
target_dir <- file.path("..", "fire_spread_data", "focal fires data", 
                        "landscapes_ig-known_non-steppe_csv-files")

for(f in 1:length(fires)) {
  print(f)
  # save separate fire lists
  fname <- paste(names(fires)[f], ".rds", sep = "")
  # saveRDS(
  #   fires[[f]],
  #   file.path("..", "fire_spread_data", "focal fires data", 
  #             "landscapes_ig-known_non-steppe", fname)
  # )
  
  # write landscape in vector form (fill vector by column)
  l <- do.call("cbind", lapply(1:8, function(x) {
    fires[[f]]$landscape[, , x] %>% as.numeric
    }
  ))
  colnames(l) <- dimnames(fires[[f]]$landscape)[[3]][1:8]
  write.csv(l, file.path(target_dir, 
                         paste(names(fires)[f], "landscape.csv", sep = "-")),
            row.names = FALSE)
  
  # write burned_ids
  bid <- t(fires[[f]]$burned_ids)
  colnames(bid) <- c("row_id", "col_id")  
  write.csv(bid, 
            file.path(target_dir, 
                      paste(names(fires)[f], "burned_ids.csv", sep = "-")),
            row.names = FALSE)
  
  # write ignition point
  igp <- t(fires[[f]]$ig_rowcol)
  colnames(igp) <- c("row_id", "col_id")  
  write.csv(igp, 
            file.path(target_dir, 
                      paste(names(fires)[f], "ignition_points.csv", sep = "-")),
            row.names = FALSE)
  
  # metadata
  metadata <- data.frame(
    size_rows = dim(fires[[f]]$landscape)[1],
    size_cols = dim(fires[[f]]$landscape)[2], 
    counts_veg_shrubland = fires[[f]]$counts_veg[1],
    counts_veg_subalpine = fires[[f]]$counts_veg[2],
    counts_veg_wet = fires[[f]]$counts_veg[3],
    counts_veg_dry = fires[[f]]$counts_veg[4]
  )
  
  write.csv(metadata, 
            file.path(target_dir, 
                      paste(names(fires)[f], "metadata.csv", sep = "-")),
            row.names = FALSE)
}