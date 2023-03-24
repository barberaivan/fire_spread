# This script prepares the landscapes where fires are going to be simulated.

# The first step is to extract the elevation layer from all images to use as
# input in windninja. These files will be saved as wind_*fire id*.tif.

# Once the wind layers are created, the landscape arrays will be made as follows:

#´   subalpine forest, {0, 1} (shrubland goes in the intercept)
#´   wet forest,       {0, 1}
#´   dry forest,       {0, 1}
#´   fwi,              (-Inf, +Inf) (standardized anomaly at pixel level)
#´   aspect            [-1, 1] # Northwestyness
#´   wind direction    [-1, 1] (windward vs leeward)
#´   elevation,        (standardized in the code)

library(terra)
library(tidyverse)

# get paths
data_path <- file.path("..", "fire_spread_data")
gee_dir <- file.path(data_path, "focal fires data", "raw data from GEE")
windninja_dir <- file.path(data_path, "focal fires data", "wind ninja files")


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
f <- vect(file.path(data_path, "patagonian_fires_spread.shp"))
# _spread has the fires edited so every burned separate patch counts as a
# separate fire, no matter whether they correspond to the same event or not.

# fwi image
fwi <- rast(file.path(data_path, "fwi_anomalies.tif"))
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


# Make landscapes ---------------------------------------------------------


# vegetation data
#                                 class code
# 1                        Non burnable    1
# 2                    Subalpine forest    2
# 3                          Wet forest    3
# 4                          Dry forest    4
# 5                           Shrubland    5
# 6                           Grassland    6 (turn into shrubland)
# 7 Anthropogenic prairie and shrubland    7 (turned into shrubland)
# 8                          Plantation    8 (turned into dry forest)


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
  elem_names <- c("landscape", "ig_rowcol",
                  "burned_layer", "burned_ids", "counts_veg")
  # landscape and ig_rowcol are used to simulate the fire, while the remaining
  # elements are used to compare with simulated ones (these are the same outputs
  # from the fire simulation function).

  lands[[i]] <- vector(mode = "list", length = 5)
  names(lands[[i]]) <- elem_names

  # get raw image values
  v <- values(raw_imgs[[i]])

  # vegetation
  veg_codes <- v[, "veg"]
  # replace nan with 1 (non burnable)
  veg_codes[is.nan(veg_codes)] <- 1

  # make veg matrix
  veg_mat <- matrix(0, length(veg_codes), 3) # 3 cols for wet, subalpine, and dry forests
  # shrubland is the reference
  colnames(veg_mat) <- c("subalpine", "wet", "dry")

  veg_mat[veg_codes == 2, 1] <- 1 # subalpine
  veg_mat[veg_codes == 3, 2] <- 1 # wet
  veg_mat[veg_codes == 4, 3] <- 1 # dry

  # get fwi
  fwi_local <- project(fwi[[as.character(fyears[i])]],
                       raw_imgs[[i]],
                       method = "cubicspline")
  names(fwi_local) <- "fwi"

  # project wind direction to match extent
  wind_local <- project(wind_imgs[[i]],
                        raw_imgs[[i]],
                        method = "cubicspline")
  names(wind_local) <- "wind"

  # get landscape in matrix form
  land_long <- cbind(
    veg_mat,
    fwi = values(fwi_local),
    aspect = cos(v[, "aspect"] * pi / 180), # northing
    windir = values(wind_local) * pi / 180, # in radians
    elev = v[, "elev"],
    burnable = v[, "burnable"],
    burned = v[, "burned"]
  )
  # Order matters. Wind and elevation must be the last data-layers,
  # with burnable and burned in the end.

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

  # Fill landscape and burned_layer
  lands[[i]]$landscape <- land_arr[, , 1:(ncol(land_long) - 1)] # without burned layer
  lands[[i]]$burned_layer <- land_arr[, , ncol(land_long)]

  # compute burned_ids (with 0-indexing!)
  burned_cells <- which(v[, "burned"] == 1)
  # length(burned_cells) == sum(v[, "burned"])
  burned_rowcols <- rowColFromCell(raw_imgs[[i]], burned_cells)
  lands[[i]]$burned_ids <- t(burned_rowcols) - 1 # for 0-indexing!

  # burned pixels by vegetation type
  counts_veg <- numeric(4)
  for(k in 2:4) counts_veg[k] <- sum(veg_mat[, (k-1)] * v[, "burned"]) # non-shrubland
  counts_veg[1] <- sum(v[, "burned"]) - sum(counts_veg[2:4]) # the remaining is shrubland
  # sum(counts_veg) == sum(v[, "burned"])
  lands[[i]]$counts_veg <- counts_veg
}


object.size(lands) / 1e6 # 2634.7 Mb
# saveRDS(lands, "lands_temp.rds")
# sapply(lands, class) %>% unique # OK


# Ignition points ---------------------------------------------------------

points_raw <- vect(file.path(data_path, "ignition_points_checked.shp"))
points <- project(points_raw, raw_imgs[[1]])

for(i in 1:n_fires) {
  # i = 1

  # subset ignition points
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


# Check all ignition points fall in burned and burnable cells
ccc <- numeric(n_fires)
for(i in 1:n_fires) {
  # i = 45
  ig <- lands[[i]]$ig_rowcol + 1 # because of 0-indexing

  point_checks <- sapply(1:ncol(ig), function(c) {
    cbind(lands[[i]]$landscape[ig[1, c], ig[2, c], "burnable"],
          lands[[i]]$burned_layer[ig[1, c], ig[2, c]])
  }) %>% colSums %>% unique()

  ccc[i] <- point_checks

  if(point_checks != 2) {
    stop(paste("Ignition point problems,", "fire_id:", fire_ids[i], "i:", i))
  }
}
ccc # perfect (must be 2).


# Save landscapes list ----------------------------------------------------

saveRDS(lands, file.path(data_path, "landscapes_ig-known_non-steppe.rds"))



# OLD CODE BELOW, --------------------------------------------------------
# used to make cholila and another fires before.

# Cholila fire ------------------------------------------------------------

land0 <- rast(paste(data_dir, "data_cholila.tif", sep = ""))
elev <- rast(paste(data_dir, "data_cholila_elevation.tif", sep = ""))
angles_90m <- rast(paste(data_dir, "data_cholila_elevation_284_10_90m_ang.asc", sep = ""))

# interpolate windir to 30 m
windir <- project(angles_90m, elev, method = "cubicspline")



veg_codes <- values(land0$veg)
unique(veg_codes)
# replace nan with 1 (non burnable)
veg_codes[is.nan(veg_codes)] <- 1

veg_mat <- matrix(0, length(veg_codes), 3) # 3 cols for wet, subalpine, and dry forests
                                            # shrubland is the reference
colnames(veg_mat) <- c("subalpine", "wet", "dry")

veg_mat[veg_codes == 3, 1] <- 1 # subalpine
veg_mat[veg_codes == 3, 2] <- 1 # wet
veg_mat[veg_codes == 4, 3] <- 1 # dry

burnable <- as.numeric(veg_codes != 1)

# make landscape
land0_mat <- values(land0)
colnames(land0_mat)

# IMPORTANT: ADD BURNABLE AND BURNED LAYERS IN THE END

landscape <- cbind(
  veg_mat,
  fwi = rep(0, nrow(land0_mat)),
  aspect = cos(land0_mat[, "aspect"] * pi / 180), # northing
  windir = values(windir)[, 1],
  elev = land0_mat[, "elev"],
  burnable = burnable,
  burned = land0_mat[, "burned"]
)

head(landscape)
dim(landscape)
object.size(landscape) / 1e6

sum(landscape[, "burnable"]) / nrow(landscape) # 68 % quemable

# save cholila landscape

saveRDS(landscape, paste(data_dir, "data_cholila_landscape.R", sep = ""))





# 2021-865 landscape ------------------------------------------------------

data_dir <- "/home/ivan/Insync/Fire spread modelling/data/focal fires data/"

land0 <- rast(paste(data_dir, "data_2021_865.tif", sep = ""))
# elev <- rast(paste(data_dir, "data_2021_865_elevation.tif", sep = ""))
# angles_90m <- rast(paste(data_dir, "data_2021_865_elevation_284_10_90m_ang.asc", sep = ""))

# interpolate windir to 30 m
# windir <- project(angles_90m, elev, method = "cubicspline")

# make vegetation data

#                                 class code
# 1                        Non burnable    1
# 2                    Subalpine forest    2
# 3                          Wet forest    3
# 4                          Dry forest    4
# 5                           Shrubland    5
# 6                           Grassland    6 (turn into shrubland)
# 7 Anthropogenic prairie and shrubland    7 (turned into shrubland)
# 8                          Plantation    8 (turned into dry forest)

veg_codes <- values(land0$veg)
unique(veg_codes)
# replace nan with 1 (non burnable)
veg_codes[is.nan(veg_codes)] <- 1

veg_mat <- matrix(0, length(veg_codes), 3) # 3 cols for wet, subalpine, and dry forests
# shrubland is the reference
colnames(veg_mat) <- c("subalpine", "wet", "dry")

veg_mat[veg_codes == 3, 1] <- 1 # subalpine
veg_mat[veg_codes == 3, 2] <- 1 # wet
veg_mat[veg_codes == 4, 3] <- 1 # dry

burnable <- as.numeric(veg_codes != 1)

# make landscape
land0_mat <- values(land0)
colnames(land0_mat)

# IMPORTANT: ADD BURNABLE AND BURNED LAYERS IN THE END

landscape <- cbind(
  veg_mat,
  fwi = rep(0, nrow(land0_mat)),
  aspect = rep(0, nrow(land0_mat)), #cos(land0_mat[, "aspect"] * pi / 180), # northing
  windir = rep(0, nrow(land0_mat)), #values(windir)[, 1],
  elev = land0_mat[, "elev"],
  burnable = burnable,
  burned = land0_mat[, "burned"]
)

head(landscape)
dim(landscape)
object.size(landscape) / 1e6

sum(landscape[, "burnable"]) / nrow(landscape) # 84 % quemable

# save 2021_865 landscape
saveRDS(landscape, paste(data_dir, "data_2021_865_landscape.rds", sep = ""))



# Export csv files to import them in C++ ----------------------------------

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