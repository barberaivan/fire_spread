# Compute FWI summer-level anomalies.

# 1) Compute summer mean for every fire season (December-March).
# 2) Get summers mean and sd considering the whole period.
# 3) Turn summer values into standardized anomalies: [x - mean(x)] / sd(x)

library(terra)
library(tidyverse)

# Define directories

data_path <- file.path("..", "fire_spread_data")
source_dir <- file.path(data_path, "dataset-cems-fire-historical-fwi-1998-1998-dic-march-all-days")
target_dir <- data_path

# select files to import

sfiles <- list.files(source_dir)
length(sfiles)
head(sfiles)

ymd <- sapply(sfiles, function(x) {
  strsplit(x, "_")[[1]][4]
})

fdata <- data.frame(name = sfiles, 
                    ymd = as.character(ymd))
fdata$y <- sapply(fdata$ymd, function(x) {
  substring(x, c(1, 5, 7), c(4, 6, 7))[1]
})

fdata$m <- sapply(fdata$ymd, function(x) {
  substring(x, c(1, 5, 7), c(4, 6, 7))[2]
})

fdata$m <- as.numeric(fdata$m)
fdata$y <- as.numeric(fdata$y)

fdata$fseason <- fdata$y
fdata$fseason[fdata$m == 12] <- fdata$fseason[fdata$m == 12] + 1

# remove extra data from 1998
out <- which(fdata$y == 1998 & fdata$m < 4)
fdata <- fdata[-out, ]


# Create yearly summer averages -------------------------------------------

years <- unique(fdata$fseason)
n_years <- length(years)
avg_list <- vector(mode = "list", length = n_years)

for(i in 1:n_years) {
  y <- years[i]
  print(y)
  ids <- which(fdata$fseason == y)
  
  daily_list <- vector(mode = "list", length = length(ids))
  
  # import daily images
  for(d in 1:length(ids)) {
    file_name <- file.path(source_dir, fdata$name[ids[d]])
    daily_list[[d]] <- rast(file_name)
  }
  
  # set as raster with a layer by day
  focal_year <- rast(daily_list)
  
  # compute focal_year mean
  avg_list[[i]] <- mean(focal_year)
}

# check:
for(i in 1:n_years) plot(avg_list[[i]], main = years[i])


# Turn into anomalies -----------------------------------------------------

global_mean <- mean(rast(avg_list))
global_sd <- stdev(rast(avg_list))
anom_list <- avg_list 

for(i in 1:n_years) {
  anom_list[[i]] <- (avg_list[[i]] - global_mean) / global_sd
}
anom_raster <- rast(anom_list)
names(anom_raster) <- years
anom_raster
writeRaster(anom_raster, 
            file.path(target_dir, "fwi_anomalies.tif"))
