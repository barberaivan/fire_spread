# Packages ----------------------------------------------------------------

library(terra)
source(file.path("weather", "fortnight_functions.R"))

# Load PNNH to compare ----------------------------------------------------

# Nahuel Huapi National Park
pnnh <- vect(file.path("data", "protected_areas", "apn_limites.shp"))
pnnh <- pnnh[pnnh$nombre == "Nahuel Huapi", ]
# pnnh <- project(pnnh, "EPSG:5343")

# Detrend FWI: get only temporal anomalies from original raster -----------

# daily FWI stack
# fwi_day <- rast(file.path("data", "fwi_daily_1998-2022", "52km",
#                           "fwi_daily_19980101_20230630.tif"))
fwi_day <- rast(file.path("data", "fwi_daily_1998-2022", "24km",
                          "fwi_daily_19970701_20230630.tif"))

# compute temporal mean and sd using the [1998-07-01, 2023-06-30] time range
# (25 years) (it used since 1999-07-01 before)
date_begin <- as.Date("1998-07-01", format = "%Y-%m-%d")
date_end <- as.Date("2023-06-30", format = "%Y-%m-%d")

# dates_keep <- which(time(fwi_day) >= date_begin & time(fwi_day) <= date_end + 1)
dates_avg <- which(time(fwi_day) >= date_begin & time(fwi_day) <= date_end + 1)
# plus one because it removed the last day

# Filter dates to have only complete years (99-22)
# fwi_day_use <- fwi_day[[dates_keep]]
fwi_day_use <- fwi_day[[dates_avg]]
# round(nlyr(fwi_day_use) / 365) # OK, 24 years

# Compute mean and sd across layers
fwi_mean_rast <- mean(fwi_day_use)
fwi_sd_rast <- app(fwi_day_use, sd)

# standardize
fwiz_day <- (fwi_day - fwi_mean_rast) / fwi_sd_rast

# writeRaster(
#   fwiz_day,
#   file.path("data", "fwi_daily_1998-2022", "52km",
#             "fwi_daily_19980101_20230630_standardized.tif")
# )
writeRaster(
  fwiz_day,
  file.path("data", "fwi_daily_1998-2022", "24km",
            "fwi_daily_19970701_20230630_standardized.tif")
)

# # check:
# plot(fwiz_day, c(15, 170), range = c(-1.5, 2.5))


# Aggregate FWI by fortnight ----------------------------------------------

# daily FWI stack (standardized)
# fwi_day <- rast(file.path("data", "fwi_daily_1998-2022", "52km",
#                           "fwi_daily_19980101_20230630_standardized.tif"))
fwi_day <- rast(file.path("data", "fwi_daily_1998-2022", "24km",
                          "fwi_daily_19970701_20230630_standardized.tif"))
fwi_day_dates <- time(fwi_day)
fwi_day_dates <- as.Date(fwi_day_dates, format = "%Y-%m-%d")
fwi_day_forts <- date2fort(fwi_day_dates)
fwi_forts <- unique(fwi_day_forts)
# aggregate raster by fortnight
fwi_fort <- do.call("c", lapply(fwi_forts, function(f) {
  print(f)
  ids <- which(fwi_day_forts == f)
  return(mean(fwi_day[[ids]]))
})) # this takes long
# Assign the first day to every fortnight, so the fortnight number does not
# depend on the reference table.
time(fwi_fort) <- fort2date(fwi_forts)
names(fwi_fort) <- fwi_forts

# writeRaster(
#   fwi_fort,
#   file.path("data", "fwi_daily_1998-2022", "52km",
#             "fwi_fortnights_19980101_20230630_standardized.tif")
# )
writeRaster(
  fwi_fort,
  file.path("data", "fwi_daily_1998-2022", "24km",
            "fwi_fortnights_19970701_20230630_standardized.tif")
)
