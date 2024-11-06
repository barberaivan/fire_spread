# Functions to compute VFI and TFI (vegetation and topographic flammability
# indices). It also includes the function to detrend NDVI.

# Data --------------------------------------------------------------------

# parameters to compute flammability indices
fi_params <- readRDS(file.path("data", "flammability indices",
                               "flammability_indices.rds"))

# model to detrend ndvi
mdetrend <- readRDS(file.path("data", "flammability indices",
                              "ndvi_detrender_model.rds"))

# Functions ---------------------------------------------------------------

# detrend NDVI, convert to 2022-equivalent.
# Remember that FireSpread takes the NDVI from the previous summer, so this
# function should take year = fire_year - 1 if the output is desired for
# fire_year.
ndvi_detrend <- function(ndvi_focal, year, ndvi_22) {
  if(year == 2022) return(ndvi_focal)

  dpred_ndvi <- data.frame(
    ndvi_dyn_logit = qlogis((ndvi_focal + 1) / 2),
    ndvi01_22 = (ndvi_22 + 1) / 2,
    year = year
  )

  dpred_ndvi$diff_logit <- mgcv::predict.gam(mdetrend, dpred_ndvi, se.fit = F)
  ndvi_dt <- plogis(dpred_ndvi$ndvi_dyn_logit - dpred_ndvi$diff_logit) * 2 - 1
  return(ndvi_dt)
}

# Compute vegetation flammability index, from vegetation vector (integer in 1:5)
# and the detrended NDVI. Vegetation classes are
# {
#   1: Wet forest,
#   2: Subalpine forest,
#   3: Dry forest,
#   4: Shrubland,
#   5: Grassland
# }
vfi_calc <- function(vegetation, ndvi) {
  vfi_raw <- numeric(length(ndvi))

  for(v in unique(vegetation)) {
    rows <- vegetation == v
    vfi_raw[rows] <-
      fi_params$a[v] +
      fi_params$b[v] *
      (ndvi[rows] - fi_params$o[v]) ^ 2
  }

  # standardize
  return((vfi_raw - fi_params$vfi_mean) / fi_params$vfi_sd)
}


# Northing importance function, so this is weighted by slope (°).
northing_weight <- function(slope) {
  plogis(-5 + 0.35 * slope)
}

# Compute northing term, from aspect and slope, both in degrees
northing_calc <- function(aspect, slope) {
  cos(aspect * pi / 180) * northing_weight(slope)
}

# Compute topographic flammability index, from elevation (m), aspect (°), and
# slope (°)
tfi_calc <- function(elevation, aspect, slope) {
  northing <- northing_calc(aspect, slope)
  tfi_raw <- fi_params$b_elev_ori * elevation +
             fi_params$b_north_ori * northing
  return((tfi_raw - fi_params$tfi_mean) / fi_params$tfi_sd)
}