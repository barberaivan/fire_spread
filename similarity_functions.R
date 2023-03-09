library(terra)

# Similarity/discrepancy functions ---------------------------------------

compare_fires_r <- function(fire1, fire2, lscale = 0.2) {

  # Extract list elements ------------------------------------------------

  burned1 <- fire1[["burned_layer"]]
  burned2 <- fire2[["burned_layer"]]

  burned_ids1 <- fire1[["burned_ids"]] + 1 # important to add 1 to use R indexing
  burned_ids2 <- fire2[["burned_ids"]] + 1

  size1 <- ncol(burned_ids1)
  size2 <- ncol(burned_ids2)

  counts1 <- fire1[["counts_veg"]]
  counts2 <- fire2[["counts_veg"]]

  # overlap_sp -----------------------------------------------------------

  # compute common pixels only in the smaller fire
  if(size1 < size2) {
    bb <- numeric(size1)
    for(i in 1:size1) {
      bb[i] <- burned2[burned_ids1[1, i], burned_ids1[2, i]]
    }
  } else {
    bb <- numeric(size2)
    for(i in 1:size2) {
      bb[i] <- burned1[burned_ids2[1, i], burned_ids2[2, i]]
    }
  }

  common <- sum(bb)
  overlap_sp <- common / (size1 + size2 - common)

  # overlap_vd -----------------------------------------------------------

  # Get vegetation distribution by fire (normalized burned areas)
  veg_types <- length(counts1)

  burned_dist_1 <- numeric(veg_types)
  burned_dist_2 <- numeric(veg_types)

  for(v in 1:veg_types) {
    burned_dist_1[v] <- counts1[v] / sum(counts1)
    burned_dist_2[v] <- counts2[v] / sum(counts2)
  }

  # compute vegetation distribution overlap
  overlap_vd <- sum(pmin(burned_dist_1, burned_dist_2))

  # deltas by veg_type ---------------------------------------------------

  # normalized difference using absolute difference. The difference by veg_type
  # is in [0, 1]. So, if we divide delta_norm by veg_num, it will be in [0, 1].
  sum_area <- counts1 + counts2
  use <- sum_area > 0

  delta_norm <- sum(abs((counts1[use] - counts2[use]) / sum_area[use]))

  # Scale to [0, 1]
  delta_norm_unit <- delta_norm / veg_types

  # Transform to similarities
  overlap_norm <- 1.0 - delta_norm_unit
  overlap_expquad <- exp(-delta_norm_unit ^ 2 / lscale) # 0.2 is the Gaussian SD.
  overlap_quad <- 1 - delta_norm_unit ^ 2

  # ---------------------------------------------------------------------

  indexes <- c(
    # pure indices
    "overlap_sp"      = overlap_sp,
    "overlap_vd"      = overlap_vd,
    "overlap_norm"    = overlap_norm,
    "overlap_expquad" = overlap_expquad,
    "overlap_quad"    = overlap_quad,

    # mixture indices
    "sp_norm_5050"    = 0.50 * overlap_sp + 0.50 * overlap_norm,
    "sp_norm_7525"    = 0.75 * overlap_sp + 0.25 * overlap_norm,
    "sp_expquad_5050" = 0.50 * overlap_sp + 0.50 * overlap_expquad,
    "sp_expquad_7525" = 0.75 * overlap_sp + 0.25 * overlap_expquad,
    "sp_quad_5050"    = 0.50 * overlap_sp + 0.50 * overlap_quad,
    "sp_quad_7525"    = 0.75 * overlap_sp + 0.25 * overlap_quad
  )

  return(indexes)
}