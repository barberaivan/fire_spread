# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(deeptime)
library(rstan)
library(terra)
library(tidyterra)
library(FireSpread)
library(foreach)       # parallelization
library(doMC)          # parallelization
library(stringr)       # add leading zeroes to file

source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for rast_from_mat and a few constants

source(file.path("flammability indices",
                 "flammability_indices_functions.R"))

# Figure size settings ----------------------------------------------------

a4h <- 29.7
a4w <- 21.0
margins <- 2.5
fig_width_max <- a4w - margins * 2
fig_height_max <- a4h - margins * 2

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

normalize <- function(x) x / sum(x)

mean_ci <- function(x) {
  qq <- quantile(x, probs = c(0.025, 0.975), method = 8) %>% unname
  return(c("mean" = mean(x), "lower" = qq[1], "upper" = qq[2]))
}

nice_theme <- function() {
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),

    axis.line = element_line(linewidth = 0.3),

    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),

    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "white", color = "white")
  )
}

# The inverse of unconstrain (constrain spread parameters to bounded support)
constrain <- function(xun, support) {
  xc <- xun

  names_logit <- colnames(xun)[colnames(xun) != "steps"]
  for(j in names_logit) {
    xc[, j] <- plogis(xun[, j]) * (support[2, j] - support[1, j]) + support[1, j]
  }

  xc[, "steps"] <- exp(xun[, "steps"])

  return(xc)
}

# Make clipped landscape to simulate a fire using only the necessary RAM.
# It takes a {row, col} vector (1-indexing) and the steps parameter.
# Returns the clipped landscape with row_shift and col_shift as attributes.
# Adding these values to the burned_ids translates the row-col to the original
# landscape, but remember to add 1 to the FireSpread result (burned_ids),
# as FireSpread uses zero-indexing.
clip_landscape <- function(rc, steps, land) {
  row_lwr <- max(rc[1] - steps, 1)
  row_upr <- min(rc[1] + steps, land_rows)
  col_lwr <- max(rc[2] - steps, 1)
  col_upr <- min(rc[2] + steps, land_cols)

  # clip landscape
  land_clipped <- land[row_lwr:row_upr, col_lwr:col_upr, ]

  # define shifts for translation
  attr(land_clipped, "row_shift") <- row_lwr-1
  attr(land_clipped, "col_shift") <- col_lwr-1

  # shift ignition cell
  rc_new <- matrix(c(
    rc[1] - attr(land_clipped, "row_shift"),
    rc[2] - attr(land_clipped, "col_shift")
  ), ncol = 1)

  attr(land_clipped, "ig_rowcol") <- rc_new

  return(land_clipped)
}

# Function to simulate a single-ignition fire in the reduced landscape.
# coef: parameters vector, including steps.
# ig_rowcol: matrix[{row, col}, 1].
# land: spread landscape.
# returns the burned ids in the original (large) landscape.
simulate_one_fire <- function(coef, ig_rowcol, land) {
  ### testo
   # coef <- coef[1, ]
   # ig_rowcol <- ig_rowcols[, 1, drop = F]
  ###

  # prepare parameters
  coef_intercepts <- rep(coef[1], n_veg)
  coef_nd <- coef[2:3]
  coef_terrain <- coef[4:(n_coef-1)]
  steps <- floor(max(coef[n_coef], 1))

  # clip landscape, first updating the vegetation layer
  land_clipped <- clip_landscape(ig_rowcol[, 1], steps, land)

  # simulate fire
  fire_sim <- simulate_fire_compare(
    layer_vegetation = land_clipped[, , "veg"],
    layer_nd = land_clipped[, , nd_variables],
    layer_terrain = land_clipped[, , terrain_variables],
    coef_intercepts = coef_intercepts,
    coef_nd = coef_nd,
    coef_terrain = coef_terrain,
    ignition_cells = attr(land_clipped, "ig_rowcol") - 1, # zero-indexing
    upper_limit = upper_limit,
    steps = steps
  )

  # get burned cells row-col and translate to original landscape scale
  burned_ids <- fire_sim$burned_ids + 1 # result has zero-indexing
  burned_ids[1, ] <- burned_ids[1, ] + attr(land_clipped, "row_shift")
  burned_ids[2, ] <- burned_ids[2, ] + attr(land_clipped, "col_shift")

  return(burned_ids)
}

# Function to extract the FWI data for a given year
get_fwi <- function(year) {
  # Get FWI data
  fort_table_y <- fort_table[fort_table$year == year, ]

  fort_lwr <- fort_table_y$fort[1] - 10 # max lag used by fire models
  f_num_seq <- seq(max(fort_table_y$fort), fort_lwr, by = -1)
  fwi_local <- fwi_fort[[as.character(f_num_seq)]]

  # names(fwi_local)[1:26]  # focal year
  # names(fwi_local)[27:36] # before focal year
  return(fwi_local)
}

# Function to simulate one fire season (year, centred at summer).
# fwi_raster: FWI raster for PNNH, with the climate of the year to be simulated.
#   Its layers are the FWI for all pixels (8) in each fortnight, in decreasing
#   order. As FWI values are used lagged and accumulated, its 10 last layers
#   correspond to the 10 fortnights before the first one in the focal year.
# steps_int_shift: constant to add to the log-mean of the log-steps distribution.
#   This is aimed at correcting the fire size distribution, which is biased
#   towards large fires in the spred model.
# ig_rate_factor: factor to multiply the ignition rate lambdas. If the spread
#   model is expected to simulate the whole fire size distribution (including
#   one-pixel fires), the ignition rate has to be increased, by a factor of
#   1 / effective_ignition_rate. The effective ignition rate is the proportion
#   of ignitions that reach the simulator (CHECK this).

# It returns a list with the following elements:
#   1) ignitions_count;
#   2) table of fire_size and FWI (focal fortnight), only simulated fires, the
#      ones that escaped;
#   3) list of burned_ids by fire.

# To make the burn probability map, a large burned_layer will be created, and
# the list of burned_ids will be used to add ones to the burned pixels
simulate_fire_season <- function(fwi_raster,
                                 steps_int_shift = 0,
                                 ig_rate_factor = 1) {
  ## Placeholders to save results
  ig_count_h <- 0
  ig_count_l <- 0
  size_fwi_table <- NULL
  burned_ids_list <- NULL

  # standardize fwi values for ignition and escape models
  fwi_raster_z <- fwi_raster
  values(fwi_raster_z) <- (values(fwi_raster_z) - fwi_mean) / fwi_sd

  # Duplicate landscape, to record burned pixels as non-burnable in the veg layer
  pnnh_land_dyn <- pnnh_land

  # Duplicate cell ids from pnnh, to remove the ones already burned
  cells_pnnh_dyn <- cells_pnnh

  ## Loop over fortnights
  for(f in 1:nfy) {
    # Subset FWI for focal fortnight
    fbegin <- nfy - f + 1 # fortnights begin in the most recent to present
    fend <- fbegin + 10
    fwi_local <- fwi_raster[[fbegin:fend]]
    fwi_local_z <- fwi_raster_z[[fbegin:fend]]

    ## Ignition rate

    # sample parameters from ignition model
    iss <- sample(1:ipost, 1)
    # weight FWI temporally
    igw <- exp(-(tseq/imod[iss, "ls"])^2) |> normalize()
    fwi_pixels <- values(fwi_local_z) %*% igw

    # compute lambdas by pixel, and get weigthed sum
    lambdas_h <- exp(imod[iss, "a[1]"] + imod[iss, "b[1]"] * fwi_pixels)
    lambdas_l <- exp(imod[iss, "a[2]"] + imod[iss, "b[2]"] * fwi_pixels)
    lambda_h <- sum(lambdas_h * pix_weights) * ig_rate_factor
    lambda_l <- sum(lambdas_l * pix_weights) * ig_rate_factor

    n_ig_h <- rnegbin(1, lambda_h, imod[iss, "phi[1]"])
    n_ig_l <- rnegbin(1, lambda_l, imod[iss, "phi[2]"])

    # If there are no ignitions, jump to next fortnight
    if(n_ig_h + n_ig_l == 0) next

    candidate_cells <- sample(cells_pnnh_dyn, size = 1000 * 2,
                              replace = F)
    cand_vals <- pnnh_vals[candidate_cells, c("vfi", "tfi", "drz", "dhz")]

    # Sample location for human ignitions
    if(n_ig_h > 0) {
      ig_count_h <- ig_count_h + n_ig_h

      iprob_h <- exp(
        cand_vals[1:1000, "vfi"] * imod[iss, "c_vfi_h"] +
          cand_vals[1:1000, "tfi"] * imod[iss, "c_tfi_h"] +
          cand_vals[1:1000, "drz"] * imod[iss, "c_dist_r"] +
          cand_vals[1:1000, "dhz"] * imod[iss, "c_dist_h"]
      )
      # replace NA with 0
      iprob_h[is.na(iprob_h)] <- 0

      cells_ig_h <- sample(candidate_cells[1:1000], size = n_ig_h, replace = F,
                           prob = iprob_h)
    } else {
      cells_ig_h <- NULL
    }

    # Sample location for lightning ignitions
    if(n_ig_l > 0) {
      ig_count_l <- ig_count_l + n_ig_l

      iprob_l <- exp(
        cand_vals[1001:2000, "vfi"] * imod[iss, "c_vfi_l"] +
          cand_vals[1001:2000, "tfi"] * imod[iss, "c_tfi_l"]
      )

      # replace NA with 0
      iprob_l[is.na(iprob_l)] <- 0

      cells_ig_l <- sample(candidate_cells[1001:2000], size = n_ig_l, replace = F,
                           prob = iprob_l)
    } else {
      cells_ig_l <- NULL
    }

    ## Escape from ignited cells
    cells_ig <- c(cells_ig_h, cells_ig_l)


# BUG ACÁ: siempre estaba simulando una ig. SIEMPRE. ----------------------
#   (ESTO NO ESTABA COMENTADO)
    # ## test if only 1 cell
    # cells_ig <- cells_ig[1]; n_ig_l = 0; n_ig_h = 1
    # causes <- rep(c("human", "lightning"), c(n_ig_h, n_ig_l))

    vals_ig <- pnnh_vals[cells_ig, c("vfi", "tfi", "drz", "dhz"), drop = F]
    # get FWI at ignited cells
    crds_ig <- xyFromCell(pnnh_rast, cells_ig)
    fwi_ig_ts <- terra::extract(fwi_local_z, crds_ig) |> as.matrix()

    # cumulative FWI for escape model
    ess <- sample(1:epost, 1)
    escw <- exp(-(tseq/emod[ess, "ls"])^2) |> normalize()
    fwi_ig <- fwi_ig_ts %*% escw

    # Design matrix for escape probability
    Xesc <- cbind(vals_ig, fwi_ig)

    # Compute escape probability
    esc_betas <- emod[ess, c("b_vfi", "b_tfi", "b_drz", "b_dhz", "b_fwi")]
    escprobs <- 1 - plogis(emod[ess, "a[1]"] - Xesc %*% esc_betas)

    # Simulate escape
    escape_ids <- runif(length(cells_ig)) < escprobs
    cells_sim <- cells_ig[escape_ids]
    causes <- causes[escape_ids]
    nsimf <- length(cells_sim)

    if(nsimf == 0) next

    # Simulate spread from escaped cells

    # Translate cell id from PNNH raster to spread raster, which is larger
    crds_sim <- xyFromCell(pnnh_rast, cells_sim)
    cells_sim_land <- cellFromXY(pnnh_land_rast, crds_sim)
    ig_rowcols <- rowColFromCell(pnnh_land_rast, cells_sim_land) |> t()# 1-indexing

    # compute accumulated FWI for spread. (Use raw values, as spread uses another
    # mean and sd to standardize)
    fwi_spread_ts <- terra::extract(fwi_local, crds_sim) |> as.matrix()
    spreadw <- exp(-(tseq/ls_fwi_spread)^2) |> normalize()
    fwi_spread <- fwi_spread_ts %*% spreadw
    fwi_spread_z <- (fwi_spread - fwi_mean_spread) / fwi_sd_spread

    # Simulate spread parameters using fwi_spread_z
    Xspread <- cbind(rep(1, nrow(fwi_spread_z)), fwi_spread_z)
    sss <- sample(1:spost, 1)
    ab <- smod$fixef[1:n_coef, 1:2, sss]
    s <- smod$fixef[1:n_coef, 3, sss] |> sqrt()
    rho <- smod$rho[, , sss]
    V <- diag(s) %*% rho %*% diag(s)
    mu <- Xspread %*% t(ab)
    mu[, n_coef] <- mu[, n_coef] + steps_int_shift
    coefs_raw <- mgcv::rmvn(nrow(Xspread), mu, V)
    coefs <- constrain(coefs_raw, support)

    # simulate each fire separately, but randomizing order
    idsfire <- sample(1:nsimf, nsimf)
    for(j in idsfire) {
      # get burned ids
      if(coefs[j, n_coef] <= 1) {
        burned <- ig_rowcols[, j, drop = F]
      } else {
        burned <- simulate_one_fire(coefs[j, ],
                                    ig_rowcols[, j, drop = F],
                                    pnnh_land_dyn)
      }
      # complete size_fwi_table
      fire_size <- ncol(burned)
      mmat <- data.frame(size = fire_size,
                         steps = max(floor(coefs[j, n_coef]), 1), # avoid zeroes
                         steps_log_mu = mu[j, n_coef],
                         steps_log_sigma = s[n_coef],
                         fwi_focal = fwi_spread_ts[j, 1],
                         fwi_cum_spread = fwi_spread[j],
                         cause = causes[j])
      size_fwi_table <- rbind(size_fwi_table, mmat)

      # expand burned_ids_list
      burned_ids_list <- c(burned_ids_list, list(burned))

      # set burned cells as non-burnable in the vegetation layer of pnnh_land_dyn
      for(cell in 1:fire_size) {
        pnnh_land_dyn[burned[1, cell], burned[2, cell], "veg"] <- 99
      }

      # remove burned cells from candidate_cells in pnnh. First, translate cells
      # from large landscape to small (pnnh) one.
      cells_large <- cellFromRowCol(pnnh_land_rast,
                                    row = burned[1, ],
                                    col = burned[2, ])
      crds_ <- xyFromCell(pnnh_land_rast, cells_large)
      cells_small <- cellFromXY(pnnh_rast, crds_)
      # include the ones burned by ignitions!
      cells_burned <- c(cells_small, cells_ig) |> unique()

      keep <- !(cells_pnnh_dyn %in% cells_burned)
      cells_pnnh_dyn <- cells_pnnh_dyn[keep]
    } # end loop over escaped cells
  } # end loop over fortnights

  # save results in a list
  out <- list(
    ig_count = c("human" = ig_count_h, "lightning" = ig_count_l),
    size_fwi_table = size_fwi_table,
    burned_ids_list = burned_ids_list,
    steps_int_shift = steps_int_shift,
    ig_rate_factor = ig_rate_factor
  )
  # reminder: burned_ids correspond to the large landscape.
  rm(pnnh_land_dyn, cells_pnnh_dyn)
  gc()
  return(out)
}

# the same function, but simulating the year inside
simulate_fire_year <- function(year,
                               steps_int_shift = 0,
                               ig_rate_factor = 1) {
  simulate_fire_season(get_fwi(year), steps_int_shift, ig_rate_factor)
}

# Simulate in parallel
simulate_fire_year_parallel <- function(nsim = 100, cores = 2,
                                        steps_int_shift = 0,
                                        ig_rate_factor = 1) {
  registerDoMC(cores)
  years <- sample(1999:2022, nsim, T) |> as.list()
  # simulate in parallel
  result <- foreach(yy = years) %dopar% {
    simulate_fire_year(yy, steps_int_shift, ig_rate_factor)
  }
  return(result)
}

get_sizeclass <- function(x, cuts = c(0.09, 10, 100)) {
  out <- sapply(cuts, function(cut) {
    as.numeric(x > cut)
  }) |> rowSums() + 1
  return(out)
}

# Get fire size-class distribution
size_dist <- function(x, cuts = c(0.09, 10, 100)) {
  sizeclass <- get_sizeclass(x, cuts)
  tt <- table(sizeclass)
  tt <- tt / sum(tt)
  names(tt) <- c("(0, 0.09]", "(0.09, 10]", "(10, 100]", "(100, ...)")
  return(tt)
}

# Emulator of (ordinal) fire-size distribution, considering corrections for
# steps parameters: a_shift = steps_int_shift, s_factor = multiplier of steps
# sigma.
sim_sizeclass <- function(a_shift = -1, s_factor = 0.5, nsim = 10000) {
  # a_shift = -1; s_factor = 0.5; nsim = 10000

  # resample FWI for spread
  rows <- sample(1:nrow(size_table), size = nsim, replace = T)
  fwi_ <- (size_table$fwi_cum_spread[rows] - fwi_mean_spread) / fwi_sd_spread

  # get intercepts
  ids <- sample(1:spost, size = nsim, replace = ifelse(nsim>spost, T, F))
  abs <- smod$fixef[n_coef, c("a", "b"), ids]
  mus <- smod$fixef[n_coef, "a", ids] + smod$fixef[n_coef, "b", ids] * fwi_
  sigmas <- smod$fixef[n_coef, "s2", ids] |> sqrt()

  # simulated modified steps
  steps_sim <- rnorm(nsim, mus + a_shift, sigmas * s_factor)

  # predict fire size distribution
  pp <- predict(sizeclass_model, newdata = data.frame(steps_log = steps_sim),
                type = "response")
  return(colMeans(pp))
}

# function of the parameters to compare distributions
size_distance <- function(x) {
  dsim <- sim_sizeclass(x[1], x[2])
  dd <- sum((dsim - ref_dist) ^ 2)
}

# Bhattacharyya distance between empirical densities, computed from samples x and y.
# n is the number of points to compute the initial density, then it is
# multiplied by 3. lower is the lower limit to fit the density.
bhattacharyya_distance <- function(x, y, n = 2 ^ 11, lower = NULL) {
  # Initial densities, to check limits
  if(is.null(lower)) {
    dx0 <- density(x, n = n)
    dy0 <- density(y, n = n)
  } else {
    dx0 <- density(x, from = lower, n = n)
    dy0 <- density(y, from = lower, n = n)
  }

  if(is.null(lower)) {
    lower <- min(c(dx0$x, dy0$x))
  }
  upper <- max(c(dx0$x, dy0$x))
  # refit densities with lower and upper
  dx <- density(x, from = lower, to = upper, n = n * 3)
  dy <- density(y, from = lower, to = upper, n = n * 3)

  sqrt_prod <- function(x) {
    dens_x <- approx(dx$x, dx$y, xout = x, rule = 2)$y
    dens_y <- approx(dy$x, dy$y, xout = x, rule = 2)$y
    return(sqrt(dens_x * dens_y))
  }

  bc <- integrate(sqrt_prod, lower, upper, subdivisions = 1e5)$value
  return(-log(bc))
}

bhattacharyya_distance_cat <- function(x, y) {
  -log(sum(sqrt(x * y)))
}


# Constants for spread ----------------------------------------------------

# constants for fire spread simulation
upper_limit <- 1
n_veg <- 5
veg_names <- c("wet", "subalpine", "dry", "shrubland", "grassland")
veg_levels <- c("Wet forest", "Subalpine forest", "Dry forest", "Shrubland", "Grassland")

n_terrain <- 2
terrain_names <- c("slope", "wind")
terrain_variables <- c("elevation", "wdir", "wspeed")
n_nd <- n_fi <- 2        # flammability indices
nd_variables <- c("vfi", "tfi")

par_names <- c("intercept", nd_variables, terrain_names, "steps")
n_coef <- length(par_names)

par_names_all <- c(par_names, "area")

n_par <- length(par_names_all)
n_pt <- 3 # b0, b1, s2

# support for parameters
ext_alpha <- 50
ext_beta <- 30

slope_sd <- fi_params$slope_term_sd

params_lower <- c(-ext_alpha, rep(0, n_coef-2), 5)
params_upper <- c(ext_alpha, rep(ext_beta, n_coef-2), NA)
names(params_lower) <- names(params_upper) <- par_names
params_upper["slope"] <- ext_beta / slope_sd

support <- rbind(params_lower, params_upper)
colnames(support) <- names(params_lower) <- names(params_upper) <- par_names
support_width <- apply(support, 2, diff)

fwi_spread_mean_sd <- readRDS(
  file.path("files", "hierarchical_model", "fwi_mean_sd_spread.rds")
)

fwi_mean_spread <- fwi_spread_mean_sd$fwi_mean
fwi_sd_spread <- fwi_spread_mean_sd$fwi_sd

# FWI lengthscale (cumulative effect)
ls_fwi_spread <- 2.751998
# fortnight_mod$ls_hat # 2.751998 fortnights
# from <weather temporal scale/FWI temporal scale.R>, line ~ 576

# Data -------------------------------------------------------------------

# FWI
fwi_fort <- rast(file.path("data", "fwi_daily_1998-2022",
                           "fwi_fortnights_19980101_20230630.tif"))
fwi_fort <- project(fwi_fort, "EPSG:5343")

# Nahuel Huapi National Park
pnnh <- vect(file.path("data", "protected_areas", "apn_limites.shp"))
pnnh <- pnnh[pnnh$nombre == "Nahuel Huapi", ]
pnnh <- project(pnnh, "EPSG:5343")

# Nahuel Huapi raster (with distance to humans, to simulate ignitions)
# pnnh_rast <- rast(file.path("data", "pnnh_images",
#                             "pnnh_data_30m.tif"))
## Load file with standardized variables (code in prepare PNNH raster)

# Nahuel Huapi raster (without distance to humans!)
pnnh_land_rast <- rast(file.path("data", "pnnh_images",
                                 "pnnh_data_spread_buffered_30m_smaller.tif"))

# Spread landscape
pnnh_land <- readRDS(file.path("data", "pnnh_images",
                               "pnnh_spread_landscape_smaller.rds"))

land_rows <- dim(pnnh_land)[1]
land_cols <- dim(pnnh_land)[2]

# ncell(pnnh_land_rast) == prod(dim(pnnh_land)[1:2]) # OK

# metrics to standardize predictors
pnnh_data_summary <- readRDS(file.path("data", "pnnh_images",
                                       "pnnh_data_summary.rds"))

dr_mean <- pnnh_data_summary$dr_mean  # distance from roads
dr_sd <- pnnh_data_summary$dr_sd
dh_mean <- pnnh_data_summary$dh_mean  # distance from human settlements
dh_sd <- pnnh_data_summary$dh_sd
fwi_mean <- pnnh_data_summary$fwi_mean
fwi_sd <- pnnh_data_summary$fwi_sd

## Observed fires data, in PNNH, to compare size distribution with simulations
igdata <- read.csv(file.path("data", "ignition", "ignition_size_data.csv"))
igdata$area_impute2 <- igdata$area_impute
igdata$area_impute2[igdata$area_impute < 0.09] <- 0.09

# Import vegetation class transforms --------------------------------------

dveg <- readxl::read_excel("/home/ivan/Insync/Mapa vegetación WWF - Lara et al. 1999/clases de vegetacion y equivalencias.xlsx",
                           sheet = "Sheet2")
# Make urban as forest, and get zero-indexing
dveg$cnum2[dveg$class1 == "Urban"] <- 1
dveg$class2[dveg$class1 == "Urban"] <- "Wet forest"
# urban is taken as forest so its burn probability changes markedly with NDVI.
n_veg_types <- V <- 5

# rename to use left_join
dveg$veg_focal <- dveg$cnum1
dveg$veg_num <- dveg$cnum2


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

# Import fire models -------------------------------------------------------

igmod <- readRDS(file.path("files", "ignition", "ignition_model_samples.rds"))
escmod <- readRDS(file.path("files", "ignition", "sizeclass_model_samples.rds"))
smod <- readRDS(file.path("files", "hierarchical_model",
                          "spread_model_samples.rds"))

# extract parameters from stan models
imod <- as.matrix(
  igmod,
  pars = c("a", "b", "phi", "ls",
           "c_fi", "c_dist", "human_prop",
           "c_vfi_h", "c_vfi_l",
           "c_tfi_h", "c_tfi_l",
           "c_dist_h", "c_dist_r")
)

emod <- as.matrix(
  escmod,
  pars = c("a", "b_fwi", "b_vfi", "b_tfi", "b_drz", "b_dhz", "ls")
)

# number of posterior samples
ipost <- nrow(imod)
epost <- nrow(imod)
spost <- dim(smod$fixef)[3]

# Prepare PNNH raster -----------------------------------------------------

# pnnh_rast$vegetation <- subst(pnnh_rast[["veg"]], dveg$cnum1, dveg$cnum2) # cnum2 has 1:5
# pnnh_rast$vfi <- vfi_calc(values(pnnh_rast$vegetation),
#                           values(pnnh_rast$ndvi))
# pnnh_rast$tfi <- tfi_calc(values(pnnh_rast$elevation),
#                           values(pnnh_rast$aspect),
#                           values(pnnh_rast$slope))
# pnnh_rast$drz <- (pnnh_rast$dist_roads / 1000 - dr_mean) / dr_sd
# pnnh_rast$dhz <- (pnnh_rast$dist_humans / 1000 - dh_mean) / dh_sd
# writeRaster(pnnh_rast,
#             file.path("data", "pnnh_images", "pnnh_data_30m_std.tif"))
pnnh_rast <- rast(file.path("data", "pnnh_images", "pnnh_data_30m_std.tif"))

# Define ignitable cells (inside park and burnable)
pnnh_rast <- mask(pnnh_rast, pnnh)
cells_pnnh0 <- cells(pnnh_rast) # all cells inside park
pnnh_vals <- values(pnnh_rast)
pnnh_veg <- pnnh_vals[cells_pnnh0, "vegetation"] # get values for cells inside park
cells_pnnh <- cells_pnnh0[!is.na(pnnh_veg)] # burnable (ignitable) cells

# Simulation settings -----------------------------------------------------

years <- 1999:2022
fort_table <- dtable[!duplicated(dtable[, c("fort", "fort_focal", "year")]), ]
nforty <- max(fort_table$fort_focal)
tseq <- 0:10 # temporal lagg for FWI (from focal to -10 fortnight)
export_dir <- file.path("files", "fire_regime_simulation")

# Simulate 10000 years ---------------------------------------------------

# mb6 <- microbenchmark::microbenchmark(
#   simulate_fire_year_parallel(1000, 6),
#   times = 1
# )
## 366.8561 s / 1000 sim
## 10000 * 366.8561 / 1000 / (60 * 60) = 1.01 h / 10000 iter
## ~ 6 min / 1000 sim
## more cores might compromise the RAM

nsim <- 10000
batch_size <- 1000
nbatch <- nsim / batch_size

# for(i in 1:nbatch) {
#   print(i)
#   out <- simulate_fire_year_parallel(batch_size, 6)
#   num <- str_pad(i, 2, pad = "0")
#   nn <- paste("sim_", num, ".rds", sep = "")
#   saveRDS(out, file.path(export_dir, nn))
#   rm(out); gc()
# }

ff <- list.files(export_dir, pattern = "sim_")

# load simulations
sims <- do.call("c", lapply(1:length(ff), function(i) {
  readRDS(file.path(export_dir, ff[i]))
}))

# get fire size table
size_table <- do.call("rbind", lapply(sims, function(x) {
  x$size_fwi_table
}))
# nro de fires en 10000 años:
nrow(size_table)
# 13982

sort(size_table$area_ha, decreasing = T)

barplot(size_dist(size_table$area_ha))
plot(hist(log10(size_table$area_ha[size_table$area_ha > 10])))
table(size_table$area_ha)
str(sims[[2]])


# how does it change if area is multiplied by a factor?
barplot(size_dist(size_table$area_ha * 0.005))



# Simulate 5000 years recording more variables ---------------------------

# # assess required sample size
# yy <- log10(size_table$area_ha[size_table$area_ha > 10])
# prop = length(yy) / nrow(size_table)
# par(mfrow = c(1, 2))
# hist(yy, main = "all", xlab = "log10 ha")
# ids <- sample(1:length(yy), floor(5000 * prop), replace = F)
# hist(yy[ids], main = "subsample", xlab = "log10 ha")
# par(mfrow = c(1, 1))
## 5000 may be OK

# To better understand the steps-size relationship
nsim <- 5000
batch_size <- 1000
nbatch <- nsim / batch_size

# for(i in 1:nbatch) {
#   print(i)
#   out <- simulate_fire_year_parallel(batch_size, 6)
#   num <- str_pad(i, 2, pad = "0")
#   nn <- paste("sim_02_", num, ".rds", sep = "")
#   saveRDS(out, file.path(export_dir, nn))
#   rm(out); gc()
# }

ff <- list.files(export_dir, pattern = "sim_02_")

# load simulations
sims <- do.call("c", lapply(1:length(ff), function(i) {
  readRDS(file.path(export_dir, ff[i]))
}))

ig_counts <- do.call("rbind", lapply(sims, function(x) {
  x$ig_count
}))


# get fire size table
size_table <- do.call("rbind", lapply(sims, function(x) {
  x$size_fwi_table
}))
# nro de fires en 10000 años:
nrow(size_table)
# 7039

size_table$area_ha <- size_table$size * 0.09 # pixels to ha
barplot(size_dist(size_table$area_ha))
barplot(size_dist(igdata$area_impute))

# ignitions that were simulated:
nrow(size_table) / sum(ig_counts) * 100 # 10.72 %
# (this is the simulated escape proportion)

qqplot(log10(igdata$area_impute), log10(size_table$area_ha))
abline(0, 1)

# ignitions that effectively spread more than one pixels:
nrow(size_table[size_table$size > 1, ]) / sum(ig_counts) * 100 # 8.41 %
# (this is the effective simulated escape proportion, i.e., considering
# the spread model sometimes makes smaller fires.)

# probability of the spread model burning more than one pixel
sum(size_table$size > 1) / nrow(size_table) * 100 # 78.43 %

##

# Fit steps parameters to match the PNNH fire-size dist -------------------

# Tune at least the intercept and sigma of the steps (spread) model (~ FWI)
# to match the fire size distribution at PNNH.
# To do so, we will fit an emulator of simulated_fire_size = f(steps).
# Then, by simulation, we will get the correction-parameters that minimize the
# KL divergence between observed and emulated fire size dist.


### sim_size ~ steps
ggplot(size_table, aes(log10(steps), log10(area_ha))) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = c(log10(0.09), log10(10)),
             linetype = "dashed", color = "red")

ggplot(size_table, aes(log10(steps), log10(size))) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = c(log10(1), log10(2), log10(10 / 0.09)),
             linetype = "dashed", color = "red") +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "blue") +
geom_abline(intercept = 0, slope = 2.5,
            linetype = "dashed", color = "green")
# there is a maximum size, which depends almost linearly of the steps

ggplot(size_table, aes(log(steps), log(size))) +
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = c(log(1), log(2), log(10 / 0.09)),
             linetype = "dashed", color = "red") +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "blue") +
  geom_abline(intercept = 0, slope = 2.5,
              linetype = "dashed", color = "green") +

  geom_point(aes(log))


ggplot(size_table, aes(steps, area_ha)) +
  geom_point(alpha = 0.4)


# La distribución de size_sim ~ steps es reee compleja. Está difícil
# aproximarla bien. Veamos primero con una aproximación ordinal, a ver si
# podemos lograr algo.


# Categorical fire size emulator -----------------------------------------
# Ordinal was not good

library(mgcv)
library(DHARMa)

size_table$sizeclass <- get_sizeclass(size_table$area_ha)
size_table$sizeclass_gam <- get_sizeclass(size_table$area_ha) - 1
size_table$steps[size_table$steps == 0] <- 1
size_table$steps_log <- log(size_table$steps)

catmod <- gam(list(sizeclass_gam
                   ~ s(steps_log, k = 6, bs = "cr"),
                   ~ s(steps_log, k = 6, bs = "cr"),
                   ~ s(steps_log, k = 6, bs = "cr")),
              family = multinom(K = 4-1), data = size_table)

ordmod <- gam(sizeclass ~ s(steps_log, k = 20, bs = "cr"),
              family = ocat(R = 4), data = size_table)

sizeclass_model <- ordmod

fitprobs <- predict(sizeclass_model, type = "response")
ysim <- sapply(1:500, function(i) {
  apply(fitprobs, 1, function(x) sample(1:4, 1, prob = x))
})

res <- createDHARMa(simulatedResponse = ysim,
                    observedResponse = size_table$sizeclass)
# plot(res)
plotResiduals(res, form = log(size_table$steps), rank = F)

barplot(table(as.vector(ysim)))
barplot(table(size_table$sizeclass))
barplot(colMeans(fitprobs))

# Optimize steps intercept shift using fire-size distribution emulator --------

ref_dist <- size_dist(igdata$area_impute)

# Compute difference in a gridpar
pargrid <- expand.grid(a_shift = seq(-8, 0.01, length.out = 50),
                       s_factor = seq(0.01, 1, length.out = 20))
parmat <- as.matrix(pargrid)
pargrid$loss <- sapply(1:nrow(pargrid), function(i) {
  if(i %% 100 == 0) print(i)
  size_distance(parmat[i, ])
})

plot(pargrid$loss ~ pargrid$a_shift)
plot(pargrid$loss ~ pargrid$s_factor)

# approximate with GAM
lossgam <- gam(loss ~ s(a_shift, k = 15, bs = "cr"),
               family = Gamma(link = "log"),
               data = pargrid)
opfun <- function(x) {
  predict(lossgam, data.frame(a_shift = x), type = "response")
}
opt_a_shift <- optim(pargrid$a_shift[which.min(pargrid$loss)], opfun,
                     method = "Brent", lower = -7, upper = 0)

steps_int_shift <- opt_a_shift$par # -4.512744

size_best <- sim_sizeclass(a_shift = steps_int_shift,
                           s_factor = 1)

par(mfrow = c(1, 2))
barplot(size_best, main = "fitted")
barplot(ref_dist, main = "observed")
par(mfrow = c(1, 1))


# Kernel regression emulator ----------------------------------------------
# HARD TO SIMULATE FROM HERE
# library(hdrcde)
# size_table$steps_log <- log(size_table$steps)
# size_table$area_log <- log(size_table$size * 0.09)
# kde <- cde(size_table$steps_log, size_table$area_log, degree = 2,
#            x.name = "steps", y.name = "area", nxmargin = 100,
#            x.margin = seq(0, 12, length.out = 100),
#            b = 1, a = 0.5)
# plot(kde)
# kde$z


# Simulate fires with emulator-optimized shift ----------------------------

# To better understand the steps-size relationship
nsim <- 5000
batch_size <- 1000
nbatch <- nsim / batch_size

steps_int_shift <- -4.512744 # opt_a_shift$par

# for(i in 1:nbatch) {
#   print(i)
#   out <- simulate_fire_year_parallel(batch_size, 6,
#                                      steps_int_shift = steps_int_shift)
#   num <- str_pad(i, 2, pad = "0")
#   nn <- paste("sim_03_", num, ".rds", sep = "")
#   saveRDS(out, file.path(export_dir, nn))
#   rm(out); gc()
# }

ff <- list.files(export_dir, pattern = "sim_03_")

# load simulations
sims <- do.call("c", lapply(1:length(ff), function(i) {
  readRDS(file.path(export_dir, ff[i]))
}))

ig_counts <- do.call("rbind", lapply(sims, function(x) {
  x$ig_count
}))

# get fire size table
size_table <- do.call("rbind", lapply(sims, function(x) {
  x$size_fwi_table
}))
# nro de fires en 10000 años:
nrow(size_table)
# 7071

size_table$area_ha <- size_table$size * 0.09 # pixels to ha
barplot(size_dist(size_table$area_ha))
barplot(size_dist(igdata$area_impute))

# ignitions that were simulated:
nrow(size_table) / sum(ig_counts) * 100 # 10.75 %
# (this is the simulated escape proportion)

probs <- ppoints(floor(nrow(igdata) / 2))
qobs <- quantile(log10(igdata$area_impute2), probs = probs, method = 8)
qsim <- quantile(log10(size_table$area_ha), probs = probs, method = 8)
plot(qsim ~ qobs)
abline(0, 1)


ks <- ks.test(log10(size_table$area_ha), log10(igdata$area_impute2),
              alternative = "less")
ks <- ks.test(log10(size_table$area_ha), log10(igdata$area_impute2),
              alternative = "greater")

ks <- ks.test(log10(igdata$area_impute2), log10(size_table$area_ha),
              alternative = "two.sided")
ks <- ks.test(log10(igdata$area_impute2), log10(igdata$area_impute2),
              alternative = "two.sided")
ks$statistic


# Parece que la KL es mejor (usar FNN)
library(FNN)
y <- log10(size_table$area_ha)
ref <- log10(igdata$area_impute2)
KL.divergence(y, ref, k = 5, algorithm = "kd_tree")



plot(density(log10(size_table$area_ha), from = log10(0.09)),
     xlim = c(log10(0.09), 4))
lines(density(log10(igdata$area_impute2), from = log10(0.09)), col = 2)

# ignitions that effectively spread more than one pixels:
nrow(size_table[size_table$size > 1, ]) / sum(ig_counts) * 100 # 1.21 %
# (this is the effective simulated escape proportion, i.e., considering
# the spread model sometimes makes smaller fires.)

# probability of the spread model burning more than one pixel
sum(size_table$size > 1) / nrow(size_table) * 100 # 11.27 %




# Search over intercept-shift space ---------------------------------------

# 10 values between -4.512744 and 0 (these were aready evaluated, at files
# sim_03_ and before, respectively).

shiftvals <- c(
  seq(-4.512744, -4.512744/2, length.out = 9)[2:8],
  seq(-4.512744/2, 0, length.out = 4)[-4]
)
# plot(c(-4.512744, 0)~ c(0, 11))
# points(shiftvals, col = 2)

ns <- length(shiftvals)

# not running in batches, all at once.
for(i in 1:ns) {
  # i = 1
  print(i)
  out <- simulate_fire_year_parallel(5000, 6,
                                     steps_int_shift = shiftvals[i])
  num <- str_pad(i+3, 2, pad = "0") # starting at file id number 4
  nn <- paste("sim_", num, "_00.rds", sep = "")
  saveRDS(out, file.path(export_dir, nn))
  rm(out); gc()
}

# TAREA: Calcular la distancia en distribución discreta y en continua,
# usando la bhattacharyya distance.



# Compare fit varying steps_intercerpt_shift ------------------------------

## Load all full simulations

# Load shift = 0
ff <- list.files(export_dir, pattern = "sim_02_")
sims0 <- do.call("c", lapply(1:length(ff), function(i) {
  readRDS(file.path(export_dir, ff[i]))
}))
# get fire size table
size_table_12 <- do.call("rbind", lapply(sims0, function(x) {
  x$size_fwi_table
}))

# Load shift = -4.512744
ff <- list.files(export_dir, pattern = "sim_03_")
# load simulations
sims_minshift <- do.call("c", lapply(1:length(ff), function(i) {
  readRDS(file.path(export_dir, ff[i]))
}))
# get fire size table
size_table_1 <- do.call("rbind", lapply(sims_minshift, function(x) {
  x$size_fwi_table
}))

# All (12)
shiftvals <- c(
  seq(-4.512744, -4.512744/2, length.out = 9)[1:8],
  seq(-4.512744/2, 0, length.out = 4)
)
ns <- length(shiftvals)

stable_list <- vector("list", ns)
stable_list[[1]] <- size_table_1   # starts from most negative shift
stable_list[[ns]] <- size_table_12

for(i in 2:(ns-1)) {
  num <- str_pad(i+2, 2, pad = "0") # starting at file id number 4
  nn <- paste("sim_", num, "_00.rds", sep = "")
  x <- readRDS(file.path(export_dir, nn))
  stable_list[[i]] <- do.call("rbind", lapply(x, function(x) {
    x$size_fwi_table
  }))
}





## Compare distributions considering all fires ----------------------------

comptab <- data.frame(
  shift = shiftvals,
  b_cat = NA,
  b_cont = NA,
  sumq = NA
)

# reference distribution
ref_cat <- size_dist(igdata$area_impute2)
ref_cont <- log10(igdata$area_impute2)
lower <- log10(0.09)

# Compute metrics
for(i in 1:ns) {
  dist_cat <- size_dist(stable_list[[i]]$size * 0.09) # pixels to ha
  dist_cont <- log10(stable_list[[i]]$size * 0.09)

  comptab$b_cat[i] <- bhattacharyya_distance_cat(dist_cat, ref_cat)
  comptab$b_cont[i] <- bhattacharyya_distance(dist_cont, ref_cont,
                                                     lower = lower)
  comptab$sumq[i] <- sum(10 ^ dist_cont) / sum(10 ^ ref_cont)
}

par(mfrow = c(2, 2))
plot(b_cat ~ shift, dat = comptab, pch = 19,
     main = "categorical size distribution")
points(b_cat ~ shift, dat = comptab[which.min(comptab$b_cat), ],
       pch = 19, col = "blue")

plot(b_cont ~ shift, dat = comptab, pch = 19,
     main = "continuous size distribution")
points(b_cont ~ shift, dat = comptab[which.min(comptab$b_cont), ],
       pch = 19, col = "blue")

plot(log(sumq) ~ shift, dat = comptab, pch = 19,
     main = "sum quotiente (sim/obs")
abline(h = log(1))
points(log(sumq) ~ shift, dat = comptab[which.min(abs(comptab$sumq - 1)), ],
       pch = 19, col = "blue")

plot(sumq ~ shift, dat = comptab, pch = 19,
     main = "sum quotiente (sim/obs")
abline(h = 1)
points(sumq ~ shift, dat = comptab[which.min(abs(comptab$sumq - 1)), ],
       pch = 19, col = "blue")
par(mfrow = c(1, 1))


# QQplots
par(mfrow = c(3, 4))
for(i in 1:ns) {
  # i = 1
  dist_cont <- log10(stable_list[[i]]$size * 0.09)
  qqplot(ref_cont, dist_cont, xlab = "observed", ylab = "simulated",
         main = comptab$shift[i])
  abline(0, 1)
}
par(mfrow = c(1, 1))

# Densities
par(mfrow = c(3, 4))
for(i in 1:ns) {
  # i = 1
  dist_cont <- log10(stable_list[[i]]$size * 0.09)
  plot(density(dist_cont, from = lower, n = 2 ^ 11), col = "blue",
       main = comptab$shift[i], xlab = "fire size (log10 ha)")
  lines(density(ref_cont, from = lower, n = 2 ^ 11), col = "black")
}
par(mfrow = c(1, 1))


# Barplots
par(mfrow = c(3, 4))
for(i in 1:ns) {
  # i = 1
  dist_cat <- size_dist(stable_list[[i]]$size * 0.09) # pixels to ha
  barplot(dist_cat, col = "blue",
          main = comptab$shift[i], xlab = "fire size class")
  barplot(ref_cat, col = rgb(0, 0, 0, 0.4), add = T)
}
par(mfrow = c(1, 1))


## Compare distributions considering fires > 0.09 ha ---------------------

comptab <- data.frame(
  shift = shiftvals,
  b_cat = NA,
  b_cont = NA,
  meanq = NA,
  sdq = NA,
  meanq_l1 = NA,
  sdq_l1 = NA
)

# reference distribution
ref_cat <- size_dist(igdata$area_impute2)[-1] |> normalize()
ref_cont <- log10(igdata$area_impute2[igdata$area_impute2 > 0.09])
lower <- log10(0.09)

# Compute metrics
for(i in 1:ns) {
  dist_cat <- size_dist(stable_list[[i]]$size * 0.09)[-1] |> normalize()
  tmp <- stable_list[[i]]$size * 0.09
  dist_cont <- log10(tmp[tmp > 0.09])

  comptab$b_cat[i] <- bhattacharyya_distance_cat(dist_cat, ref_cat)
  comptab$b_cont[i] <- bhattacharyya_distance(dist_cont, ref_cont,
                                              lower = lower)
  comptab$meanq[i] <- mean(10 ^ dist_cont) / mean(10 ^ ref_cont)
  comptab$sdq[i] <- sd(10 ^ dist_cont) / sd(10 ^ ref_cont)

  # the same metrics but removing the single largest
  dist2 <- sort(dist_cont, decreasing = T)[-1]
  ref2 <- sort(ref_cont, decreasing = T)[-1]
  comptab$meanq_l1[i] <- mean(10 ^ dist2) / mean(10 ^ ref2)
  comptab$sdq_l1[i] <- sd(10 ^ dist2) / sd(10 ^ ref2)
}

par(mfrow = c(2, 2))
plot(b_cat ~ shift, dat = comptab, pch = 19,
     main = "categorical size distribution")
points(b_cat ~ shift, dat = comptab[which.min(comptab$b_cat), ],
       pch = 19, col = "blue")

plot(b_cont ~ shift, dat = comptab, pch = 19,
     main = "continuous size distribution")
points(b_cont ~ shift, dat = comptab[which.min(comptab$b_cont), ],
       pch = 19, col = "blue")

plot(meanq ~ shift, dat = comptab, pch = 19,
     main = "means quotient (sim/obs)")
abline(h = 1)
points(meanq ~ shift, dat = comptab[which.min(abs(comptab$meanq - 1)), ],
       pch = 19, col = "blue")

# plot(log(meanq) ~ shift, dat = comptab, pch = 19,
#      main = "means quotient (sim/obs)")
# abline(h = log(1))
# points(log(meanq) ~ shift, dat = comptab[which.min(abs(comptab$meanq - 1)), ],
#        pch = 19, col = "blue")

# plot(sdq ~ shift, dat = comptab, pch = 19,
#      main = "sd quotient (sim/obs), log")
# abline(h = 1)
# points(sdq ~ shift, dat = comptab[which.min(abs(comptab$sdq - 1)), ],
#        pch = 19, col = "blue")

plot(meanq_l1 ~ shift, dat = comptab, pch = 19,
     main = "means quotient (sim/obs)")
abline(h = 1)
points(meanq_l1 ~ shift, dat = comptab[which.min(abs(comptab$meanq_l1 - 1)), ],
       pch = 19, col = "blue")
par(mfrow = c(1, 1))

# QQplots
par(mfrow = c(3, 4))
for(i in 1:ns) {
  # i = 1
  tmp <- stable_list[[i]]$size * 0.09
  dist_cont <- log10(tmp[tmp > 0.09])
  qqplot(ref_cont, dist_cont, xlab = "observed", ylab = "simulated",
         main = comptab$shift[i])
  abline(0, 1)
}
par(mfrow = c(1, 1))

# Densities
par(mfrow = c(3, 4))
for(i in 1:ns) {
  # i = 1
  tmp <- stable_list[[i]]$size * 0.09
  dist_cont <- log10(tmp[tmp > 0.09])

  dsim <- density(dist_cont, from = lower, n = 2 ^ 11)
  dobs <- density(ref_cont, from = lower, n = 2 ^ 11)
  yy <- range(c(dsim$y, dobs$y)) * 1.05
  plot(dsim, col = "blue", ylim = yy,
       main = comptab$shift[i], xlab = "fire size (log10 ha)")
  lines(dobs, col = "black")
}
par(mfrow = c(1, 1))

# Densities notlog
par(mfrow = c(3, 4))
for(i in 1:ns) {
  # i = 9
  tmp <- stable_list[[i]]$size * 0.09
  dist_cont <- tmp[tmp > 0.09]
  sort(dist_cont, decreasing = T) |> summary()
  sort(dist_cont, decreasing = T)[-1] |> summary()

  dsim <- density(dist_cont, from = lower, n = 2 ^ 11)
  dobs <- density(10 ^ ref_cont, from = lower, n = 2 ^ 11)

  yy <- range(c(dsim$y, dobs$y)) * 1.05
  xx <- range(c(dsim$x, dobs$x))

  plot(dsim, col = "blue", ylim = yy, xlim = xx,
       main = comptab$shift[i], xlab = "fire size (ha)")
  lines(dobs, col = "black")
}
par(mfrow = c(1, 1))


# Barplots
par(mfrow = c(3, 4))
for(i in 1:ns) {
  # i = 1
  dist_cat <- size_dist(stable_list[[i]]$size * 0.09)[-1] |> normalize()# pixels to ha
  barplot(dist_cat, col = "blue",
          main = comptab$shift[i], xlab = "fire size class")
  barplot(ref_cat, col = rgb(0, 0, 0, 0.4), add = T)
}
par(mfrow = c(1, 1))



# Ideas: ------------------------------------------------------------------

# Comparar la distribución únicamente usando los fuegos > 0.09 ha (1 pix).
# Opciones:
#   (1) Tunear el intercept sólo en base a eso, dejando como no quemado los
#       píxeles afectados por ese tipo de ignición.
#       Luego calcular la tasa efectiva de escape, como la proporción de incendios
#       simulados que salen de 1 pix. Para corregir el nro de incendios,
#       simular
#   (2) Simular steps de una normal truncada en log(2), para que nunca simule
#       menos de un pixel.
#       Pero esto no lo va a resolver. Demasiadas igniciones ocurren con
#       steps altos que de todos modos no superan el pixel.



# Calcular burn probability en el paisaje ---------------------------------

# Import full simulations

# Load shift = 0
ff <- list.files(export_dir, pattern = "sim_02_")
sims0 <- do.call("c", lapply(1:length(ff), function(i) {
  readRDS(file.path(export_dir, ff[i]))
}))

# Load shift = -4.512744
ff <- list.files(export_dir, pattern = "sim_03_")
sims_minshift <- do.call("c", lapply(1:length(ff), function(i) {
  readRDS(file.path(export_dir, ff[i]))
}))

# All (12)
shiftvals <- c(
  seq(-4.512744, -4.512744/2, length.out = 9)[1:8],
  seq(-4.512744/2, 0, length.out = 4)
)
ns <- length(shiftvals)

simlist <- vector("list", ns)
simlist[[1]] <- sims_minshift   # starts from most negative shift
simlist[[ns]] <- sims0          # zero shift

for(i in 2:(ns-1)) {
  num <- str_pad(i+2, 2, pad = "0") # starting at file id number 4
  nn <- paste("sim_", num, "_00.rds", sep = "")
  x <- readRDS(file.path(export_dir, nn))
  simlist[[i]] <- x
}

nsim <- length(simlist[[1]])
rprob_list <- vector("list", ns)
prob_means <- numeric(ns)
for(i in 1:ns) {
  i = 12
  print(i)
  bp <- matrix(0, nrow(pnnh_land), ncol(pnnh_land))

  # set non-burnable as -1
  bp[pnnh_land[, , "veg"] == 99] <- NA

  # add 1 in burned pixels
  for(y in 1:nsim) {
    print(y)
    # i = 1; y = 1
    nf <- simlist[[i]][[y]]$burned_ids_list |> length() # number of fires that year
    if(nf < 1) next
    for(f in 1:nf) {
      # f = 1
      idsmat <- simlist[[i]][[y]]$burned_ids_list[[f]]
      for(cell in 1:ncol(idsmat)) {
        bp[idsmat[1, cell], idsmat[2, cell]] <- bp[idsmat[1, cell], idsmat[2, cell]] + 1
      }
    }
  }

  # divide by nsim (years) to get annual prob
  bp <- bp / nsim
  r <- rast_from_mat(bp, pnnh_land_rast)
  r <- mask(r, pnnh)
  rprob_list[[i]] <- r
  prob_means[i] <- mean(values(r), na.rm =T)
}

rprob <- do.call("c", rprob_list)
plot(rprob)
