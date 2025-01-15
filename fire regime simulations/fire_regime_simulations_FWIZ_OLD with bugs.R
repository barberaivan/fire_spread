# Notes -------------------------------------------------------------------

# Simulation of fire regime in the Parque Nacional Nahuel Huapi (PNNH). To 
# reduce border effects, I use a 10 km buffer around the park. Ignitions will be
# simulated only inside that region (vector pnnh_buff). The raster including
# that polygon is pnnh_rast. For spread I provide a slightly larger raster: 
# pnnh_land_rast.
 
# The ignition model was fitted using as area the whole area of the park,
# but ignitions were possible only in its burnable portion. The factor to 
# correct this when using a new area with different burnable portion is to
# mutlitply the new burnable area by the total / burnable ration in the park,
# which is ~ 1.273858.

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
library(brms)

source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for rast_from_mat and a few constants

source(file.path("flammability indices",
                 "flammability_indices_functions.R"))

source(file.path("weather",
                 "fortnight_functions.R"))

# Figure size settings ----------------------------------------------------

a4h <- 29.7
a4w <- 21.0
margins <- 2.5
fig_width_max <- a4w - margins * 2
fig_height_max <- a4h - margins * 2

# Functions ---------------------------------------------------------------

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

# Inverse-logit scaled between L and U. If x is a matrix with as many columns as
# elements in L and U, the transform is done column-wise. The same happens if 
# x is a vector as long as L and U. 
# If x is an array and L and U are scalars, the same transform is made to the whole 
# array.
# It's named "2" to avoid stepping onto a simpler function used in the mcmc code.
invlogit_scaled2 <- function(x, L, U) {
  if(is.matrix(x)) {
    out <- x
    for(i in 1:ncol(x)) out[, i] <- plogis(x[, i]) * (U[i] - L[i]) + L[i]
    return(out)
  } else {
    return(plogis(x) * (U - L) + L)
  }
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
# ig_rowcol: matrix[{row, col}, 1], with 1-indexing.
# land: spread landscape.
# steps_min: minimum steps value. Default is 2, so escaped ignitions really 
#   escape. Otherwise, they may be 1, and that burns only the ignition cell.
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
  steps <- floor(coef[n_coef])

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

# Function to extract the FWI data for a given year, given a time-series raster
# with more than one year (fwi_raster_ts). 
#   The FWI raster has the climate of the year to be simulated. 
#   Its layers are the FWI for all pixels in each fortnight, in decreasing order. 
#   As FWI values are used lagged and accumulated, its 10 last layers correspond 
#   to the 10 fortnights before the first one in the focal year.
#   Note that for the ignition model, the spatial average must be computed.
get_fwi_raster <- function(year, fwi_raster_ts) {
  
  fort_max <- date2fort(as.Date(paste(as.character(year), 06, 30, sep = "-")))
  fort_min <- fort_max - nfy - nlags + 1
  f_num_seq <- seq(fort_max, fort_min, by = -1)
  
  # get raster
  fwi_raster <- fwi_raster_ts[[as.character(f_num_seq)]]
  
  return(fwi_raster)
}

# Average spatially the FWI in the PNNH. Returns a vector with the same lenght
# as nlyr(fwi_raster).
spatial_avg_fwi <- function(fwi_raster) {
  fwivals <- terra::extract(fwi_raster, pnnh_buff, exact = T, 
                            weights = T, raw = T)
  outcols <- which(colnames(fwivals) %in% c("ID", "weight"))
  pix_weights <- normalize(fwivals[, "weight"])
  fwi_vals <- as.vector(pix_weights %*% fwivals[, -outcols])
  names(fwi_vals) <- names(fwi_raster)
  return(fwi_vals)
}

# Fucntion to update the spread landscape (land_wrapper$pnnh_land_dyn) after a fire occurs.
# row_cols: matrix[{row, col}, ncells] with the rows and columns of the cells
#   that must be replaced by the value replace.
# replace: 99, which is non-burnable in the spread model.
vegetation_update <- function(row_cols, replace = 99, d = 1) {
  ncell <- ncol(row_cols)
  for(cell in 1:ncell) {
    land_wrapper$pnnh_land_dyn[row_cols[1, cell], row_cols[2, cell], d] <<- replace
  }
}

# Function to simulate one fire season (year, centred at summer).
# fwi_raster: The FWI raster has the climate of the year to be simulated. 
#   Its layers are the FWI for all pixels in each fortnight, in decreasing order. 
#   As FWI values are used lagged and accumulated, its 10 last layers correspond 
#   to the 10 fortnights before the first one in the focal year.
#   Note that for the ignition model, the spatial average must be computed.
#      
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
  
  ### TEST
  # fwi_raster = get_fwi_raster(2098, fwi_source_list[[1]])
  # steps_int_shift = -0.95
  # ig_rate_factor = 1
  # f = 16
  ###
  
  ## Placeholders to save results
  ig_count_h <- 0
  ig_count_l <- 0
  size_fwi_table <- NULL
  burned_ids_list <- NULL
  
  # extract fwi data
  fwi_vec <- spatial_avg_fwi(fwi_raster) # average FWI for pnnh
  
  # Duplicate landscape, to record burned pixels as non-burnable in the veg layer
  land_wrapper$pnnh_land_dyn <- pnnh_land

  # Duplicate cell ids from pnnh, to remove the ones already burned
  cells_pnnh_dyn <- cells_pnnh

  ## Loop over fortnights
  for(f in 1:nfy) {
    # Subset FWI for focal fortnight
    fbegin <- nfy - f + 1 # fortnights begin in the most recent to present
    fend <- fbegin + 10
    
    fwi_local <- fwi_raster[[fbegin:fend]]
    fwi_area <- fwi_vec[fbegin:fend]   # spatial average in the whole area
    
    ## Ignition rate

    # sample parameters from ignition model
    iss <- sample(1:ipost, 1)
    # weight FWI temporally
    igw <- exp(-(tseq/imod[iss, "ls"])^2) |> normalize()
    
    # FWI temporal cumulative weighted mean in the whole area
    fwi_ig <- t(fwi_area) %*% igw

    # compute lambdas by pixel, and get weigthed sum
    lambdas_h <- 
      plogis(imod[iss, "a[1]"] + imod[iss, "b[1]"] * fwi_ig) *
      imod[iss, "U[1]"] * study_area_size
    lambdas_l <- 
      plogis(imod[iss, "a[2]"] + imod[iss, "b[2]"] * fwi_ig) *
      imod[iss, "U[2]"] * study_area_size
    lambda_h <- lambdas_h * ig_rate_factor
    lambda_l <- lambdas_l * ig_rate_factor
  
    # compute phi, which increases with lambda
    phis_h <- imod[iss, "phi_a[1]"] + imod[iss, "phi_b[1]"] * lambdas_h
    phis_l <- imod[iss, "phi_a[2]"] + imod[iss, "phi_b[2]"] * lambdas_l
    
    n_ig_h <- rnegbin(1, lambda_h, phis_h)
    n_ig_l <- rnegbin(1, lambda_l, phis_l)

    # If there are no ignitions, jump to next fortnight
    if(n_ig_h + n_ig_l == 0) next

    candidate_cells <- sample(cells_pnnh_dyn, size = nss * 2,
                              replace = F)
    cand_vals <- pnnh_vals[candidate_cells, c("vfi", "tfi", "drz", "dhz")]

    # Sample location for human ignitions
    if(n_ig_h > 0) {
      ig_count_h <- ig_count_h + n_ig_h

      iprob_h <- exp(
        cand_vals[ssh, "vfi"] * imod[iss, "c_vfi_h"] +
        cand_vals[ssh, "tfi"] * imod[iss, "c_tfi_h"] +
        cand_vals[ssh, "drz"] * imod[iss, "c_dist_r"] +
        cand_vals[ssh, "dhz"] * imod[iss, "c_dist_h"]
      )
      # replace NA with 0
      iprob_h[is.na(iprob_h)] <- 0

      cells_ig_h <- sample(candidate_cells[ssh], size = n_ig_h, replace = F,
                           prob = iprob_h)
    } else {
      cells_ig_h <- NULL
    }

    # Sample location for lightning ignitions
    if(n_ig_l > 0) {
      ig_count_l <- ig_count_l + n_ig_l

      iprob_l <- exp(
        cand_vals[ssl, "vfi"] * imod[iss, "c_vfi_l"] +
        cand_vals[ssl, "tfi"] * imod[iss, "c_tfi_l"]
      )

      # replace NA with 0
      iprob_l[is.na(iprob_l)] <- 0

      cells_ig_l <- sample(candidate_cells[ssl], size = n_ig_l, replace = F,
                           prob = iprob_l)
    } else {
      cells_ig_l <- NULL
    }

    ## Escape from ignited cells
    cells_ig <- c(cells_ig_h, cells_ig_l)
    causes <- rep(c("human", "lightning"), c(n_ig_h, n_ig_l))

    vals_ig <- pnnh_vals[cells_ig, c("vfi", "tfi", "drz", "dhz"), drop = F]
    # get FWI at ignited cells
    crds_ig <- xyFromCell(pnnh_rast, cells_ig)
    points <- vect(crds_ig, "points")
    fwi_esc_ts <- terra::extract(fwi_local, points, method = "bilinear",
                                 raw = T)[, -1, drop = F]

    # cumulative FWI for escape model
    ess <- sample(1:epost, 1)
    escw <- exp(-(tseq/emod[ess, "ls"])^2) |> normalize()
    fwi_esc <- fwi_esc_ts %*% escw

    # Design matrix for escape probability
    Xesc <- cbind(vals_ig, fwi_esc)

    # Compute escape probability (logistic regression)
    esc_betas <- emod[ess, c("b_vfi", "b_tfi", "b_drz", "b_dhz", "b_fwi")]
    escprobs <- plogis(emod[ess, "a"] + Xesc %*% esc_betas)

    # Simulate escape
    escape_ids <- runif(length(cells_ig)) < escprobs
    nsimf <- sum(escape_ids)
    if(nsimf == 0) next
    cells_sim <- cells_ig[escape_ids]
    causes <- causes[escape_ids]
    
    # Simulate spread from escaped cells

    # Translate cell id from PNNH raster to spread raster, which is larger
    crds_sim <- xyFromCell(pnnh_rast, cells_sim)
    cells_sim_land <- cellFromXY(pnnh_land_rast, crds_sim)
    ig_rowcols <- rowColFromCell(pnnh_land_rast, cells_sim_land) |> t()# 1-indexing
    
    # evaluate whether the escaped cell is within the pnnh (it might be outside)
    points_esc <- vect(crds_sim, "points")
    inside_park <- relate(points_esc, pnnh, "intersects") |> as.vector()
    
    # compute accumulated FWI for spread
    fwi_spread_ts <- fwi_esc_ts[escape_ids, , drop = F]
    # standardize
    fwi_spread_ts <- (fwi_spread_ts - fwi_mean_spread) / fwi_sd_spread
    spreadw <- exp(-(tseq/ls_fwi_spread)^2) |> normalize()
    fwi_spread <- fwi_spread_ts %*% spreadw

    # Simulate spread parameters using fwi_spread
    Xspread <- cbind(rep(1, nrow(fwi_spread)), fwi_spread)
    sss <- sample(1:spost, 1)
    ab <- smod$fixef[1:n_coef, 1:2, sss]
    s <- smod$fixef[1:n_coef, 3, sss] |> sqrt()
    rho <- smod$rho[, , sss]
    V <- diag(s) %*% rho %*% diag(s)
    mu <- Xspread %*% t(ab)
    mu[, n_coef] <- mu[, n_coef] + steps_int_shift
    coefs_raw <- mgcv::rmvn(nrow(Xspread), mu, V)
    params_upper[n_coef] <- smod$stepsU[sss]
    coefs <- invlogit_scaled2(coefs_raw, params_lower, params_upper)
    
    # simulate each fire separately, but randomizing order
    idsfire <- sample(1:nsimf, nsimf)
    for(j in idsfire) {
      
      #### 
      #### acá faltaría un if, que haga que si el punto de ignición
      #### ya se quemó por otro fuego en la misma quincena, no se simule.
      #### 
      
      # get burned ids
      burned <- simulate_one_fire(coefs[j, ],
                                  ig_rowcols[, j, drop = F],
                                  land_wrapper$pnnh_land_dyn)
    
      # complete size_fwi_table
      fire_size <- ncol(burned)
      mmat <- data.frame(size = fire_size,
                         steps = floor(coefs[j, n_coef]), 
                         fwi_focal_spread = fwi_spread_ts[j, 1],
                         fwi_cum_spread = fwi_spread[j],
                         cause = causes[j],
                         inside_park = inside_park[j],
                         crds_x = crds_sim[j, "x"],
                         crds_y = crds_sim[j, "y"])
      size_fwi_table <- rbind(size_fwi_table, mmat)

      # expand burned_ids_list
      burned_ids_list <- c(burned_ids_list, list(burned))
      
      # set burned cells as non-burnable in the vegetation layer of 
      # land_wrapper$pnnh_land_dyn
      # (in-place modification is the most efficient)
      vegetation_update(burned)
            
      # remove burned cells from candidate_cells in pnnh. First, translate cells
      # from large landscape to small (pnnh) one.
      cells_large <- cellFromRowCol(pnnh_land_rast,
                                    row = burned[1, ],
                                    col = burned[2, ])
      crds_ <- xyFromCell(pnnh_land_rast, cells_large)
      cells_small <- cellFromXY(pnnh_rast, crds_)
      # include the ones burned by ignitions!
      cells_burned <- c(cells_small, cells_ig) |> unique()
      
      # Use bitmask for efficiency
      max_value <- max(c(cells_pnnh_dyn, cells_burned), na.rm = T)
      bitmask <- rep(FALSE, max_value)
      bitmask[cells_burned] <- TRUE
      keep <- !bitmask[cells_pnnh_dyn]
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
  rm(cells_pnnh_dyn)
  gc()
  return(out)
}

# the same function, but simulating the year inside
simulate_fire_year <- function(year, fwi_source,
                               steps_int_shift = 0,
                               ig_rate_factor = 1) {
  
  # # TEST
  # year = 2098
  # fwi_source = fwi_source_list[[1]]
  # steps_int_shift = -0.95
  # ig_rate_factor = 1
  # # 
  simulate_fire_season(get_fwi_raster(year, fwi_source), 
                       steps_int_shift, ig_rate_factor)
}

# Simulate in parallel
simulate_fire_year_parallel <- function(nsim = 100, cores = 2,
                                        fwi_source, scenario, decade,
                                        steps_int_shift = 0,
                                        ig_rate_factor = 1) {
  
  # ## Test
  # nsim = 20
  # cores = 2
  # fwi_source = fwi_modern
  # steps_int_shift = 0
  # ig_rate_factor = 1
  # scenario = "ssp126"
  # decade = "2040"
  # ## 
  
  registerDoMC(cores)
  
  # routine for modern period
  if(!missing(fwi_source)) {
    message("Simulating modern period.")
    years <- sample(modern_years, nsim, T) |> as.list()
    # simulate in parallel
    result <- foreach(yy = years) %dopar% {
      simulate_fire_year(yy, fwi_source, steps_int_shift, ig_rate_factor)
    }
    return(result)
  } else {
    
    message(paste("Decade ", decade, "; scenario ", scenario, sep = ""))
    
  # routine for future years (uses scenario and decade) 
    years <- sample(dec_years[[decade]], nsim, T)
    simulations <- sample(1:nclimsim, nsim, T, prob = modmem$weight)
    arglist <- lapply(1:nsim, function(i) {
      list("year" = years[i], "simulation" = simulations[i])
    })
    
    message("Loading rasters.")
    # load rasters for target decade and scenario.
    fwi_source_list <- vector("list", nclimsim)
    for(i in 1:nclimsim) {
      key_words <- c(modmem$model[i], scenario, modmem$emember[i], decade)
      patt <- paste(key_words, collapse = ".*")
      row <- grep(patt, proj_files)
      fwi_source_list[[i]] <- rast(file.path(clim_dir, proj_files[row]))
    }
    
    message("Simulating fires.")
    # simulate in parallel
    result <- foreach(arg = arglist) %dopar% {
      simulate_fire_year(year = arg$year, 
                         fwi_source = fwi_source_list[[arg$simulation]], 
                         steps_int_shift, ig_rate_factor)
    }
    rm(fwi_source_list)
    gc()
    return(result)
  }
}

# Simulate sequentially. 
simulate_fire_year_sequential <- function(nsim = 100, sims_before = 0,
                                          fwi_source, scenario, decade,
                                          steps_int_shift = 0,
                                          ig_rate_factor = 1) {
  
  # ## Test
  # nsim = 20
  # fwi_source = fwi_modern
  # steps_int_shift = 0
  # ig_rate_factor = 1
  # scenario = "ssp126"
  # decade = "2040"
  # ## 

  # routine for modern period
  if(!is.list(fwi_source)) {
    message("Simulating modern period.")
    years <- sample(modern_years, nsim, T) |> as.list()
    # simulate in parallel
    result <- foreach(yy = years) %dopar% {
      simulate_fire_year(yy, fwi_source, steps_int_shift, ig_rate_factor)
    }
    return(result)
  } else {
    
    # routine for future years (uses scenario and decade) 
    years <- sample(dec_years[[decade]], nsim, T)
    simulations <- sample(1:nclimsim, nsim, T, prob = modmem$weight)
    
    # simulate years sequentially
    result <- vector("list", nsim)
      
    for(i in 1:nsim) {
      message(i + sims_before)
      
      result[[i]] <- simulate_fire_year(
        year = years[i], 
        fwi_source = fwi_source[[simulations[i]]], 
        steps_int_shift, ig_rate_factor
      )
    }
    return(result)
  }
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

# Make qqplot comparing the size distribution of simulated and observed fires,
# at log10 scale. Its argument is a list of size tables (st), one for each
# simulation experiment. titles is a vector of the same length indicating the 
# titles to give to the plots.
qqcompare <- function(stlist, titles, area_lower = 0.09,
                      xobs = igdata$area_impute2, rc = c(4, 3)) {
  xobs <- log10(xobs[xobs > area_lower])
  bdd <- numeric(length(stlist))
  ll <- length(stlist)
  par(mfrow = c(rc[1], rc[2]))
  for(i in 1:ll) {
    xsim <- stlist[[i]]$size * 0.09
    xsim <- log10(xsim[xsim > area_lower])
    
    # get batt distance
    bd <- round(bhattacharyya_distance(xsim, xobs), 3)
    bdd[i] <- bd
    tit <- round(titles[i], 3)
    tt <- paste("par = ", tit, "; bd = ", bd, sep = "")
    qqplot(xobs, xsim, main = tt, pch = 19); abline(0, 1)
  }
  par(mfrow = c(1, 1))
  return(bdd)
}

# Function to compute the annual burn probability in the PNNH, based on many
# simulated years.
burnprob <- function(simlist, titles, refval = 0.0018) {
  
  ll <- length(simlist)
  tab <- data.frame(
    par = titles,
    bp = numeric(ll)
  )
  
  # list of rasters 
  rlist <- vector("list", ll)
  
  # matrix to fill
  bp0 <- matrix(0, nrow(pnnh_land), ncol(pnnh_land))
  # set non-burnable as -1
  bp0[pnnh_land[, , "veg"] == 99] <- NA
  
  for(i in 1:ll) {
    print(i)
    bpmat <- bp0
    
    # add 1 in burned pixels
    nsim <- length(simlist[[i]])
    for(y in 1:nsim) {
      # print(y)
      # i = 1; y = 1
      nf <- simlist[[i]][[y]]$burned_ids_list |> length() # number of fires that year
      if(nf < 1) next
      for(f in 1:nf) {
        # f = 1
        idsmat <- simlist[[i]][[y]]$burned_ids_list[[f]]
        for(cell in 1:ncol(idsmat)) {
          bpmat[idsmat[1, cell], idsmat[2, cell]] <- bpmat[idsmat[1, cell], idsmat[2, cell]] + 1
        }
      }
    }
    
    # divide by nsim (years) to get annual prob
    bpmat <- bpmat / nsim
    r0 <- rast_from_mat(bpmat, pnnh_land_rast)
    r <- mask(r, pnnh)
    tab$bp[i] <- mean(values(r), na.rm = T)
    rlist[[i]] <- r
  }
  
  out <- list(
    bprast = do.call("c", rlist),
    bpdf = tab
  )
  
  plot(bp ~ par, data = tab, pch = 19)
  abline(h = refval, col = 2)
  
  return(out)
}

map_burnprob <- function(sims) {
  # matrix to fill
  bpmat <- matrix(0, nrow(pnnh_land), ncol(pnnh_land))
  # set non-burnable as -1
  bpmat[pnnh_land[, , "veg"] == 99] <- NA
  
  # add 1 in burned pixels
  nsim <- length(sims)
  for(y in 1:nsim) {
    # print(y)
    # i = 1; y = 1
    nf <- sims[[y]]$burned_ids_list |> length() # number of fires that year
    if(nf < 1) next
    for(f in 1:nf) {
      # f = 1
      idsmat <- sims[[y]]$burned_ids_list[[f]]
      for(cell in 1:ncol(idsmat)) {
        bpmat[idsmat[1, cell], idsmat[2, cell]] <- bpmat[idsmat[1, cell], idsmat[2, cell]] + 1
      }
    }
  }
  
  # divide by nsim (years) to get annual prob
  bpmat <- bpmat / nsim
  r0 <- rast_from_mat(bpmat, pnnh_land_rast)
  rm(bpmat); gc()
  return(r0)
}

average_burnprob <- function(r) {
  r <- mask(r, pnnh)
  return(mean(values(r), na.rm = T))
}





# Simulation settings -----------------------------------------------------

modern_years <- 1999:2022
tseq <- 0:10 # temporal lagg for FWI (from focal to -10 fortnight)
nt <- length(tseq)
export_dir <- file.path("files", "fire_regime_simulation_FWIZ")
# number of candidate cells to spatially sample ignitions
nss <- 4000
ssh <- 1:nss
ssl <- (nss+1):(nss*2)

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

params_lower <- c(-ext_alpha, rep(0, n_coef-2), 2)
params_upper <- c(ext_alpha, rep(ext_beta, n_coef-2), NA)
names(params_lower) <- names(params_upper) <- par_names
params_upper["slope"] <- ext_beta / slope_sd

support <- rbind(params_lower, params_upper)
colnames(support) <- names(params_lower) <- names(params_upper) <- par_names
support_width <- apply(support, 2, diff)

fwi_spread_mean_sd <- readRDS(
  file.path("files", "hierarchical_model_FWIZ", "fwi_mean_sd_spread.rds")
)

fwi_mean_spread <- fwi_spread_mean_sd$fwi_mean
fwi_sd_spread <- fwi_spread_mean_sd$fwi_sd

# FWI lengthscale (cumulative effect)
ls_fwi_spread <- 2.820506
# from <weather/FWI fortnight matrix for spread and lengthscale estimation.R>, line ~ 313

nlags <- 10

# Spatial data ------------------------------------------------------------

# Nahuel Huapi National Park (strict)
pnnh <- vect(file.path("data", "protected_areas", "apn_limites.shp"))
pnnh <- pnnh[pnnh$nombre == "Nahuel Huapi", ]
pnnh <- project(pnnh, "EPSG:5343")

# Buffer around PNNH (10 km buffer, to simulate ignitions; posgar 2007)
pnnh_buff <- vect(file.path("data", "protected_areas", "pnnh_buff_10000.shp"))

# Nahuel Huapi raster (with distance to humans, to simulate ignitions)
pnnh_rast <- rast(file.path("data", "pnnh_images",
                            "pnnh_data_30m_buff_10000.tif"))
## Load file with standardized variables (code in prepare PNNH raster)

# Nahuel Huapi raster (without distance to humans!)
pnnh_land_rast <- rast(file.path("data", "pnnh_images",
                                 "pnnh_data_spread_buffered_30m.tif"))

# Spread landscape
pnnh_land <- readRDS(file.path("data", "pnnh_images",
                               #"pnnh_spread_landscape.rds"))
                               "pnnh_spread_landscape_urban-nonburnable.rds"))

# turn pnnh_land into raster to evaluate which are the burnable cells in the 
# landscape
pnnh_rast_veg <- rast_from_mat(pnnh_land[, , "veg"], pnnh_land_rast$veg)

land_rows <- dim(pnnh_land)[1]
land_cols <- dim(pnnh_land)[2]

# wrapper so the pnnh_land_dyn is inside the global environment
land_wrapper <- list(
  pnnh_land_dyn = pnnh_land
)

# ncell(pnnh_land_rast) == prod(dim(pnnh_land)[1:2]) # OK

# metrics to standardize predictors
pnnh_data_summary <- readRDS(file.path("data", "pnnh_images",
                                       "pnnh_data_summary.rds"))

dr_mean <- pnnh_data_summary$dr_mean  # distance from roads
dr_sd <- pnnh_data_summary$dr_sd
dh_mean <- pnnh_data_summary$dh_mean  # distance from human settlements
dh_sd <- pnnh_data_summary$dh_sd

## Observed fires data, in PNNH, to compare size distribution with simulations
igdata <- read.csv(file.path("data", "ignition", "ignition_size_data.csv"))
igdata$area_impute2 <- igdata$area_impute
igdata$area_impute2[igdata$area_impute < 0.09] <- 0.09

## Mappaed fires, also to compare fire size distribution.
mapfires <- read.csv(file.path(
  "data",
  "climatic_data_by_fire_fwi-fortnight-cumulative_FWIZ.csv"
))

# Annual burn prob in PNNH ------------------------------------------------

# To compare with the simulated one. In the four national parks (Lanin-Alerces),
# it was 0.0018 in 1999-2022.

# firesmap <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp")
# firesmap <- project(firesmap, "EPSG:5343")
# firesmap <- terra::intersect(firesmap, pnnh)
# ba_pnnh <- sum(expanse(firesmap, unit = "ha"))
# pnnh_rast_mask <- mask(pnnh_rast, pnnh)
# pix_burnable <- sum(values(pnnh_rast_mask$veg) < 9, na.rm = T)
# area_burnable <- sum(pix_burnable * 0.09)
# bp_pnnh <- (ba_pnnh / area_burnable) / length(1999:2022)
# 0.0009896007
bp_ref <-  c("pnnh" = 0.0009896007, "parks" = 0.0018, "all" = 0.0025)

# Correct area for ignition rate in new landscape -------------------------

# pnnh_rast_buff_mask <- mask(pnnh_rast, pnnh_buff)
# burnable_buff <- sum(values(pnnh_rast_buff_mask$veg) < 9, na.rm = T) * 0.09 / 100 # km2
# total_pnnh <- expanse(pnnh, unit = "km")
# pnnh_rast_mask <- mask(pnnh_rast, pnnh)
# burnable_pnnh <- sum(values(pnnh_rast_mask$veg) < 9, na.rm = T) * 0.09 / 100
# total_burnable_r <- total_pnnh / burnable_pnnh

total_burnable_r <- 1.273858 # ratio total / burnable in pnnh
burnable_buff <- 9982.515    # burnable in buffer around pnnh, km2

study_area_size <- burnable_buff * total_burnable_r

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

# Import fire models -------------------------------------------------------

igmod <- readRDS(file.path("files", "ignition_FWIZ", "ignition_model_samples.rds"))
escmod <- readRDS(file.path("files", "ignition_FWIZ", "escape_model_samples.rds"))
smod <- readRDS(file.path("files", "hierarchical_model_FWIZ",
                          "spread_model_samples.rds"))

# extract parameters from stan models
imod <- as.matrix(
  igmod,
  pars = c("a", "b", "U", "phi_a", "phi_b", "ls",
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
#             file.path("data", "pnnh_images", "pnnh_data_30m_buff_10000_std.tif"))
pnnh_rast_full <- rast(file.path("data", "pnnh_images", "pnnh_data_30m_buff_10000_std.tif"))

# Define ignitable cells (inside park's buffer and burnable)
pnnh_rast <- mask(pnnh_rast_full, pnnh_buff)
cells_pnnh0 <- cells(pnnh_rast) # all cells inside park
pnnh_vals <- values(pnnh_rast)  # all values, inside or not (does not make a diff)

keep <- pnnh_vals[, "veg"] != 8 &         # remove urban
        pnnh_vals[, "vegetation"] < 6 &   # remove unburnable
        !is.na(pnnh_vals[, "vegetation"]) # remove NA (outside)

cells_pnnh <- cells_pnnh0[keep] 
ncells_pnnh <- length(cells_pnnh)

# Climatic data -----------------------------------------------------------

# FWI
fwi_modern <- rast(file.path("data", "fwi_daily_1998-2022", "24km",
                             "fwi_fortnights_19970701_20230630_standardized_pnnh.tif"))
# already in posgar 2007 and cropped to the park, to save memory.


# weights and emember names for the climatic models
clim_dir <- file.path("data", "fwi_projections", "fwi_fortnights_standardized")
clim_dir0 <- file.path("data", "fwi_projections")

mtable <- read.csv(file.path(clim_dir0, "models_table.csv"))
modmem <- read.csv(file.path(clim_dir0, "models_emembers_table.csv"))
models_n <- aggregate(emember ~ model, modmem, length)

# divide each model's weight between its ensemble members (flat)
modmem$weight <- rep(mtable$weight / models_n$emember, models_n$emember)
nclimsim <- nrow(modmem)

# files with future climate
proj_files <- list.files(clim_dir, pattern = ".tif$")

scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")
ns <- length(scenarios)
decades <- c("2040", "2090")
nd <- length(decades)

dec_years <- list(
  "2040" = 2040:2049,
  "2090" = 2090:2099
)


# Small trial with current climate ----------------------------------------

out <- simulate_fire_year_parallel(2000, 3,
                                   fwi_source = fwi_modern,
                                   steps_int_shift = -0.95,
                                   ig_rate_factor = 1)

mm <- map_burnprob(out)
bp <- average_burnprob(mm)

plot(mm)
plot(pnnh_buff, add = T)
plot(pnnh, add = T)

# guardado como 
# simulation modern trial -1
# bp = 0.001028144
# está muy bien eso, pero se quema muchísimo la ciudad. 
# Y si ponemos no quemable lo urbano cerca de las ciudades? Si, hacer eso.


## ELEVATION < MIN(PARQUE) SETEARLA EN MIN(PARQUE), PARA QUE EL
## MANSO NO SE QUEME TANTO.

## URBANO EN BARILOCHE Y VLA SETEARLO COMO NO QUEMABLE.

# Ahora corrí con urban = non-burnable
# bp = 0.0008827982

# Pruebo lo mismo con shift = -0.95
# bp = 0.0009557399

# Pruebo lo mismo con shift = -0.9
# bp = 0.001202773

# Small search over intercept-shift space ---------------------------------
# [RAN OVER THE STRICT PNNH POLYGON, SO IT'S OUTDATED]

# Decrease intercept, leave ignition rate as fitted, and set the minimum steps
# at 2. That's because the simulator takes steps = 1 as ignition.
shiftvals <- seq(0, -2, by = -0.2)
ns <- length(shiftvals)

# not running in batches, all at once.
Sys.time()
for(i in 7:ns) {
  # i = 1
  print(i)
  out <- simulate_fire_year_parallel(3000, 3,
                                     steps_int_shift = shiftvals[i],
                                     ig_rate_factor = 1)
  num <- str_pad(i, 2, pad = "0")
  nn <- paste("sim_search_", num, ".rds", sep = "")
  saveRDS(out, file.path(export_dir, nn))
  rm(out); gc()
}
Sys.time()
# First 6 runs (0 to -1):
# start: "2024-11-27 01:10:04 -03"
# end:   "2024-11-27 03:52:23 -03"

# Remaining 5 runs (-1.2 to -2):
# start: "2024-11-27 13:35:21 -03"
# end:   "2024-11-27 15:46:04 -03"

# load files
simss <- lapply(list.files(export_dir, pattern = "sim_search_"), function(tt) {
  readRDS(file.path(export_dir, tt))
})

# list of size tables
stlist <- lapply(simss, function(x) {
  do.call("rbind", lapply(x, function(sim) sim[["size_fwi_table"]]))
})

## qqplots
# against fire size distribution in PNNH
bbb <- qqcompare(stlist, shiftvals, rc = 4:3)
# Fire size distribution of mapped fires > 10 ha
bbb2 <- qqcompare(stlist, shiftvals, area_lower = 10, xobs = mapfires$area_ha,
                  rc = 4:3)

dfcomp <- data.frame(shift = shiftvals, distance = bbb, distance2 = bbb2)

plot(distance ~ shift, dfcomp, pch = 19, 
     ylab = "Bhattacharyya distance (log10 size)",
     xlab = "steps intercept shift")

plot(distance2 ~ shift, dfcomp, pch = 19, 
     ylab = "Bhattacharyya distance (log10 size)",
     xlab = "steps intercept shift")

## Compare annual burn probability
# bplist <- burnprob(simss, shiftvals, refval = 0.0018)
# Error: C stack usage  7974356 is too close to the limit

pfit <- data.frame(
  par = shiftvals,
  bp = NA
)
ll <- length(simss)
rlist <- vector("list", ll)

for(i in 1:ll) {
  print(i)
  mm <- map_burnprob(simss[[i]])
  rlist[[i]] <- mm
  pfit$bp[i] <- average_burnprob(mm)
}

bprast <- do.call("c", rlist)

plot(bp ~ par, pfit, pch = 19,
     ylab = "annual burn probability",
     xlab = "steps intercept shift",
     ylim = c(0, max(pfit$bp, bp_ref) * 1.1))
abline(h = bp_ref, col = c("blue", "red", "red"), lty  = 2)

plot(bprast[[shiftvals == -1]])
plot(pnnh, add = T)

# Usar -1?


# Simulate modern period --------------------------------------------------

nbatch <- 4
nsim <- 2500
shift <- -0.95

for(b in 2:nbatch) {
  message(paste("Batch", b))
  out <- simulate_fire_year_parallel(nsim, 3,
                                     fwi_source = fwi_modern,
                                     steps_int_shift = shift,
                                     ig_rate_factor = 1)
  nn <- paste("frs_modern_batch_", b, ".rds", sep = "")
  saveRDS(out, file.path(export_dir, nn))
  rm(out); gc()
  message("\n")
}

# Simulate future periods -------------------------------------------------

for(d in 1:nd) {
  message("------------------------------------------------------")
  for(s in 1:ns) {
    
    d = 2
    s = 4
    for(b in 2:nbatch) {
    # for(b in 1:nbatch) {
      message(paste("Batch", b))
      out <- simulate_fire_year_parallel(nsim, cores = 1, #3,
                                         decade = decades[d],
                                         scenario = scenarios[s],
                                         steps_int_shift = shift,
                                         ig_rate_factor = 1)
      nn <- paste("frs_", decades[d], "_", 
                  scenarios[s], "_batch_", b, ".rds", sep = "")
      saveRDS(out, file.path(export_dir, nn))
      rm(out); gc()
      message("\n")
    }
    message("____________________")
  }
}

# 2090 ssp585 used too much RAM, so we run that with 1 core.


# Simulate 2090 ssp585 ----------------------------------------------------

# out <- readRDS(file.path(export_dir, "frs_2090_ssp585_batch_1.rds"))
# length(out)
# lapply(out, is.null) |> unlist() |> sum()
# # already have 1667 simulations

10000-1600
# Simulate 8400 years, in 84 batches of 100
nb <- 84*2 # 168
nsim <- 100/2

decade = "2090"
scenario = "ssp585"

# load rasters for target decade and scenario.
fwi_source_list <- vector("list", nclimsim)
for(i in 1:nclimsim) {
  key_words <- c(modmem$model[i], scenario, modmem$emember[i], decade)
  patt <- paste(key_words, collapse = ".*")
  row <- grep(patt, proj_files)
  fwi_source_list[[i]] <- rast(file.path(clim_dir, proj_files[row]))
}

sims_before <- 0
for(b in 1:nb) {
  message(paste("Batch", b, "---------------------------------------------"))
  
  tmp <- simulate_fire_year_sequential(
    nsim = nsim, batch = b, sims_before = sims_before,
    fwi_source = fwi_source_list, scenario = scenario, decade = decade,
    steps_int_shift = -0.95, ig_rate_factor = 1
  )
  
  bnum <- str_pad(b+1, 3, pad = "0")
  nn <- paste("frs_", decade, "_", 
              scenario, "_batch_", bnum, ".rds", sep = "")
  saveRDS(tmp, file.path(export_dir, nn))

  rm(out); gc()

  sims_before <- nsim * b
}
# eensd


# Check if there are missing simulations ----------------------------------

filenames <- list.files(export_dir)
files_modern <- filenames[grep("frs_modern_batch_", filenames)]

sapply(files_modern, function(name) {
  x <- readRDS(file.path(export_dir, name))
  nulls <- lapply(x, is.null) |> unlist()
  sum(!nulls)
})  # OK

# future scenarios

counts <- matrix(NA, ns, nd)

for(d in 1:nd) {
  message(paste("Decade", decades[d]))
  for(s in 1:ns) {
    message(paste("Scenario", scenarios[s]))
    if(!(d == 2 & s == 4)) {
      patt <- paste("frs_", decades[d], "_", scenarios[s], "_batch_", sep = "")
      files_focal <- filenames[grep(patt, filenames)]
      counts[s, d] <- sapply(files_focal, function(name) {
        x <- readRDS(file.path(export_dir, name))
        nulls <- lapply(x, is.null) |> unlist()
        return(sum(!nulls))
      }) |> sum()
    }
  }
}

# Simulate more years sequentially (recover missing years) ----------------

# rewrite to avoide saving on disk
counts <- matrix(c(10000, 10000, 10000, 8333, 6667, 6667, 3334, NA),
                 ns, nd)

getmore <- data.frame(decade = c("2040", "2090", "2090", "2090"),
                      scenario = c(scenarios[4], scenarios[1:3]),
                      get = 10000 - c(counts[4, 1], counts[1:3, 2]))

G <- nrow(getmore)

for(g in 1:G) {
  
  nget <- getmore$get[g]
  nb0 <- floor(nget / 50)
  rem <- nget %% 50
  
  nsims <- c(rep(50, nb0), rem)
  nb <- length(nsims)
  
  decade <- getmore$decade[g]
  scenario <- getmore$scenario[g]
  
  # load rasters for target decade and scenario.
  fwi_source_list <- vector("list", nclimsim)
  for(i in 1:nclimsim) {
    key_words <- c(modmem$model[i], scenario, modmem$emember[i], decade)
    patt <- paste(key_words, collapse = ".*")
    row <- grep(patt, proj_files)
    fwi_source_list[[i]] <- rast(file.path(clim_dir, proj_files[row]))
  }
  
  sims_before <- 10000 - nget
  
  for(b in 1:nb) {
    message(paste("Decade ", decade, "; scenario ", scenario,
                  "; batch ", b, " ---------------------", sep = ""))
    
    tmp <- simulate_fire_year_sequential(
      nsim = nsims[b], sims_before = sims_before,
      fwi_source = fwi_source_list, scenario = scenario, decade = decade,
      steps_int_shift = -0.95, ig_rate_factor = 1
    )
    
    bnum <- str_pad(b+4, 3, pad = "0")
    nn <- paste("frs_", decade, "_", 
                scenario, "_batch_", bnum, ".rds", sep = "")
    saveRDS(tmp, file.path(export_dir, nn))
    gc()
    
    sims_before <- sims_before + nsims[b]
  }
}



# Summarise simulations results -------------------------------------------
# Along with annual burned proportion in PNNH


filenames <- list.files(export_dir, pattern = "frs")
ftab <- data.frame(filename = filenames)
ftab$decade <- factor("2040", levels = c("2040", "2090", "modern"))
ftab$scenario <- factor("ssp585", levels = c(scenarios, "modern"))

scdec <- expand.grid(
  scenario = factor(scenarios, levels = c(scenarios, "modern")), 
  decade = factor(decades, levels = c(decades, "modern"))
)
scdec$scdec <- paste(scdec$decade, scdec$scenario, sep = "_")
ncomb <- nrow(scdec)

ftab$decade[grep("modern", filenames)] <- "modern"
ftab$scenario[grep("modern", filenames)] <- "modern"

for(i in 1:ncomb) {
  rows <- grep(scdec$scdec[i], filenames)
  ftab$decade[rows] <- scdec$decade[i]
  ftab$scenario[rows] <- scdec$scenario[i]
}

nsim <- 10000
nexp <- ncomb + 1

exp_names <- c("modern", scdec$scdec)

# Crear la bitmask del parque nacional (celdas dentro del parque)
pnnh_land_rast_template <- rast(ext = ext(pnnh_land_rast), 
                                resolution = res(pnnh_land_rast), 
                                crs = crs(pnnh_land_rast),
                                vals = 1)
rmask <- mask(pnnh_land_rast_template, pnnh)

cells_in_pnnh <- !is.na(values(rmask)) & values(pnnh_rast_veg) != 99
# Vector lógico con TRUE para celdas del parque

# Crear bitmask del parque nacional (TRUE para celdas del parque)
bitmask_pnnh <- rep(FALSE, ncell(pnnh_land_rast))
bitmask_pnnh[which(cells_in_pnnh)] <- TRUE
ncell_pnnh <- sum(bitmask_pnnh) # only burnables
ncols_rast <- ncol(pnnh_land_rast)


# Loop over experiments (scenarios * decades + modern = 9)
for(e in 1:nexp) {
  print(e)
  
  if(e == 1) {
    fn <- ftab$filename[ftab$decade == "modern"]
  } else {
    rows <- ftab$decade == scdec$decade[e-1] & 
      ftab$scenario == scdec$scenario[e-1]
    fn <- ftab$filename[rows]
  }

  # Vector para acumulación de quemas (inicializado en cero)
  burned_count <- rep(0, ncell(pnnh_land_rast))
  
  # Vector para proporción quemada en el parque por iteración
  burned_prop <- numeric(nsim)
  
  # record iteration with cells that burned twice
  duplicate <- numeric(nsim)
  
  iter <- 0
  
  # Loop sobre los archivos de simulación
  for (f in 1:length(fn)) {
    # f = 3
    # print(paste("file", f))
    
    # Leer simulaciones
    sims <- readRDS(file.path(export_dir, fn[f]))
    sims <- sims[sapply(sims, function(x) !is.null(x))]
    ns <- length(sims)
    
    # Loop sobre los años simulados
    for (s in 1:ns) {
      # s = 38
      # message(paste(exp_names[e], ", iteration ", s + iter, sep = ""))
      # print(s)
      
      # Verificar si hubo incendios
      if (length(sims[[s]]$burned_ids_list) == 0) next
      
      # Colapsar todas las coordenadas quemadas en una matriz
      burned_ids_all <- do.call("cbind", sims[[s]]$burned_ids_list) |> t()
      
# LOOPEAR VIENDO CUANDO APARECEN PIXELES DUPLICADOS -----------------------
# CUÁNTOS CASOS SON? MAPEARLO, DIBUJARLO. sON CASOS EN QUE EL VEGETATION_UPDATE 
# NO FUNCIONÓ
      
      # Calcular los índices de celdas quemadas
      rows <- burned_ids_all[, 1]
      cols <- burned_ids_all[, 2]
      bcells <- (rows - 1) * ncols_rast + cols  # Índices de celdas
      
      bcells_u <- unique(bcells)
      
      # Actualizar el conteo de quemas en el vector
      burned_count[bcells_u] <- burned_count[bcells_u] + 1
      
      # Calcular la proporción quemada en el parque
      pp <- sum(bitmask_pnnh[bcells_u]) / ncell_pnnh
      
      burned_prop[s + iter] <- sum(bitmask_pnnh[bcells_u]) / ncell_pnnh
      if(length(bcells) > length(bcells_u)) {
        duplicate[s + iter] <- 1
      }
    }
    
    iter <- ns + iter
  }
  
  # Crear el raster de proporción quemada
  burn_prob_rast <- rast(ext = ext(pnnh_land_rast),
                         resolution = res(pnnh_land_rast),
                         crs = crs(pnnh_land_rast),
                         vals = burned_count / nsim)

  nnrast <- paste("burn_prob_map", "-", exp_names[e], ".tif", sep = "")
  nnvec <- paste("burn_prop_distribution", "-", exp_names[e], ".rds", sep = "")
  nndup <- paste("burn_prop_duplicate", "-", exp_names[e], ".rds", sep = "")
  writeRaster(burn_prob_rast, file.path(export_dir, nnrast), overwrite = T)
  saveRDS(burned_prop, file.path(export_dir, nnvec))
  saveRDS(duplicate, file.path(export_dir, nndup))
}

# último escenario, file 3, iter 38, tiene duplicados

# Burn probability maps ---------------------------------------------------

rastnames <- list.files(export_dir, "burn_prob_map")
rastnames <- c(rastnames[9], rastnames[1:8])
ne <- length(rastnames)

bpmaps <- do.call("c", lapply(1:ne, function(i) {
  rast(file.path(export_dir, rastnames[i])) 
}))

plot(bpmaps * 100)

# Burned proportion distribution ------------------------------------------

vecnames <- list.files(export_dir, "burn_prop_distribution")
vecnames <- c(vecnames[9], vecnames[1:8])

dupnames <- list.files(export_dir, "burn_prop_duplicate")
dupnames <- c(dupnames[9], dupnames[1:8])

ne <- length(vecnames)

levs <- c("modern_modern", 
          "2040_ssp126", "2040_ssp245", "2040_ssp370", "2040_ssp585",
          "2090_ssp126", "2090_ssp245", "2090_ssp370", "2090_ssp585")

bprop <- data.frame(
  exp = factor(rep(levs, each = 10000), levels = levs),
  bprop = do.call("c", lapply(1:ne, function(i) {
    readRDS(file.path(export_dir, vecnames[i])) 
  })),
  dup = do.call("c", lapply(1:ne, function(i) {
    readRDS(file.path(export_dir, dupnames[i])) 
  }))
)

bprop$bprop |> range() # ahora sí

aggregate(dup ~ exp, bprop, mean)
# hay duplicados everywhere.

matrix(bprop$dup, ncol = 10)[1, ]

# SIGO CON PROP > 1. METERSE AL LOOP Y VER QUÉ PASA CON SIMULACION --------



sort(bprop$bprop, dec = T)

bprop <- tidyr::separate(bprop, 1, into = c("decade", "scenario"), 
                         sep = "_", remove = F)
head(bprop)

# Plots are too bad
bpmeans <- aggregate(bprop * 100 ~ exp, bprop, mean)
names(bpmeans)[2] <- "bperc"
bpmeans$q <- bpmeans$bperc / bpmeans$bperc[1]

bpmeans <- aggregate(bprop ~ exp, bprop, max)

View(bprop[bprop$bprop > 1, ])

# library(HDInterval) # HDI
# bphdis <- aggregate(bprop * 100 ~ exp, bprop, hdi)

# Fit stan model to estimate posterior distribution of the mean.
# lowval <- 1 / sum(cells_in_pnnh)
# bprop$bprop2 <- bprop$bprop
# bprop$bprop2[bprop$bprop == 0] <- lowval

bprop_model <- brm(bf(bprop ~ exp, phi ~ exp, zi ~ exp), 
                   family = zero_inflated_beta(), data = bprop,
                   chains = 4, cores = 4)
probs_summ <- fitted(bprop_model, newdata = bpmeans, summary = T)




# Tareas ---------------------------------------------------------------



# Dibujar mapas lindos: 
#   burn prob de modelos separados, 
#   burn prob de escenarios simulados,
#   burn prob de fuegos puntuales

# IN PLACE UPDATE

X <- matrix(rnorm(10000 ^ 2), 10000, 10000)

inplace_replace <- function(rows, cols) {
  for(i in 1:500) {
    X[rows[i], cols[i]] <<- NA
  }
}

microbenchmark::microbenchmark(
  "inplace" = {
    rows <- sample(1:10000, 500)
    cols <- sample(1:10000, 500)
    inplace_replace(rows, cols)
  },
  "ninplace" = {
    rows <- sample(1:10000, 500)
    cols <- sample(1:10000, 500)
    for(i in 1:500) {
      X[rows[i], cols[i]] <- NA
    }
  },
  times = 30
) # much faster inplace


# Copia de lo que hice (funciona!)

ipr <- function(r, c) {
  ll$X[r, c] <<- NA
}

ll <- list()
ll$X <- matrix(0, 3, 3)

for(b in 1:3) {
  print(paste("batch", b))
  ll$X <- matrix(1:9, 3, 3)

  for(i in 1:3) {
    for(j in 1:3) {
      ipr(i, j)
      print(ll$X)
    }
  }
}


# Probar de simular un mismo fuego dos veces, -----------------------------
# a ver qué pasa, a ver si funciona el veg update

