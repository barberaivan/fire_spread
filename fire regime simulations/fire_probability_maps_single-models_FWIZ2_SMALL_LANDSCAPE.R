# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(deeptime)
library(rstan)
library(terra)
library(tidyterra)

source(file.path("flammability indices",
                 "flammability_indices_functions.R"))
## Add vegetation recoding here.

# Figure size settings ----------------------------------------------------

a4h <- 29.7
a4w <- 21.0
margins <- 2.5
fig_width_max <- a4w - margins * 2
fig_height_max <- a4h - margins * 2

# Functions ---------------------------------------------------------------

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
    out <- sapply(1:ncol(x), function(i) {
      plogis(x[, i]) * (U[i] - L[i]) + L[i]
    })
    return(out)
  } else {
    return(plogis(x) * (U - L) + L)
  }
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

params_lower <- c(-ext_alpha, rep(0, n_coef-2), 2) 
params_upper <- c(ext_alpha, rep(ext_beta, n_coef-2), 2000)
# dont worry about stepsL and stepsU, they will not be used
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

# Data -------------------------------------------------------------------

# Nahuel Huapi National Park
pnnh <- vect(file.path("data", "protected_areas", "apn_limites.shp"))
pnnh <- pnnh[pnnh$nombre == "Nahuel Huapi", ]
pnnh <- project(pnnh, "EPSG:5343")

# Nahuel Huapi raster
pnnh_rast <- rast(file.path("data", "pnnh_images", "pnnh_data_120m.tif"))

# metrics to standardize predictors
pnnh_data_summary <- readRDS(file.path("data", "pnnh_images",
                                       "pnnh_data_summary.rds"))

dr_mean <- pnnh_data_summary$dr_mean  # distance from roads
dr_sd <- pnnh_data_summary$dr_sd
dh_mean <- pnnh_data_summary$dh_mean  # distance from human settlements
dh_sd <- pnnh_data_summary$dh_sd
fwi_mean <- pnnh_data_summary$fwi_mean
fwi_sd <- pnnh_data_summary$fwi_sd

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

# Models ------------------------------------------------------------------

igmod <- readRDS(file.path("files", "ignition_FWIZ", "ignition_model_samples.rds"))
escmod <- readRDS(file.path("files", "ignition_FWIZ", "escape_model_samples.rds"))
smod <- readRDS(file.path("files", "hierarchical_model_FWIZ",
                          "spread_model_samples.rds"))

# Prepare raster for maps -------------------------------------------------

pnnh_rast$vegetation <- pnnh_rast$veg

pnnh_df <- as.data.frame(values(pnnh_rast))
names(pnnh_df)[names(pnnh_df) == "veg"] <- "veg_focal"
pnnh_df <- left_join(pnnh_df, dveg[, c("veg_focal", "veg_num")],
                     by = "veg_focal")

pnnh_rast$vegetation <- pnnh_df$veg_num
pnnh_rast$vfi <- vfi_calc(values(pnnh_rast$vegetation),
                          values(pnnh_rast$ndvi))
pnnh_rast$tfi <- tfi_calc(values(pnnh_rast$elevation),
                          values(pnnh_rast$aspect),
                          values(pnnh_rast$slope))
pnnh_rast$drz <- (pnnh_rast$dist_roads / 1000 - dr_mean) / dr_sd
pnnh_rast$dhz <- (pnnh_rast$dist_humans / 1000 - dh_mean) / dh_sd
pnnh_rast$fwi_z <- (pnnh_rast$fwi - fwi_mean) / fwi_sd

pnnh_df <- as.data.frame(pnnh_rast)

# Ignition probability -----------------------------------------------------

igpar_h <- as.matrix(igmod, pars = c("c_vfi_h", "c_tfi_h",
                                     "c_dist_r", "c_dist_h")) |> t()
igpar_l <- as.matrix(igmod, pars = c("c_vfi_l", "c_tfi_l")) |> t()
npost <- ncol(igpar_l)

Xig_h <- model.matrix.lm(~ vfi + tfi + drz + dhz - 1, data = pnnh_df,
                         na.action = "na.pass")
Xig_l <- model.matrix.lm(~ vfi + tfi - 1, data = pnnh_df,
                         na.action = "na.pass")

na_h <- apply(Xig_h, 1, anyNA)
na_l <- apply(Xig_l, 1, anyNA)

# if predictions are computed for all samples at once, the RAM is not enough.
# Compute posterior mean with a loop
prob_pred_h <- numeric(nrow(Xig_h) - sum(na_h))
prob_pred_l <- numeric(nrow(Xig_l) - sum(na_l))
weight <- 1 / npost # to avoid sum overflow

Xig_h_complete <- Xig_h[!na_h, ]
Xig_l_complete <- Xig_l[!na_l, ]

for(k in 1:npost) {
  if(k %% 100 == 0) print(k)
  prob_pred_h <- prob_pred_h + exp(Xig_h_complete %*% igpar_h[, k]) * weight
  prob_pred_l <- prob_pred_l + exp(Xig_l_complete %*% igpar_l[, k]) * weight
}
# scale to max = 1
prob_pred_h <- prob_pred_h / max(prob_pred_h)
prob_pred_l <- prob_pred_l / max(prob_pred_l)

prob_pred_h_full <- numeric(nrow(Xig_h))
prob_pred_l_full <- numeric(nrow(Xig_l))
prob_pred_h_full[na_h] <- NA
prob_pred_l_full[na_l] <- NA
prob_pred_h_full[!na_h] <- prob_pred_h
prob_pred_l_full[!na_l] <- prob_pred_l

# add probabilities to raster
pnnh_rast$igprob_h <- prob_pred_h_full
pnnh_rast$igprob_l <- prob_pred_l_full

plot(pnnh_rast$igprob_h, range = c(0, 0.6))
plot(pnnh_rast$igprob_l, range = c(0, 0.5))

# Escape probability (1pix) ------------------------------------------------

# Escape from 1pix (0.09 ha) is 1 - p(size <= 0.09 ha). Such probability is
# computed as plogis(a[1] - linear predictor), because a[1] is the logit
# cumulative probability for the first size class.

esc_par <- as.matrix(escmod, pars = c("a[1]",
                                      "b_vfi", "b_tfi",
                                      "b_drz", "b_dhz")) |> t()
# no fwi effect included

# negate linear predictor
esc_par[2:nrow(esc_par), ] <- esc_par[2:nrow(esc_par), ] * (-1)

Xesc <- model.matrix.lm(~ vfi + tfi + drz + dhz, data = pnnh_df,
                        na.action = "na.pass")
na_esc <- apply(Xesc, 1, anyNA)

# if predictions are computed for all samples at once, the RAM is not enough.
# Compute posterior mean with a loop
prob_esc <- numeric(nrow(Xesc) - sum(na_esc))
npost <- ncol(esc_par)
weight <- 1 / npost # to avoid sum overflow
Xesc_complete <- Xesc[!na_esc, ]

for(k in 1:npost) {
  if(k %% 100 == 0) print(k)
  prob_esc <- prob_esc + (1 - plogis(Xesc_complete %*% esc_par[, k])) * weight
}

prob_esc_full <- numeric(nrow(Xesc))
prob_esc_full[na_esc] <- NA
prob_esc_full[!na_esc] <- prob_esc

# add probabilities to raster
pnnh_rast$escprob <- prob_esc_full

plot(pnnh_rast[[c("igprob_l", "igprob_h", "escprob")]])

# Spread probabiity -------------------------------------------------------

wind_sd <- 1.464333 # from landscapes_preparation.R
wind_mps <- 4.0 # 14.4 km/h = 4 m/s * 3.6
wind_spread <- wind_mps / wind_sd

pnnh_rast$slope_spread <- sin(pnnh_rast$slope * pi / 180)
pnnh_df <- as.data.frame(pnnh_rast)

Xspread <- model.matrix.lm(~ vfi + tfi - 1, data = pnnh_df,
                           na.action = "na.pass")

na_spread <- apply(Xspread, 1, anyNA)

# if predictions are computed for all samples at once, the RAM is not enough.
# Compute posterior mean with a loop
prob_spread <- numeric(nrow(Xspread) - sum(na_spread))
npost <- dim(smod$fixef)[3]
weight <- 1 / npost # to avoid sum overflow
Xspread_complete <- Xspread[!na_spread, ]

for(i in 1:npost) {
  if(i %% 100 == 0) print(i)
  # i = 1
  # mu at unconstrained scale
  mu <- smod$fixef[1:n_coef, "a", i] # fwi asumed = 0, i.e., the mean for spread.

  # Compute choleski factor of vcov matrix for random effects
  sds <- smod$fixef[1:n_coef, "s2", i] |> sqrt()
  V <- diag(sds) %*% smod$rho[, , i] %*% diag(sds)
  Vchol_U <- chol(V)

  # unconstrained random effects
  ranef_unc <- mgcv::rmvn(1, mu, V)
  ranef_unc <- matrix(ranef_unc, nrow = 1)
  colnames(ranef_unc) <- par_names
  ranef <- invlogit_scaled2(ranef_unc, support[1, ], support[2, ])

  # static spread probabiity
  pp <- plogis(
    ranef["intercept"] +
    Xspread_complete %*% ranef[c("vfi", "tfi")]
  )

  prob_spread <- prob_spread + pp * weight
}

prob_spread_full <- numeric(nrow(Xspread))
prob_spread_full[na_spread] <- NA
prob_spread_full[!na_spread] <- prob_spread

# add probabilities to raster
pnnh_rast$spreadprob <- prob_spread_full

plot(pnnh_rast[[c("igprob_l", "igprob_h", "escprob", "spreadprob")]])

# writeRaster(pnnh_rast, file.path("data", "pnnh_images", "pnnh_data_120m_ig-esc-spread-prob_FWIZ.tiff"))