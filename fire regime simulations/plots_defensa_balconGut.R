library(terra)
library(ggplot2)
library(tidyterra)
library(ggspatial)
library(viridis)
library(patchwork)
library(brms)
library(bayestestR) # hdi interval
library(FireSpread)
library(extrafont)
# extrafont::font_import()

library(rstan)

source(file.path("flammability indices",
                 "flammability_indices_functions.R"))
## Add vegetation recoding here.

source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for rast_from_mat and a few constants

theme_set(theme_classic()) # for non maps
theme_set(theme_minimal()) # for maps

# Functions ---------------------------------------------------------------

normalize <- function(x) x / sum(x)

mean_ci <- function(x) {
  qq <- quantile(x, probs = c(0.025, 0.975), method = 8) %>% unname
  return(c("mean" = mean(x), "lower" = qq[1], "upper" = qq[2]))
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

# Custom theme ------------------------------------------------------------

nice_theme <- function() {
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    axis.line = element_line(linewidth = 0.3),
    
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    
    plot.title = element_text(size = 11),
    
    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "white", color = "white")
  )
}

map_theme <- function() {
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    
    axis.text = element_text(size = 12),
    axis.title = element_blank(),
    
    plot.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    
    strip.text = element_text(size = 12),
    strip.background = element_rect(fill = "white", color = "white")
  )
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

# metrics to standardize predictors
pnnh_data_summary <- readRDS(file.path("data", "pnnh_images",
                                       "pnnh_data_summary.rds"))

dr_mean <- pnnh_data_summary$dr_mean  # distance from roads
dr_sd <- pnnh_data_summary$dr_sd
dh_mean <- pnnh_data_summary$dh_mean  # distance from human settlements
dh_sd <- pnnh_data_summary$dh_sd
fwi_mean <- pnnh_data_summary$fwi_mean
fwi_sd <- pnnh_data_summary$fwi_sd

# Load files --------------------------------------------------------------

# vegetation transform data:
dveg <- readxl::read_excel(
  "/home/ivan/Insync/Mapa vegetación WWF - Lara et al. 1999/clases de vegetacion y equivalencias.xlsx",
  sheet = "Sheet2"
)

balcon_raw <- rast(
  "/home/ivan/Insync/Fire spread modelling/fire_spread/data/focal fires data/defensa_tesis_extra/fire_data_raw_with_distance_1999_25j.tif"
)

# Aplica la reclasificación
veg_map0 <- balcon_raw[["veg"]]
veg_map <- as.factor(veg_map0)

# Convierte el raster reclasificado en un factor con etiquetas
levels(veg_map) <- list(data.frame(
  ID = 1:11, 
  Class = c(
    "Bosque húmedo", 
    "Bosque subalpino", 
    "Bosque subalpino",
    "Bosque seco",
    "Matorral",
    "Pastizal",
    "Plantación",
    "Áreas urbanas",
    "Lagos",
    "Altoandino",
    "Altoandino"
  )
))
# plot(veg_map)

# Load fires and ignition points
fires <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp")
fire <- fires[fires$fire_id == "1999_25j", ]

igpoints <- vect("/home/ivan/Insync/Fire spread modelling/fire_spread/data/ignition_points_checked_with_date-fort-fwiz2.shp")
igpoint <- igpoints[igpoints$fire_id_ig == "1999_25j", ]

# Data and predictors map --------------------------------------------------

veg_map_burn <- mask(veg_map, veg_map0, maskvalue = 8:11)
veg_map_unburn <- mask(veg_map, veg_map0, maskvalue = 8:11, inverse = T)
ndvi_lay <- mask(balcon_raw$ndvi_prev, veg_map0, maskvalue = 8:11)
elev_lyr <- mask(balcon_raw$elevation, veg_map0, maskvalue = 8:11)

dist_h_lyr <- mask(balcon_raw$dist_humans, veg_map0, maskvalue = 8:11)
dist_r_lyr <- mask(balcon_raw$dist_roads, veg_map0, maskvalue = 8:11)
values(dist_h_lyr) <- values(dist_h_lyr) / 1000
values(dist_r_lyr) <- values(dist_r_lyr) / 1000

# settings gráficos
bp_vir <- "F"
bp_begin <- 1
bp_end <- 0.1
bp_name <- "Probabilidad\nde quema (%)"
maxcell <- 100000 * 2

park_contour <- "black"
park_lwd <- 0.2

## Vegetation map (burnable areas)

veg_burn <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = veg_map_burn, maxcell = maxcell) +
  scale_fill_viridis(discrete = T, option = "D", 
                     na.value = "transparent", na.translate = F,
                     begin = 0, end = 1,
                     name = "Vegetación") +
  
  # Fire
  geom_spatvector(data = fire, fill = NA, color = park_contour, 
                  linewidth = 0.4) + 
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
# veg_burn

## Vegetation map (non-burnable areas)

veg_unburn <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = veg_map_unburn, maxcell = maxcell * 2) +
  scale_fill_manual(values = viridis(3, option = "E", begin = 0.15, end = 0.85)[c(3, 1, 2)],
                    na.value = "transparent", na.translate = F,
                    name = "Superficie\nno quemable") +
  
  # National park
  geom_spatvector(data = fire, fill = NA, color = park_contour, 
                  linewidth = 0.4) + 
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
# veg_unburn

ndvimap <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = ndvi_lay, maxcell = maxcell) +
  scale_fill_viridis(option = "G", na.value = "transparent",
                     begin = 0, end = 1, direction = -1,
                     name = "NDVI") +
  
  # National park
  geom_spatvector(data = fire, fill = NA, color = "white", 
                  linewidth = 0.4) + 
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
# ndvimap

elev <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = elev_lyr, maxcell = maxcell) +
  scale_fill_viridis(option = "A", na.value = "transparent",
                     begin = 0, end = 1,
                     name = "Altitud\n(m s.n.m.)") +
  
  # National park
  geom_spatvector(data = fire, fill = NA, color = "white", 
                  linewidth = 0.4) + 
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(linewidth = 0.3, color = "gray60")) +
  ggspatial::annotation_scale(
    location = "tr", height = unit(1.5, "mm"),
    bar_cols = c("grey60", "white"),
    text_col = "white",
    text_cex = 1,
    line_width = 0.4
  ) +
  annotation_north_arrow(
    location = "tl", which_north = "true",  
    style = north_arrow_orienteering(
      line_width = 0.2,
      line_col = "black",
      fill = c("grey60", "white"),
      text_col = "white",
      text_family = "serif",
      text_face = NULL,
      text_size = 0,
      text_angle = 0
    ),
    width = unit(6, "mm"),
    height = unit(6, "mm")
  )

# elev 


dr_map <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = dist_r_lyr, maxcell = maxcell) +
  scale_fill_viridis(option = "E", na.value = "transparent",
                     begin = 0, end = 0.9, direction = -1,
                     name = "Distancia a\ncaminos (km)") +
  
  # National park
  geom_spatvector(data = fire, fill = NA, color = "black", 
                  linewidth = 0.4) + 
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
# dr_map

dh_map <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = dist_h_lyr, maxcell = maxcell) +
  scale_fill_viridis(option = "E", na.value = "transparent",
                     begin = 0, end = 0.9, direction = -1,
                     name = "Distancia a\npoblados (km)") +
  
  # National park
  geom_spatvector(data = fire, fill = NA, color = "black", 
                  linewidth = 0.4) + 
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
# dh_map



preds <- (
  veg_burn + veg_unburn + ndvimap + elev + dr_map + dh_map 
) + plot_layout(
  nrow = 2, byrow = T
) &
  theme(legend.key.width = unit(2.5, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.box.just = "left",
        legend.box.margin = margin(0, 0, 0, -1, unit = "mm"),
        legend.justification = "left",
        
        legend.margin = margin(0, 0, 0, 0, unit = "mm"), # Márgenes alrededor de la leyenda
        legend.spacing = unit(10, "mm"), 
        legend.spacing.y = unit(10, "mm"), # Espaciado entre filas de la leyenda
        
        # legend.title = element_text(size = 8),
        # legend.text = element_text(size = 7),
        
        panel.spacing.x = unit(0, "mm"),
        plot.margin = margin(1, 1, 1, 1, unit = "mm"))

ggsave("fire regime simulations/figures_defensa/balcon_predictors.png",
       plot = preds, bg = "white",
       width = 38, height = 20, units = "cm")

# Data and predictors map - no fire --------------------------------------

## Vegetation map (burnable areas)

veg_burn <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = veg_map_burn, maxcell = maxcell) +
  scale_fill_viridis(discrete = T, option = "D", 
                     na.value = "transparent", na.translate = F,
                     begin = 0, end = 1,
                     name = "Vegetación") +
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

## Vegetation map (non-burnable areas)

veg_unburn <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = veg_map_unburn, maxcell = maxcell * 2) +
  scale_fill_manual(values = viridis(3, option = "E", begin = 0.15, end = 0.85)[c(3, 1, 2)],
                    na.value = "transparent", na.translate = F,
                    name = "Superficie\nno quemable") +
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )


ndvimap <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = ndvi_lay, maxcell = maxcell) +
  scale_fill_viridis(option = "G", na.value = "transparent",
                     begin = 0, end = 1, direction = -1,
                     name = "NDVI") +
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

elev <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = elev_lyr, maxcell = maxcell) +
  scale_fill_viridis(option = "A", na.value = "transparent",
                     begin = 0, end = 1,
                     name = "Altitud\n(m s.n.m.)") +
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(linewidth = 0.3, color = "gray60")) +
  ggspatial::annotation_scale(
    location = "tr", height = unit(1.5, "mm"),
    bar_cols = c("grey60", "white"),
    text_col = "white",
    text_cex = 1,
    line_width = 0.4
  ) +
  annotation_north_arrow(
    location = "tl", which_north = "true",  
    style = north_arrow_orienteering(
      line_width = 0.2,
      line_col = "black",
      fill = c("grey60", "white"),
      text_col = "white",
      text_family = "serif",
      text_face = NULL,
      text_size = 0,
      text_angle = 0
    ),
    width = unit(6, "mm"),
    height = unit(6, "mm")
  )


dr_map <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = dist_r_lyr, maxcell = maxcell) +
  scale_fill_viridis(option = "E", na.value = "transparent",
                     begin = 0, end = 0.9, direction = -1,
                     name = "Distancia a\ncaminos (km)") +
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )

dh_map <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = dist_h_lyr, maxcell = maxcell) +
  scale_fill_viridis(option = "E", na.value = "transparent",
                     begin = 0, end = 0.9, direction = -1,
                     name = "Distancia a\npoblados (km)") +
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() + 
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )


preds2 <- (
  veg_burn + veg_unburn + ndvimap + elev + dr_map + dh_map 
) + plot_layout(
  nrow = 2, byrow = T#, 
  # widths = rep(1, 3), heights = rep(1.3, 2)
) &
  theme(legend.key.width = unit(2.5, "mm"),
        legend.key.height = unit(5, "mm"),
        legend.box.just = "left",
        legend.box.margin = margin(0, 0, 0, -1, unit = "mm"),
        legend.justification = "left",
        
        legend.margin = margin(0, 0, 0, 0, unit = "mm"), # Márgenes alrededor de la leyenda
        legend.spacing = unit(10, "mm"), 
        legend.spacing.y = unit(10, "mm"), # Espaciado entre filas de la leyenda
        
        # legend.title = element_text(size = 8),
        # legend.text = element_text(size = 7),
        
        panel.spacing.x = unit(0, "mm"),
        plot.margin = margin(1, 1, 1, 1, unit = "mm"))


ggsave("fire regime simulations/figures_defensa/balcon_predictors_no-fire.png",
       plot = preds2, bg = "white",
       width = 38, height = 20, units = "cm")



# Predictions --------------------------------------------------------------

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




# Burn probability maps ---------------------------------------------------

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

# Prepare raster for maps -------------------------------------------------

balcon_raw <- rast(
  "/home/ivan/Insync/Fire spread modelling/fire_spread/data/focal fires data/defensa_tesis_extra/fire_data_raw_with_distance_1999_25j.tif"
)
balcon_rast <- balcon_raw

names(balcon_rast) <- c(
  "veg", "ndvi", "ndvi_22", "elevation", "slope",
  "aspect", "burned", "dist_humans", "dist_roads"
)

balcon_rast$vegetation <- balcon_raw$veg

balcon_df <- as.data.frame(values(balcon_rast))
names(balcon_df)[names(balcon_df) == "veg"] <- "veg_focal"
balcon_df <- left_join(balcon_df, dveg[, c("veg_focal", "veg_num")],
                     by = "veg_focal")

balcon_rast$vegetation <- balcon_df$veg_num
balcon_rast$vfi <- vfi_calc(values(balcon_rast$vegetation),
                          values(balcon_rast$ndvi))
balcon_rast$tfi <- tfi_calc(values(balcon_rast$elevation),
                          values(balcon_rast$aspect),
                          values(balcon_rast$slope))
balcon_rast$drz <- (balcon_rast$dist_roads / 1000 - dr_mean) / dr_sd
balcon_rast$dhz <- (balcon_rast$dist_humans / 1000 - dh_mean) / dh_sd

# Get FWI from jan-feb
fwi_modern <- rast(file.path("data", "fwi_daily_1998-2022", "24km",
                             "fwi_fortnights_19970701_20230630_standardized_pnnh.tif"))
# already in posgar 2007 and cropped to the park, to save memory.
months <- lubridate::month(time(fwi_modern))
use <- which(months == 2)
fwi_modern_mean <- 2

balcon_rast$fwi_z <- fwi_modern_mean
balcon_df <- as.data.frame(balcon_rast)

# Ignition probability -----------------------------------------------------

igpar_h <- as.matrix(igmod, pars = c("c_vfi_h", "c_tfi_h",
                                     "c_dist_r", "c_dist_h")) |> t()
igpar_l <- as.matrix(igmod, pars = c("c_vfi_l", "c_tfi_l")) |> t()
npost <- ncol(igpar_l)

Xig_h <- model.matrix.lm(~ vfi + tfi + drz + dhz - 1, data = balcon_df,
                         na.action = "na.pass")
Xig_l <- model.matrix.lm(~ vfi + tfi - 1, data = balcon_df,
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
balcon_rast$igprob_h <- prob_pred_h_full
balcon_rast$igprob_l <- prob_pred_l_full

plot(balcon_rast$igprob_h, range = c(0, 0.6))
plot(balcon_rast$igprob_l, range = c(0, 0.5))

# Escape probability (1pix) ------------------------------------------------

# Escape from 1pix (0.09 ha) probability (logistic regression)

esc_par <- as.matrix(escmod, pars = c("a",
                                      "b_vfi", "b_tfi",
                                      "b_drz", "b_dhz",
                                      "b_fwi")) |> t()

Xesc <- model.matrix.lm(~ vfi + tfi + drz + dhz + fwi_z, data = balcon_df,
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
  prob_esc <- prob_esc + plogis(Xesc_complete %*% esc_par[, k]) * weight
}

prob_esc_full <- numeric(nrow(Xesc))
prob_esc_full[na_esc] <- NA
prob_esc_full[!na_esc] <- prob_esc

# add probabilities to raster
balcon_rast$escprob <- prob_esc_full

# writeRaster(balcon_rast, file.path("data", "balcon_images", "balcon_data_120m_ig-esc-prob.tiff"))

plot(balcon_rast[[c("igprob_l", "igprob_h", "escprob")]])

# Spread probabiity -------------------------------------------------------

# balcon_rast <- rast(file.path("data", "balcon_images", "balcon_data_120m_ig-esc-prob.tiff"))

wind_sd <- 1.46
wind_mps <- 4.0 # 14.4 km/h = 4 m/s * 3.6
wind_spread <- wind_mps / wind_sd

balcon_rast$slope_spread <- sin(balcon_rast$slope * pi / 180)
balcon_df <- as.data.frame(balcon_rast)

Xspread <- model.matrix.lm(~ vfi + tfi - 1, data = balcon_df,
                           na.action = "na.pass")

na_spread <- apply(Xspread, 1, anyNA)

# if predictions are computed for all samples at once, the RAM is not enough.
# Compute posterior mean with a loop
prob_spread <- numeric(nrow(Xspread) - sum(na_spread))
npost <- dim(smod$fixef)[3]
weight <- 1 / npost # to avoid sum overflow
Xspread_complete <- Xspread[!na_spread, ]

# standardize at the spraed scale (expquad fortnight)
fwi_spread <- (mean(values(balcon_rast$fwi)) - fwi_mean_spread) / fwi_sd_spread

for(i in 1:npost) {
  if(i %% 100 == 0) print(i)
  # i = 1
  # mu at unconstrained scale
  mu <- smod$fixef[1:n_coef, "a", i] + smod$fixef[1:n_coef, "b", i] * fwi_spread
  
  # Compute choleski factor of vcov matrix for random effects
  sds <- smod$fixef[1:n_coef, "s2", i] |> sqrt()
  V <- diag(sds) %*% smod$rho[, , i] %*% diag(sds)
  Vchol_U <- chol(V)
  
  # unconstrained random effects
  ranef_unc <- mgcv::rmvn(1, mu, V)
  ranef_unc <- matrix(ranef_unc, nrow = 1)
  colnames(ranef_unc) <- par_names
  ranef <- constrain(ranef_unc, support = support) |> t()
  
  # static spread probabiity
  pp <- plogis(
    ranef["intercept", ] +
      Xspread_complete %*% ranef[c("vfi", "tfi"), ]
  )
  
  prob_spread <- prob_spread + pp * weight
}

prob_spread_full <- numeric(nrow(Xspread))
prob_spread_full[na_spread] <- NA
prob_spread_full[!na_spread] <- prob_spread

# add probabilities to raster
balcon_rast$spreadprob <- prob_spread_full

plot(balcon_rast[[c("igprob_l", "igprob_h", "escprob", "spreadprob")]])
writeRaster(balcon_rast, file.path("data", "focal fires data", "defensa_tesis_extra", "balcon_ig-esc-spread_prob.tiff"))


# Plot burn probability by model ------------------------------------------

# Unburnable layer has gaps. Fill them
window <- matrix(1, nrow = 41, ncol = 41)

# Get mode in window
filled <- focal(
  veg_map_unburn,
  w = window,
  fun = function(x, ...) {
    x <- na.omit(x)
    if (length(x) == 0) return(NA)
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  },
  na.policy = "only",    # solo procesa celdas NA
  filename = "",         # sin guardar a disco
  datatype = "INT1U"     # formato entero
)

levels(filled) <- levels(veg_map)

# Replace NA in the original
veg_map_unburn_filled <- cover(veg_map_unburn, filled)
levels(veg_map_unburn_filled) <- levels(veg_map)



balcon_rast <- rast(file.path("data", "focal fires data", "defensa_tesis_extra", "balcon_ig-esc-spread_prob.tiff"))
burn_prob_models <- balcon_rast[[c("igprob_h", "igprob_l", "escprob", "spreadprob")]]
burn_prob_models <- burn_prob_models * 100

# burn prob separate
ll <- vector("list", 4)
titles <- c(
  "Probabilidad relativa de\nignición por humanos",
  "Probabilidad relativa de\nignición por rayos",
  "Probabilidad de escape",
  "Probabilidad de\npropagación"
)

for(m in 1:4) {
  # Burn prob anual
  ll[[m]] <- 
    
    ggplot() +
    
    # Non burnable layer
    geom_spatraster(data = veg_map_unburn_filled, maxcell = maxcell * 2,
                    show.legend = F) +
    scale_fill_manual(values = viridis(3, option = "E", begin = 0.15, end = 0.85)[c(3, 1, 2)],
                      na.value = "transparent", na.translate = F) +
    
    ggnewscale::new_scale_fill() +
    
    # Burn probability map
    geom_spatraster(data = burn_prob_models[[m]], maxcell = maxcell) +
    scale_fill_viridis(option = bp_vir, na.value = "transparent",
                       begin = bp_begin, end = bp_end,
                       name = "%") +

    scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                       labels = function(x) sprintf("%.2f°", x)) +
    scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                       labels = function(x) sprintf("%.2f°", x)) +
    # edit coordinates format
    coord_sf() +
    map_theme() + 
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.key.width = unit(2, 'mm')) + 
    labs(title = titles[m]) 
}


bp1 <- (
  ll[[1]] + 
    ll[[2]] + 
    ll[[3]] + 
    ll[[4]]
) + plot_layout(nrow = 2, byrow = T)

ggsave(
  "fire regime simulations/figures_defensa/balcon_burn_prob_models.png",
  plot = bp1, width = 20, height = 20, units = "cm", bg = "white"
)


# Simulations -------------------------------------------------------------

# settings gráficos
pointsize <- 3.5
bp_vir <- "F"
bp_begin <- 1
bp_end <- 0.1
bp_name <- "Probabilidad\nde quema (%)"
park_contour <- "black"
park_lwd <- 0.2
maxcell <- 100000 * 2

# load balcon landscape for spread
landscape_file <- "/home/ivan/Insync/Fire spread modelling/fire_spread/data/focal fires data/landscapes/1999_25j.rds"
full_data <- readRDS(landscape_file)
spread_data <- full_data[c("landscape", "ig_rowcol",
                           "burned_layer", "burned_ids",
                           "counts_veg")]
land <- spread_data$landscape

fwi_spread <- 2

# ig and escape prob
balcon_df <- as.data.frame(balcon_rast)
probs_df <- balcon_df[, c("igprob_h", "igprob_l", "escprob")]
# assign zero prob to NAs
for (v in 1:3) {
  probs_df[is.na(probs_df[, v]), v] <- 0
}


## Function to simulate fires
firesim <- function(seed, n_human = 2, n_light = 2, ii = 1, save = TRUE) {
  set.seed(seed)
  
  # simulate wich cells burn
  cells_human <- sample(
    1:nrow(probs_df), n_human, prob = probs_df$igprob_h, replace = F
  )
  cells_light <- sample(
    1:nrow(probs_df), n_light, prob = probs_df$igprob_l, replace = F
  )
  
  crds_human <- terra::xyFromCell(balcon_rast, cells_human)
  points_human <- vect(crds_human, "points")
  crs(points_human) <- crs(balcon_rast)
  
  crds_light <- terra::xyFromCell(balcon_rast, cells_light)
  points_light <- vect(crds_light, "points")
  crs(points_light) <- crs(balcon_rast)
  
  # Map them
  ig_human <- 
    ggplot() +
    
    # Non burnable layer
    geom_spatraster(data = veg_map_unburn_filled, maxcell = maxcell * 2) +
    scale_fill_manual(values = viridis(3, option = "E", begin = 0.15, end = 0.85)[c(3, 1, 2)],
                      na.value = "transparent", na.translate = F) +
    
    ggnewscale::new_scale_fill() +
    
    geom_spatraster(data = burn_prob_models[[1]], maxcell = maxcell) +
    scale_fill_viridis(
      option = bp_vir, na.value = "transparent",
      begin = bp_begin, end = bp_end,
      name = "%"
    ) + 
    geom_spatvector(data = points_human, 
                    color = "black", 
                    fill = "cyan1",#viridis(1, begin = 0.1, option = "C"),
                    shape = 21, size = pointsize, stroke = 0.8) +
    map_theme() +
    scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0)) +
    theme(
      legend.position = "none",
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )

  ig_light <- 
    ggplot() +
    
    # Non burnable layer
    geom_spatraster(data = veg_map_unburn_filled, maxcell = maxcell * 2,
                    show.legend = F) +
    scale_fill_manual(values = viridis(3, option = "E", begin = 0.15, end = 0.85)[c(3, 1, 2)],
                      na.value = "transparent", na.translate = F) +
    
    ggnewscale::new_scale_fill() +
    
    geom_spatraster(data = burn_prob_models[[2]], maxcell = maxcell) +
    scale_fill_viridis(
      option = bp_vir, na.value = "transparent",
      begin = bp_begin, end = bp_end,
      name = "%"
    ) + 
    geom_spatvector(data = points_light, 
                    color = "black", 
                    fill = "springgreen2",#viridis(1, begin = 0.45, option = "C"),
                    shape = 21, size = pointsize, stroke = 0.8) +
    scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0)) +
    map_theme() + 
    theme(
      legend.position = "none",
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )

  # Escape probability
  points_try <- rbind(points_human, points_light)
  escape_prob <- probs_df$escprob[c(cells_human, cells_light)]
  escape <- ifelse(runif(nrow(points_try)) < escape_prob, TRUE, FALSE)
  if (sum(escape) > 0) {
    escmap <- 
      ggplot() +
      
      # Non burnable layer
      geom_spatraster(data = veg_map_unburn_filled, maxcell = maxcell * 2,
                      show.legend = F) +
      scale_fill_manual(values = viridis(3, option = "E", begin = 0.15, end = 0.85)[c(3, 1, 2)],
                        na.value = "transparent", na.translate = F) +
      
      ggnewscale::new_scale_fill() +
      
      geom_spatraster(data = burn_prob_models[[3]], maxcell = maxcell) +
      scale_fill_viridis(
        option = bp_vir, na.value = "transparent",
        begin = bp_begin, end = bp_end,
        name = "%"
      ) + 
      geom_spatvector(data = points_try[!escape, ], 
                      color = "black", fill = "gray",
                      shape = 21, size = pointsize, stroke = 0.8) +
      geom_spatvector(data = points_try[escape, ], 
                      color = "black",
                      fill = "orange",#viridis(1, begin = 1, option = "C"),
                      shape = 21, size = pointsize, stroke = 0.8) +
      scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0)) +
      map_theme() +
      theme(
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank()
      )

    # Simulate spread
    cells_escape <- c(cells_human, cells_light)[escape]
    ignition_rc <- rowColFromCell(balcon_rast, cells_escape) |> t()
    
    nf <- ncol(ignition_rc)
    burns <- vector("list", nf)
    
    for (i in 1:nf) {
      # rotate wind to make it variable
      land_rot <- land
      angle <- runif(1, 0, 2 * pi)
      land_rot[, , "wdir"] <- (land[, , "wdir"] + angle) %% (2 * pi)
      
      # Simulate spread parameters using fwi_spread
      Xspread <- cbind(1, fwi_spread)
      sss <- sample(1:spost, 1)
      ab <- smod$fixef[1:n_coef, 1:2, sss]
      s <- smod$fixef[1:n_coef, 3, sss] |> sqrt()
      rho <- smod$rho[, , sss]
      V <- diag(s) %*% rho %*% diag(s)
      mu <- Xspread %*% t(ab)
      mu[, n_coef] <- mu[, n_coef]
      coefs_raw <- mgcv::rmvn(nrow(Xspread), mu, V)
      params_upper[n_coef] <- smod$stepsU[sss]
      coefs <- invlogit_scaled2(coefs_raw, params_lower, params_upper)
      
      # Tidy parameters
      particle <- coefs[1, ]
      coef_intercepts <- rep(particle[1], n_veg)
      coef_nd <- particle[2:3]
      coef_terrain <- particle[4:(n_coef-1)]
      steps <- particle[n_coef]
      
      # Simulate spread
      burns[[i]] <- simulate_fire(
        layer_vegetation = land[, , "veg"],
        layer_nd = land[, , nd_variables],
        layer_terrain = land_rot[, , terrain_variables],
        coef_intercepts = coef_intercepts,
        coef_nd = coef_nd,
        coef_terrain = coef_terrain,
        ignition_cells = ignition_rc[, i, drop = F] - 1, # zero-indexing
        upper_limit = upper_limit,
        steps = steps
      )
    }
    
    burn_array <- abind::abind(burns, along = 3)
    burn_matrix <- apply(burn_array, 1:2, max)
    burn_rast <- rast_from_mat(burn_matrix, balcon_rast)
    
    # vectorize raster
    burn_mask <- burn_rast
    burn_mask[burn_rast == 0] <- NA
    burn_mask_vect <- terra::as.polygons(burn_mask, dissolve = TRUE, na.rm = TRUE)
    
    spread <- 
      ggplot() +
      
      # Non burnable layer
      geom_spatraster(data = veg_map_unburn_filled, maxcell = maxcell * 2,
                      show.legend = F) +
      scale_fill_manual(values = viridis(3, option = "E", begin = 0.15, end = 0.85)[c(3, 1, 2)],
                        na.value = "transparent", na.translate = F) +

      ggnewscale::new_scale_fill() +
      
      # Burn probability map
      geom_spatraster(data = burn_prob_models[[4]], maxcell = maxcell, inherit.aes = F) +
      scale_fill_viridis(option = bp_vir, na.value = "transparent",
                         begin = bp_begin, end = bp_end,
                         name = "%") +
      
      # Fire
      geom_spatvector(data = burn_mask_vect, 
                      fill = "black", 
                      color = NA, 
                      alpha = 0.6) +
      
      # ig points  
      geom_spatvector(data = points_try[escape, ], 
                      color = "black", fill = "orange",
                      shape = 21, size = pointsize, stroke = 0.8) +
      scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0)) +
      map_theme() +
      theme(
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank()
      ) 
    
    all <- (ig_human + ig_light + escmap + spread) + 
      plot_layout(nrow = 2, byrow = T)
    
    print(all)
    
    if (save) {
      fn <- sprintf(
        "fire regime simulations/figures_defensa/balcon_simulations_%02d.png", 
        ii
      )
      
      ggsave(
        fn, plot = all, bg = "white", width = 19, height = 20, units = "cm"
      )
    }
  } else {
    message("No ignition escaped.")
  }
}

# Settings que se ven bien (semilla, nhum, nlight:
firesim(1, 2, 2, ii = 1, save = T)
firesim(2, 2, 2, ii = 2, save = T)
firesim(3, 2, 2, ii = 3, save = T)
firesim(7, 2, 2, ii = 4, save = T)
firesim(11, 2, 2, ii = 5, save = T)
firesim(13, 2, 2, ii = 6, save = T)

# Fire spread animation ---------------------------------------------------

# Choose nice fires -------------------------------------------------------

# Import FWI data for balcon gutierrez
# FWI expquad anomaly by fortnight
fwi_data <- read.csv(file.path("data", "climatic_data_by_fire_fwi-fortnight-cumulative_FWIZ2.csv"))
fwi_fire_anom <- fwi_data$fwi_fort_expquad[fwi_data$fire_id == "1999_25j"]
# Standardize to spread scale
fmsd <- readRDS(file.path(
  "files", "hierarchical_model_FWIZ", "fwi_mean_sd_spread.rds"
))
fwi_fire_spread <- (fwi_fire_anom - fmsd$fwi_mean) / fmsd$fwi_sd

# Get Balcón del Gutiérrez fitted parameters, logit scales, except for steps.
bhat_logit <- smod$ranef[, "1999_25j", ]
nsim <- ncol(bhat_logit)

# constrain random effects
ranef_fit <- invlogit_scaled2(t(bhat_logit), params_lower, params_upper)
ranef_fit[, "steps"] <- bhat_logit["steps", ]
ranef_fit <- t(ranef_fit)

fire_rast <- balcon_rast$burned

burn_fitted <- function(column, seed = NULL, ret = F) {
  
  ss <- ifelse(is.null(seed), column, seed)
  
  # Get parameter
  particle <- ranef_fit[, column]
  
  # Tidy parameters
  coef_intercepts <- rep(particle[1], n_veg)
  coef_nd <- particle[2:3]
  coef_terrain <- particle[4:(n_coef-1)]
  steps <- particle[n_coef]
  
  set.seed(ss)
  
  out <- simulate_fire_compare(
    layer_vegetation = land[, , "veg"],
    layer_nd = land[, , nd_variables],
    layer_terrain = land[, , terrain_variables],
    coef_intercepts = coef_intercepts,
    coef_nd = coef_nd,
    coef_terrain = coef_terrain,
    ignition_cells = spread_data$ig_rowcol, # created with zero-indexing
    upper_limit = upper_limit,
    steps = steps
  )
  
  ov <- overlap_spatial(out, spread_data)
  
  out_rast <- rast_from_mat(out$burned_layer, balcon_rast)
  
  par(mfrow = c(1, 2))
  plot(fire_rast)
  plot(out_rast)
  par(mfrow = c(1, 1))
  print(ov)
  
  if (ret) {
    values(out_rast)[values(out_rast) == 0] <- NA
    out_poly <- terra::as.polygons(out_rast, na.rm = T, dissolve = T)
    out_poly$overlap <- ov
    return(out_poly)
  }
}

# use seed 2 col 2
burn_fitted(2, 2)

# The same, but simulating random effects
burn_simulate <- function(column, seed = NULL, ret = F) {
  
  sss <- ifelse(is.null(seed), column, seed)
  set.seed(sss)
  
  # Simulate spread parameters using fwi_fire_spread
  Xspread <- cbind(1, fwi_fire_spread)
  sss <- sample(1:spost, 1)
  ab <- smod$fixef[1:n_coef, 1:2, sss]
  s <- smod$fixef[1:n_coef, 3, sss] |> sqrt()
  rho <- smod$rho[, , sss]
  V <- diag(s) %*% rho %*% diag(s)
  mu <- Xspread %*% t(ab)
  mu[, n_coef] <- mu[, n_coef]
  coefs_raw <- mgcv::rmvn(nrow(Xspread), mu, V)
  params_upper[n_coef] <- smod$stepsU[sss]
  coefs <- invlogit_scaled2(coefs_raw, params_lower, params_upper)
  
  # Tidy parameters
  particle <- coefs
  coef_intercepts <- rep(particle[1], n_veg)
  coef_nd <- particle[2:3]
  coef_terrain <- particle[4:(n_coef-1)]
  steps <- particle[n_coef]
  
  out <- simulate_fire_compare(
    layer_vegetation = land[, , "veg"],
    layer_nd = land[, , nd_variables],
    layer_terrain = land[, , terrain_variables],
    coef_intercepts = coef_intercepts,
    coef_nd = coef_nd,
    coef_terrain = coef_terrain,
    ignition_cells = spread_data$ig_rowcol, # created with zero-indexing
    upper_limit = upper_limit,
    steps = steps
  )
  
  ov <- overlap_spatial(out, spread_data)
  
  out_rast <- rast_from_mat(out$burned_layer, balcon_rast)
  
  par(mfrow = c(1, 2))
  plot(fire_rast)
  plot(out_rast)
  par(mfrow = c(1, 1))
  print(ov)
  
  if (ret) {
    values(out_rast)[values(out_rast) == 0] <- NA
    out_poly <- terra::as.polygons(out_rast, na.rm = T, dissolve = T)
    out_poly$overlap <- ov
    return(out_poly)
  }
}

# use seed 2 col 2
burn_simulate(4) # nice
burn_simulate(9)
burn_simulate(15) # nice
k <- 17

# very nice
# bigger: 21 23 42
# smaller: 15 32

k <- k+1
burn_simulate(21) # nice
print(k)


# Prepare animation layers ------------------------------------------------

# all burnable
base_layer <- veg_map0
values(base_layer) <- 2
# non-burnable
vals_fill <- !is.na(values(veg_map_unburn))
values(base_layer)[vals_fill] <- 1

levels(nonburnable_layer) <- data.frame(
  value = 1,
  desc = "No quemable"
)

burnable_layer <- na_raster
vals_fill <- is.na(values(veg_map_unburn))
values(burnable_layer)[vals_fill] <- 2
levels(burnable_layer) <- data.frame(
  value = 2,
  desc = "Quemable"
)

# Simulate fire for animation  -------------------------------------------

# fitted sample 2, seed = 2

# Get parameter
column <- 2; ss <- 2
particle <- ranef_fit[, column]

# Tidy parameters
coef_intercepts <- rep(particle[1], n_veg)
coef_nd <- particle[2:3]
coef_terrain <- particle[4:(n_coef-1)]
steps <- particle[n_coef]

set.seed(ss)

anim <- simulate_fire_animate(
  layer_vegetation = land[, , "veg"],
  layer_nd = land[, , nd_variables],
  layer_terrain = land[, , terrain_variables],
  coef_intercepts = coef_intercepts,
  coef_nd = coef_nd,
  coef_terrain = coef_terrain,
  ignition_cells = spread_data$ig_rowcol, # created with zero-indexing
  upper_limit = upper_limit,
  steps = steps
)

plot(rast_from_mat(anim, balcon_rast))


## Create animation layers

simrast_anim <- rast_from_mat(anim, balcon_rast)
values(simrast_anim)[values(simrast_anim) < 1] <- NA
plot(simrast_anim)

n_frames <- 21
nsteps <- max(values(simrast_anim), na.rm = TRUE)
steps <- seq(2, nsteps, length.out = n_frames + 1) %>% as.integer

sim_steps <- lapply(
  1:n_frames,
  function(f) {
    # f = 1
    front <- simrast_anim
    values(front)[values(front) != steps[f]] <- NA
    front <- as.factor(front)
    levels(front) <- data.frame(value = steps[f], front = "Frente del incendio")
    
    back <- simrast_anim
    values(back)[values(back) < steps[f]] <- 1
    values(back)[values(back) >= steps[f]] <- NA
    back <- as.factor(back)
    levels(back) <- data.frame(value = 1, back = "Quemado")
    
    bb <- c(front, back)
    
    return(bb)
  }
)

# define colors for fire
colors <- c(
  "black",
  viridis(1, begin = 1, option = "A"),
  viridis(1, begin = 0.2, option = "A"),
  viridis(1, begin = 0.6, option = "A")
)

# Save frames
for(s in 1:n_frames) {
  
  # s = 1
  rr <- base_layer
  vals_fill <- !is.na(values(sim_steps[[s]][["back"]]))
  values(rr)[vals_fill] <- 3 # burned
  vals_fill <- !is.na(values(sim_steps[[s]][["front"]]))
  values(rr)[vals_fill] <- 4 # front
  
  levels(rr) <- data.frame(
    value = 1:4,
    desc = c("No quemable", "Quemable", "Quemado", "Frente")
  )
  
  advance <- ggplot() +
    geom_spatraster(data = rr) +
    scale_fill_manual(
      values = colors,
      na.value = "transparent",
      na.translate = F
    ) +
    
    # ignition point
    geom_spatvector(
      data = igpoint,
      fill = "white", size = 2, color = "black",
      shape = 21, stroke = 0.8
    ) +
    
    annotation_scale(
      height = unit(1, "mm"),
      text_col = "black", location = "tr", style = "ticks",
      line_col = "black"
    ) +
    
    theme(
      legend.title = element_blank(),#element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
      legend.position = "right",
      axis.text = element_text(size = 8),
      # axis.text.x = element_text(angle = 50, hjust = 1),
      plot.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.spacing.y = unit(0.1, 'mm'),
      legend.margin = margin(),
      axis.ticks = element_line(color = "gray60")
    ) +
    
    scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                       labels = function(x) sprintf("%.2f°", x)) +
    scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                       labels = function(x) sprintf("%.2f°", x)) +
    
    ggtitle("Incendio simulado")
  # advance
  
  fn <- file.path(
    "fire regime simulations", "figures_defensa", "animation_frames",
    paste("spread_anim_", stringr::str_pad(s, 2, pad = "0"), ".png", sep = "")
  )
  ggsave(
    fn, plot = advance, width = 12, height = 9, units = "cm", bg = "white",
    dpi = 300
  )
}


# Overlap maps, animated simulation ---------------------------------------

sim_mat <- anim
sim_mat[sim_mat > 0] <- 1

sim <- simrast_anim
values(sim)[values(simrast_anim) > 0] <- 1
sim_poly <- terra::as.polygons(sim, dissolve = TRUE, na.rm = TRUE)
sim_poly$type <- "Incendio simulado"

obs <- rast_from_mat(spread_data$burned_layer, balcon_rast)
values(obs)[values(obs) == 0] <- NA
obs_poly <- terra::as.polygons(obs, dissolve = TRUE, na.rm = TRUE)
obs_poly$type <- "Incendio real"

# Count intersection and union
vv <- cbind(values(sim), values(obs))
unn <- apply(vv, 1, sum, na.rm = T)
unn <- na.omit(unn)
size_union <- sum(unn >= 1)
size_inter <- sum(unn == 2)
size_inter / size_union # 0.83

## Landscape layer (burnable and not)
rr <- base_layer
levels(rr) <- data.frame(
  value = 1:4,
  desc = c("No quemable", "Quemable")
)

colors_rr <- c(
  "black",
  viridis(1, begin = 1, option = "A")
)

## Both fires together
fires_vec <- rbind(obs_poly, sim_poly)
fires_vec$type <- factor(fires_vec$type, unique(fires_vec$type))
colors_fires <- c(
  viridis(1, begin = 0.6, option = "A"),
  viridis(1, begin = 0.2, option = "A")
)

p1 <- ggplot() +
  geom_spatraster(data = rr) +
  scale_fill_manual(
    values = colors_rr,
    na.value = "transparent",
    na.translate = F
  ) +
  
  ggnewscale::new_scale_fill() +
  geom_spatvector(
    data = fires_vec, 
    mapping = aes(fill = type, alpha = type)
  ) +
  scale_fill_manual(values = colors_fires) +
  scale_alpha_manual(values = c(1, 0.6)) +
  
  # ignition point
  geom_spatvector(
    data = igpoint,
    fill = "white", size = 2, color = "black",
    shape = 21, stroke = 0.8
  ) +
  
  annotation_scale(
    height = unit(1, "mm"),
    text_col = "black", location = "tr", style = "ticks",
    line_col = "black"
  ) +
  
  theme(
    legend.title = element_blank(),#element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
    legend.position = "right",
    axis.text = element_text(size = 8),
    # axis.text.x = element_text(angle = 50, hjust = 1),
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.spacing.y = unit(0.1, 'mm'),
    legend.margin = margin(),
    axis.ticks = element_line(color = "gray60")
  ) +
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) + 
  ggtitle("Superposición")
p1


# intersection
inter_poly <- terra::intersect(sim_poly, obs_poly)

p2 <- ggplot() +
  geom_spatraster(data = rr) +
  scale_fill_manual(
    values = colors_rr,
    na.value = "transparent",
    na.translate = F
  ) +
  
  # intersection (duplicate)
  geom_spatvector(data = inter_poly, fill = colors_fires[1], alpha = 1) +
  geom_spatvector(data = inter_poly, fill = colors_fires[2], alpha = 0.6) + 
  
  # ignition point
  geom_spatvector(
    data = igpoint,
    fill = "white", size = 2, color = "black",
    shape = 21, stroke = 0.8
  ) +
  
  theme(
    legend.title = element_blank(),#element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.spacing.y = unit(0.1, 'mm'),
    legend.margin = margin()
  ) +
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) + 
  ggtitle("Intersección")
p2

# union
union_mat <- sim_mat + spread_data$burned_layer
union_mat[union_mat > 0] <- 1
union_mat[union_mat == 0] <- NA
union_rast <- rast_from_mat(union_mat, balcon_rast)
union_poly <- terra::as.polygons(union_rast, dissolve = TRUE, na.rm = TRUE)

# Color union
color_union <- viridis(1, begin = 0.4, option = "A")

p3 <- ggplot() +
  geom_spatraster(data = rr) +
  scale_fill_manual(
    values = colors_rr,
    na.value = "transparent",
    na.translate = F
  ) +
  
  # union
  geom_spatvector(data = union_poly, fill = color_union, alpha = 0.8) +
  
  # ignition point
  geom_spatvector(
    data = igpoint,
    fill = "white", size = 2, color = "black",
    shape = 21, stroke = 0.8
  ) +
  
  theme(
    legend.title = element_blank(),#element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.spacing.y = unit(0.1, 'mm'),
    legend.margin = margin()
  ) +
  
  scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) +
  scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                     labels = function(x) sprintf("%.2f°", x)) + 
  ggtitle("Unión")
p3

merged <- 
(p1 + p2 + p3) + 
  plot_layout(nrow = 1, byrow = TRUE, guides = "collect") & 
  theme(legend.position = "right")

ggsave(
  "fire regime simulations/figures_defensa/overlap_computation.png",
  bg = "white", plot = merged, width = 25, height = 9, units = "cm"
)


# A few simulations with varying overlap ---------------------------------

# fitted
fire1 <- burn_fitted(2, 2, T)

# simulated:
# very nice
# bigger: 21 23 42
# smaller: 15 32
fire2 <- burn_simulate(21, 21, T)
fire3 <- burn_simulate(15, 15, T)
fire4 <- burn_simulate(23, 23, T)
fire5 <- burn_simulate(42, 42, T)
fire6 <- burn_simulate(32, 32, T)

flist <- list(
  fire1, fire2, fire3, fire4, fire5, fire6
)

pp <- vector("list", 6)

for (i in 1:6) {
  
  ovround <- round(flist[[i]]$overlap, 3)
  
  pp[[i]] <- 
    ggplot() +
    geom_spatraster(data = rr) +
    scale_fill_manual(
      values = colors_rr,
      na.value = "transparent",
      na.translate = F
    ) +
    
    # observed and simulated fires
    geom_spatvector(data = obs_poly, fill = colors_fires[1], alpha = 1) +
    geom_spatvector(data = flist[[i]], fill = colors_fires[2], alpha = 0.6) +
    
    # ignition point
    geom_spatvector(
      data = igpoint,
      fill = "white", size = 2, color = "black",
      shape = 21, stroke = 0.8
    ) +
    
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 14, hjust = 0.5, vjust = -0.4)
    ) +
    
    scale_x_continuous(breaks = seq(-71.44, -71.40, by = 0.04), expand = c(0, 0),
                       labels = function(x) sprintf("%.2f°", x)) +
    scale_y_continuous(breaks = seq(-41.20, -41.16, by = 0.04), expand = c(0, 0),
                       labels = function(x) sprintf("%.2f°", x)) + 
    ggtitle(ovround)
  # pp[[i]]
}

mm <- (
  pp[[1]] + pp[[2]] + pp[[3]] + pp[[4]] + pp[[5]] + pp[[6]]
) + 
  plot_layout(nrow = 2, byrow = T)

ggsave(
  "fire regime simulations/figures_defensa/overlap_examples.png",
  bg = "white", plot = mm, width = 25, height = 18, units = "cm"
)
