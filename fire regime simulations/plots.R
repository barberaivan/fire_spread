library(terra)
library(ggplot2)
library(tidyterra)
library(ggspatial)
library(viridis)
library(patchwork)
library(brms)
library(bayestestR) # hdi interval

library(extrafont)
# extrafont::font_import()

theme_set(theme_classic()) # for non maps
theme_set(theme_minimal()) # for maps

# Figure size settings ----------------------------------------------------

a4h <- 29.7
a4w <- 21.0
margins <- 2.5
fig_width_max <- a4w - margins * 2
fig_height_max <- a4h - margins * 2

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
    
    axis.text = element_text(size = 7),
    axis.title = element_blank(),
    
    plot.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    
    strip.text = element_text(size = 9),
    strip.background = element_rect(fill = "white", color = "white")
  )
}


# Functions to summarise ---------------------------------------------------

hdi_lower <- function(x) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = 0.95)[["CI_low"]])
}

hdi_upper <- function(x) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = 0.95)[["CI_high"]])
}

hdmean <- function(x, ci = 0.95, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(mean(x), ci$CI_low, ci$CI_high)
  names(result) <- paste(rep(name, 3), c("mean", "lower", "upper"), sep = "_")
  result
}

# Load files --------------------------------------------------------------

# Burn probability maps
source_dir <- file.path("files", "fire_regime_simulation_FWIZ")
files <- list.files(source_dir, pattern = glob2rx("*burn_prob*.tif*"))
files <- c(files[9], files[1:8])

rlist <- lapply(1:9, function(i) {
  rast(file.path(source_dir, files[i]))
})

exp_names <- c("1999-2022",
               "2040-2049 SSP1 2.6", "2040-2049 SSP2 4.5",
               "2040-2049 SSP3 7.0", "2040-2049 SSP5 8.5",
               "2090-2099 SSP1 2.6", "2090-2099 SSP2 4.5",
               "2090-2099 SSP3 7.0", "2090-2099 SSP5 8.5")

sc_names <- c("SSP1-2.6", "SSP2-4.5",
              "SSP3-7.0", "SSP5-8.5")

bmap <- do.call("c", rlist)
bmap <- bmap * 100 # prob as percent
names(bmap) <- exp_names

# PNNH layer, with lakes and non-burnable
pnnh_rast_full <- rast(file.path("data", "pnnh_images", "pnnh_data_30m_buff_10000_std.tif"))
veg_map0 <- pnnh_rast_full$veg # 11 classes

# vegetation transform data:
dveg <- readxl::read_excel("/home/ivan/Insync/Mapa vegetación WWF - Lara et al. 1999/clases de vegetacion y equivalencias.xlsx",
                           sheet = "Sheet2")

# Aplica la reclasificación
veg_map <- as.factor(veg_map0)

# Convierte el raster reclasificado en un factor con etiquetas
levels(veg_map) <- list(data.frame(
  ID = 1:11, 
  Class = c("Bosque húmedo", 
            "Bosque subalpino", 
            "Bosque subalpino",
            "Bosque seco",
            "Matorral",
            "Pastizal",
            "Plantación",
            "Áreas urbanas",
            "Lagos",
            "Altoandino",
            "Altoandino")
))
# plot(veg_map)

bmap2 <- resample(bmap, veg_map)
bmap3 <- mask(bmap2, veg_map0, maskvalue = 8:11) # necesita nros, use raw layer

# Nahuel Huapi National Park (strict)
pnnh <- vect(file.path("data", "protected_areas", "apn_limites.shp"))
pnnh_latlong <- pnnh[pnnh$nombre == "Nahuel Huapi", ]
pnnh <- project(pnnh_latlong, "EPSG:5343")

# Buffer around PNNH (10 km buffer, to simulate ignitions; posgar 2007)
pnnh_buff <- vect(file.path("data", "protected_areas", "pnnh_buff_10000.shp"))

# PNNH map with predicted burn prob by model:
burn_prob_raster <- rast(file.path("data", "pnnh_images", 
                                   "pnnh_data_120m_buff_10000_ig-esc-spread-prob_FWIZ.tiff"))
burn_prob_models <- burn_prob_raster[[c("igprob_h", "igprob_l", 
                                        "escprob", "spreadprob")]] * 100



# Load fires and ignition points
fires <- vect("/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp")
igpoints <- vect("/home/ivan/Insync/Fire spread modelling/ignition_data/ignition_points_pnnh_bari-kitzberger.shp")
igpoints$Causa <- factor(igpoints$cause,
                         levels = c("human", "lightning", "unknown"),
                         labels = c("Humanos", "Rayos", "Desconocida"))

fires_pnnh <- fires[relate(fires, pnnh_latlong, relation = "intersects")]
fires_pnnh$Incendios <- "Quemado"
# plot(pnnh_latlong)
# plot(fires_pnnh, add = T, col = "green")
# plot(igpoints, add = T, col = "red")

# Data and predictors map --------------------------------------------------

veg_map_burn <- mask(veg_map, veg_map0, maskvalue = 8:11)
veg_map_unburn <- mask(veg_map, veg_map0, maskvalue = 8:11, inverse = T)
ndvi_lay <- mask(pnnh_rast_full$ndvi, veg_map0, maskvalue = 8:11)
elev_lyr <- mask(pnnh_rast_full$elevation,
                 veg_map0, maskvalue = 9)

# settings gráficos
bp_vir <- "F"
bp_begin <- 1
bp_end <- 0.1
bp_name <- "Probabilidad\nde quema (%)"
maxcell <- 100000

park_contour <- "black"
park_lwd <- 0.2


## National park with fires and points

datamap <-  ggplot() +
  
  # National park
  geom_spatvector(data = pnnh, fill = "white", color = park_contour, 
                  linewidth = 0.4) + 
  
  # Fires
  geom_spatvector(data = fires_pnnh, color = viridis(1, option = "C"),
                  mapping = aes(fill = Incendios), linewidth = 0.1) + 
  scale_fill_viridis(option = "C", discrete = T, name = "") +
  
  # Ignition points
  geom_spatvector(data = igpoints, size = 1.3, stroke = 0.3,
                  mapping = aes(color = Causa, shape = Causa)) +
  scale_color_viridis(option = "C", discrete = T, begin = 0.3, end = 0.8,
                      name = "Focos (causa)") +
  scale_shape_manual(values = 21:23, name = "Focos (causa)") +
  
  scale_x_continuous(breaks = seq(-72, -71, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  scale_y_continuous(breaks = seq(-41.5, -40.5, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() +
  theme(panel.grid.major = element_line(linewidth = 0.1, color = "gray20"),
        axis.text.x = element_blank(),
        axis.ticks = element_line(linewidth = 0.1, color = "gray20")) +
  ggspatial::annotation_scale(
    location = "tr", height = unit(0.7, "mm"),
    bar_cols = c("grey60", "white"),
    text_col = "gray20",
    text_cex = 0.5,
    line_width = 0.2
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
    width = unit(3, "mm"),
    height = unit(3, "mm")
  )
datamap

## Argentina MAP
arg_sf <- rnaturalearth::ne_countries(scale = "large", country = "Argentina", 
                                      returnclass = "sf")
arg_map <- ggplot() +
  geom_sf(data = arg_sf, fill = "white", color = "gray40") +
  geom_sf(data = pnnh, fill = "forestgreen", alpha = 1) +
  scale_x_continuous(breaks = seq(-70, -60, by = 10), expand = c(0, 0),
                     labels = function(x) sprintf("%+d°", x)) +
  scale_y_continuous(breaks = seq(-25, -55, by = -10), expand = c(0, 0),
                     labels = function(y) sprintf("%+d°", y)) +
  coord_sf() +
  map_theme() +
  theme_minimal()
arg_map

## Vegetation map (burnable areas)


veg_burn <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = veg_map_burn, maxcell = maxcell) +
  scale_fill_viridis(discrete = T, option = "D", 
                     na.value = "transparent", na.translate = F,
                     begin = 0, end = 1,
                     name = "Vegetación") +
  
  # National park
  geom_spatvector(data = pnnh, fill = NA, color = park_contour, 
                  linewidth = 0.4) + 
  
  scale_x_continuous(breaks = seq(-72, -71, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  scale_y_continuous(breaks = seq(-41.5, -40.5, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(linewidth = 0.1, color = "gray20"),
        axis.text.x = element_blank())
veg_burn

## Vegetation map (non-burnable areas)

veg_unburn <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = veg_map_unburn, maxcell = maxcell) +
  scale_fill_manual(values = viridis(3, option = "E", begin = 0.15, end = 0.85)[c(3, 1, 2)],
                    na.value = "transparent", na.translate = F,
                    name = "Superficie\nno quemable") +
  
  # National park
  geom_spatvector(data = pnnh, fill = NA, color = park_contour, 
                  linewidth = 0.4) + 
  
  scale_x_continuous(breaks = seq(-72, -71, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  scale_y_continuous(breaks = seq(-41.5, -40.5, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.text = element_blank())
veg_unburn

## Vegetation map (non-burnable areas)



ndvimap <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = ndvi_lay, maxcell = maxcell) +
  scale_fill_viridis(option = "G", na.value = "transparent",
                     begin = 0, end = 1, direction = -1,
                     name = "NDVI") +
  
  # National park
  geom_spatvector(data = pnnh, fill = NA, color = "white", 
                  linewidth = 0.4) + 
  
  scale_x_continuous(breaks = seq(-72, -71, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  scale_y_continuous(breaks = seq(-41.5, -40.5, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(linewidth = 0.1))
ndvimap

elev <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = elev_lyr, maxcell = maxcell) +
  scale_fill_viridis(option = "A", na.value = "transparent",
                     begin = 0, end = 1,
                     name = "Altitud\n(m s.n.m.)") +
  
  # National park
  geom_spatvector(data = pnnh, fill = NA, color = park_contour, 
                  linewidth = 0.4) + 
  
  scale_x_continuous(breaks = seq(-72, -71, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  scale_y_continuous(breaks = seq(-41.5, -40.5, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(linewidth = 0.1),
        axis.text.y = element_blank())
elev 

preds <- 
  (
    arg_map + theme(axis.text = element_text(size = 7)) +
    datamap + theme(axis.text = element_blank()) +
    veg_burn + 
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank()) +
    veg_unburn + 
      theme(axis.ticks = element_blank(),
            axis.text = element_blank()) +
    ndvimap + 
    elev +
      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank())
  ) +
  plot_layout(nrow = 3, byrow = T,
              heights = c(2, 2, 2),
              widths = c(0.7, 0.7)) &
  theme(legend.key.width = unit(2.5, "mm"),
        legend.key.height = unit(2.5, "mm"),
        legend.box.just = "left",
        legend.box.margin = margin(0, 0, 0, -1, unit = "mm"),
        legend.justification = "left",
        
        legend.margin = margin(0, 0, 0, 0, unit = "mm"), # Márgenes alrededor de la leyenda
        legend.spacing = unit(1, "mm"), 
        legend.spacing.y = unit(1.3, "mm"), # Espaciado entre filas de la leyenda
        
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        
        panel.spacing.x = unit(0, "mm"),
        plot.margin = margin(1, 1, 1, 1, unit = "mm"))
  

ggsave("fire regime simulations/figures/pnnh_landscape.pdf",
       plot = preds, 
       width = fig_width_max * 0.8, height = 18, units = "cm")
ggsave("fire regime simulations/figures/pnnh_landscape.png",
       plot = preds, 
       width = fig_width_max * 0.8, height = 18, units = "cm", bg = "white")


# Burn probability: separate models and current climate -------------------

# settings gráficos
bp_vir <- "F"
bp_begin <- 1
bp_end <- 0.1
bp_name <- "Probabilidad\nde quema (%)"
maxcell <- 100000

park_contour <- "black"
park_lwd <- 0.2

# Burn prob anual
bp_modern <- ggplot() +
  
  # Burn probability map
  geom_spatraster(data = bmap3[[1]], maxcell = maxcell) +
  scale_fill_viridis(option = bp_vir, na.value = "transparent",
                     begin = bp_begin, end = bp_end,
                     name = "%") +
  
  # National park
  geom_spatvector(data = pnnh, fill = NA, color = park_contour, 
                  linewidth = park_lwd) + 
  geom_spatvector(data = pnnh_buff, fill = NA, color = "gray20", 
                  linewidth = park_lwd, alpha = 0.6) + 
  
  scale_x_continuous(breaks = seq(-72, -71, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  scale_y_continuous(breaks = seq(-41.5, -40.5, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
                     # edit coordinates format
  coord_sf() +
  map_theme() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(linewidth = 0.1, color = "gray20"),
        legend.key.width = unit(2, 'mm')) + 
  labs(title = "E. Probabilidad de quema\nanual (1999-2022)") +
  
  ggspatial::annotation_scale(
    location = "tr", height = unit(0.7, "mm"),
    bar_cols = c("grey60", "white"),
    text_col = "black",
    text_cex = 0.5,
    line_width = 0.2
  ) +
  annotation_north_arrow(
    location = "tl", which_north = "true",  
    style = north_arrow_orienteering(
      line_width = 0.2,
      line_col = "black",
      fill = c("grey60", "white"),
      text_col = "transparent",
      text_family = "Arial",
      text_face = NULL,
      text_size = 0,
      text_angle = 0
    ),
    width = unit(3, "mm"),
    height = unit(3, "mm"))

bp_modern


# burn prob separate
ll <- vector("list", 5)
titles <- c(
  "A. Probabilidad relativa de\nignición por humanos",
  "B. Probabilidad relativa de\nignición por rayos",
  "C. Probabilidad de escape",
  "D. Probabilidad de\npropagación"
)

for(m in 1:4) {
  # Burn prob anual
  ll[[m]] <- 
    
    ggplot() +
    
    # Burn probability map
    geom_spatraster(data = burn_prob_models[[m]], maxcell = maxcell) +
    scale_fill_viridis(option = bp_vir, na.value = "transparent",
                       begin = bp_begin, end = bp_end,
                       name = "%") +
    
    # National park
    geom_spatvector(data = pnnh, fill = NA, color = park_contour, 
                    linewidth = park_lwd) + 
    # geom_spatvector(data = pnnh_buff, fill = NA, color = park_contour, 
    #                 linewidth = park_lwd) + 
    
    scale_x_continuous(breaks = seq(-72, -71, by = 0.5), expand = c(0, 0),
                       labels = function(x) sprintf("%.1f°", x)) +
    scale_y_continuous(breaks = seq(-41.5, -40.5, by = 0.5), expand = c(0, 0),
                       labels = function(x) sprintf("%.1f°", x)) +
    # edit coordinates format
    coord_sf() +
    map_theme() + 
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.key.width = unit(2, 'mm')) + 
    labs(title = titles[m]) 
}


bp1 <- ll[[1]] + ll[[2]] + ll[[3]] + ll[[4]] + bp_modern +
  plot_layout(ncol = 2, byrow = T)

ggsave("fire regime simulations/figures/burn_prob_models_modern.pdf",
       plot = bp1, width = fig_width_max * 0.7, height = fig_height_max - 3,
       units = "cm")
ggsave("fire regime simulations/figures/burn_prob_models_modern.png",
       plot = bp1, width = fig_width_max * 0.7, height = fig_height_max - 3,
       units = "cm", bg = "white")


# Burn prob 2040 ----------------------------------------------------------

bp_2040 <- bmap3[[2:5]] 
names(bp_2040) <- sc_names


# df for north arrow and scalebar
raster_df <- as.data.frame(bp_2040, xy = TRUE)
raster_df$lyr <- "SSP1-2.6"

bpmap40 <-
ggplot() +
  
  # Burn probability map
  geom_spatraster(data = bp_2040, maxcell = maxcell) +
  scale_fill_viridis(option = bp_vir, na.value = "transparent",
                     begin = bp_begin, end = bp_end,
                     name = bp_name) +
  
  facet_wrap(~ lyr) +
  
  # National park
  geom_spatvector(data = pnnh, fill = NA, color = park_contour, 
                  linewidth = park_lwd) + 
  geom_spatvector(data = pnnh_buff, fill = NA, color = park_contour, 
                  linewidth = park_lwd) + 
  
  scale_x_continuous(breaks = seq(-72, -71, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  scale_y_continuous(breaks = seq(-41.5, -40.5, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() +
  theme(panel.grid = element_blank(),
        panel.spacing.y = unit(4, "mm"),
        legend.key.width = unit(2, "mm"),
        axis.ticks = element_line(linewidth = 0.1)) +
  
  # scale and north arrow
  ggspatial::annotation_scale(
    data = raster_df, 
    location = "tr", height = unit(0.7, "mm"),
    bar_cols = c("grey60", "white"),
    text_col = "black",
    text_cex = 0.5,
    line_width = 0.2
  ) +
  annotation_north_arrow(
    data = raster_df[1, ], 
    location = "tl", which_north = "true",  
    style = north_arrow_orienteering(
      line_width = 0.2,
      line_col = "black",
      fill = c("grey60", "white"),
      text_col = "transparent",
      text_family = "arial",
      text_face = NULL,
      text_size = 0,
      text_angle = 0
    ),
    width = unit(3, "mm"),
    height = unit(3, "mm")) +
  ggtitle("2040-2049")

bpmap40

ggsave("fire regime simulations/figures/burn_prob_2040.pdf",
       plot = bpmap40, 
       width = 13, height = 18, units = "cm")
ggsave("fire regime simulations/figures/burn_prob_2040.png",
       plot = bpmap40, 
       width = 13, height = 18, units = "cm", bg = "white")


# Burn prob 2090 ----------------------------------------------------------

bp_2090 <- bmap3[[6:9]] 
names(bp_2090) <- sc_names

# df for north arrow and scalebar
raster_df <- as.data.frame(bp_2040, xy = TRUE)
raster_df$lyr <- "SSP1-2.6"

bpmap90 <- 
ggplot() +
  
  # Burn probability map
  geom_spatraster(data = bp_2090, maxcell = maxcell) +
  scale_fill_viridis(option = bp_vir, na.value = "transparent",
                     begin = bp_begin, end = bp_end,
                     name = bp_name) +
  
  facet_wrap(~ lyr) +
  
  # National park
  geom_spatvector(data = pnnh, fill = NA, color = park_contour, 
                  linewidth = park_lwd) + 
  geom_spatvector(data = pnnh_buff, fill = NA, color = park_contour, 
                  linewidth = park_lwd) + 
  
  scale_x_continuous(breaks = seq(-72, -71, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  scale_y_continuous(breaks = seq(-41.5, -40.5, by = 0.5), expand = c(0, 0),
                     labels = function(x) sprintf("%.1f°", x)) +
  # edit coordinates format
  coord_sf() +
  map_theme() +
  theme(panel.grid = element_blank(),
        panel.spacing.y = unit(4, "mm"),
        legend.key.width = unit(2, "mm"),
        axis.ticks = element_line(linewidth = 0.1)) +
  
  # scale and north arrow
  ggspatial::annotation_scale(
    data = raster_df, 
    location = "tr", height = unit(0.7, "mm"),
    bar_cols = c("grey60", "white"),
    text_col = "black",
    text_cex = 0.5,
    line_width = 0.2
  ) +
  annotation_north_arrow(
    data = raster_df[1, ], 
    location = "tl", which_north = "true",  
    style = north_arrow_orienteering(
      line_width = 0.2,
      line_col = "black",
      fill = c("grey60", "white"),
      text_col = "transparent",
      text_family = "arial",
      text_face = NULL,
      text_size = 0,
      text_angle = 0
    ),
    width = unit(3, "mm"),
    height = unit(3, "mm")) +
  
  ggtitle("2090-2099")

bpmap90

ggsave("fire regime simulations/figures/burn_prob_2090.pdf",
       plot = bpmap90, 
       width = 13, height = 18, units = "cm")
ggsave("fire regime simulations/figures/burn_prob_2090.png",
       plot = bpmap90, 
       width = 13, height = 18, units = "cm", bg = "white")



# Burned proportion in PNNH and lightning-burned proportion ---------------

exp_names2 <- c("2040-2049 SSP1 2.6", "2040-2049 SSP2 4.5",
                "2040-2049 SSP3 7.0", "2040-2049 SSP5 8.5",
                "2090-2099 SSP1 2.6", "2090-2099 SSP2 4.5",
                "2090-2099 SSP3 7.0", "2090-2099 SSP5 8.5",
                "1999-2022")

sc_names <- c("SSP1-2.6", "SSP2-4.5",
              "SSP3-7.0", "SSP5-8.5")

tab <- readRDS(file.path("files", "fire_regime_simulation_FWIZ",
                         "burn_prop_distribution_merged.rds"))

ske <- aggregate(burned_count ~ experiment + scenario + decade, 
                 tab, mean)[, 1:3]
ske$land_size <- 100
ske$burned_count <- 100

ske <- ske[c(9, 1:8), ]

ske$experiment2 <- factor(ske$experiment, levels = levels(tab$experiment),
                          labels = exp_names)

sss <- c("Moderno", rep(c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"), 2))
ske$scenario <- factor(sss, levels = unique(sss))
ppp <- c("1999-2022", rep("2040-2049", 4), rep("2090-2099", 4))
ske$period <- factor(ppp, levels = unique(ppp))

# ## Model for burned proportion
# bprop_model <- brm(
#   formula = bf(burned_count | trials(land_size) ~ experiment, 
#                phi ~ experiment),  
#   family = beta_binomial(), 
#   data = tab,
#   chains = 4, iter = 2000, warmup = 1000, refresh = 50, cores = 4
# )
# saveRDS(bprop_model, file.path("files", "fire_regime_simulation_FWIZ",
#                                "burn_prop_model.rds"))


# ## Model for lightning proportion
# lprop_model <- brm(
#   formula = bf(light_count | trials(burned_count) ~ experiment, 
#                phi ~ experiment),  
#   family = beta_binomial(), 
#   data = tab[tab$burned_count > 0, ],
#   chains = 4, iter = 2000, warmup = 1000, refresh = 50, cores = 4
# )
# saveRDS(lprop_model, file.path("files", "fire_regime_simulation_FWIZ",
#                                "burn_prop_model_lightning.rds"))


## Predicciones


bprop_model <- readRDS(file.path("files", "fire_regime_simulation_FWIZ",
                                 "burn_prop_model.rds"))
lprop_model <- readRDS(file.path("files", "fire_regime_simulation_FWIZ",
                                 "burn_prop_model_lightning.rds"))
 
pfit_samples <- fitted(bprop_model, newdata = ske, summary = F) |> t()
bprop <- cbind(
  ske,
  apply(pfit_samples, 1, hdmean, name = "p") |> t() |> as.data.frame()
)
bprop$response <- "Proporción quemada anual (%)"

ww <- 0.9
ggplot(bprop, aes(x = period, y = p_mean, ymin = p_lower, ymax = p_upper,
                  color = scenario, fill = scenario)) +
  geom_bar(stat = "identity", 
           position = position_dodge2(width = ww, preserve = "single"),
           size = 0.3, alpha = 0.6, width = ww) +
  geom_linerange(position = position_dodge2(width = ww, preserve = "single")) +
  scale_color_viridis(option = "C", discrete = T, end = 0.8, name = "Escenario\nclimático") +
  scale_fill_viridis(option = "C", discrete = T, end = 0.8, name = "Escenario\nclimático",
                     alpha = 0.5) + 
  xlab("Período") + 
  ylab("Proporción quemada anual (%)") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 3, 0.5)) +
  nice_theme() +
  theme(panel.grid.major.y = element_line(linewidth = 0.25, color = "gray80"))


ggsave("fire regime simulations/figures/burn_prob_mean_comparison.pdf",
       plot = last_plot(), 
       width = 14, height = 8, units = "cm")
ggsave("fire regime simulations/figures/burn_prob_mean_comparison.png",
       plot = last_plot(), 
       width = 14, height = 8, units = "cm", bg = "white")

write.csv(bprop, "fire regime simulations/figures/burn_prob_means.csv", 
          row.names = F)


# Lightning
lfit_samples <- fitted(lprop_model, ske, summary = F) |> t()
lprop <- cbind(
  ske,
  apply(lfit_samples, 1, hdmean, name = "p") |> t() |> as.data.frame()
)
lprop$response <- "Fracción quemada por rayos (%)"

ggplot(lprop, aes(x = scenario, y = p_mean, ymin = p_lower, ymax = p_upper,
                  color = period, fill = period, shape = period)) +
  geom_linerange() +
  geom_point(size = 3, alpha = 0.2) +
  ylab("Fracción quemada por rayos (%)") + 
  xlab("Escenario climático")


# BOTH

both_prop <- rbind(bprop, lprop)
both_prop$response <- factor(
  both_prop$response,
  levels = c("Proporción quemada anual (%)",
             "Fracción quemada por rayos (%)"),
  labels = c("A. Proporción quemada anual (%)",
             "B. Fracción quemada por rayos (%)")
)

ww <- 0.8
ggplot(both_prop, aes(x = period, y = p_mean, ymin = p_lower, ymax = p_upper,
                      color = scenario, fill = scenario)) +
  geom_bar(stat = "identity", 
           position = position_dodge2(width = ww, preserve = "single"),
           size = 0.3, alpha = 0.6, width = ww) +
  geom_linerange(position = position_dodge2(width = ww, preserve = "single")) +
  scale_color_viridis(option = "C", discrete = T, end = 0.8, name = "Escenario\nclimático") +
  scale_fill_viridis(option = "C", discrete = T, end = 0.8, name = "Escenario\nclimático",
                     alpha = 0.5) + 
  xlab("Período") + 
  facet_wrap(vars(response), nrow = 2, axes = "all", strip.position = "left",
             axis.labels = "margins", scales = "free_y") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title.y = element_blank(),
        strip.placement = "outside",
        panel.spacing.y = unit(4, "mm")) +
  nice_theme()




# TAREAS ------------------------------------------------------------------

# volver a correr los mapas