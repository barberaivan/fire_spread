# Animated maps to show

options(scipen = 999) # turn off scientific notation

library(terra)
library(tidyverse)
library(Rcpp)
library(gganimate)

library(stringr) # add leading 0 to file names
library(magick)

theme_set(theme_bw())

sourceCpp("spread_functions.cpp")
source("spread_functions.R") # for land_cube()


# data --------------------------------------------------------------------

# land <- readRDS(file.path("data", "focal fires data", 
#                           "landscapes_ig-known_non-steppe", "1999_25j.rds"))
# land_raster <- rast(file.path("data", "focal fires data", "raw data from GEE",
#                               "fire_data_raw_1999_25j.tif"))

land <- readRDS(file.path("data", "focal fires data", 
                          "landscapes_ig-known_non-steppe", "1999_28.rds"))
land_raster <- rast(file.path("data", "focal fires data", "raw data from GEE",
                              "fire_data_raw_1999_28.tif"))

n_veg_types <- 4
n_terrain <- 4

vegetation_names <- c("shrubland", "subalpine", "wet", "dry")
terrain_names <- c("northing", "elev", "wind", "slope")
par_names <- c(vegetation_names, terrain_names)

# prior simulator for graphical prior predictive checks
prior_sim <- function(mu_int = 0, sd_int = 20, r_01 = 0.05, r_z = 0.15,
                      n_veg_types = 4) {
  
  intercepts <- rnorm(n_veg_types, mu_int, sd_int)
  names(intercepts) <- vegetation_names
  
  slopes <- c(
    "northing" = rexp(1, r_01),                 # positive 
    "wind" = rexp(1, r_01),                   # positive
    "elev" = (-1) * rexp(1, r_z),        # negative
    "slope" = rexp(1, r_01)                   # positive
  )
  
  return(c(intercepts, slopes))
}


# vegetation map ----------------------------------------------------------

# Recode vegetation and lakes
# unique vegetation classes:
# shrubland, subalpine, wet, dry  // add high andean and lakes separately
#    0           1       2    3
base_mat <- land$vegetation
base_mat[land$vegetation == 99] <- 4 # high andean
# fill lakes
lower_elev <- min(land$terrain[, , "elev"])
add_elev <- diff(range(land$terrain[, , "elev"])) * 0.2
elev_lim <- lower_elev + add_elev
base_mat[(land$vegetation == 99) & (land$terrain[, , "elev"] < elev_lim)] <- 5 # lakes

pix_levels <- c("matorral", "lenga", "coihue", "ciprés", 
                "altoandino", "lagos", 
                "frente del incendio", "quemado")


# simulate ----------------------------------------------------------------

# use this to searh example parameters
# params <- prior_sim(mu_int = -3, sd_int = 0.2, r_01 = 1, r_z = 1)

fire_types <- c("small", "mid", "huge")
# example particles:
#              shrubland  subalpine        wet        dry   northing       wind       elev      slope 
params_small <- c(-4, -4, -4, -4,  0.5207415,  0.6087135, -0.5,  1.0558857)
params_mid <- c(-2.9289697, -3.3436993, -3.2605035, -2.8138192,  0.5207415,  0.6087135, -0.5,  2)
params_huge <- c(10, 10, 10, 10,  0.5207415,  0.6087135, -1.9038349,  1.0558857)

params_mat <- cbind(params_small, params_mid, params_huge)
colnames(params_mat) <- fire_types

seeds <- c(7957, 9984, 341)
names(seeds) <- fire_types

# loop to create an animation by fire type

# for(f in fire_types) {
  f = "small"
  
  s <- seeds[f]
  # s <- round(runif(1, 1, 10000))
  set.seed(s)
  fire1 <- simulate_fire_animate(
    terrain = land$terrain, 
    vegetation = land$vegetation,
    ignition_cells = land$ig_rowcol,
    coef = params_mat[, f],
    n_veg_types = n_veg_types,
    upper_limit = 0.5
  )
  f1 <- rast_from_mat(fire1, land_raster)
  plot(f1)
  # create array with vegetation and fire 
  steps <- max(fire1) + 1 # add a step with cells turned off
  fire1_arr <- array(0, dim = c(nrow(land$vegetation), ncol(land$vegetation),
                                steps),
                     dimnames = list(
                       row = 1:nrow(land$vegetation),
                       col = 1:ncol(land$vegetation), 
                       step = 1:steps
                     ))
  
  
  # include burned pixels in each step
  for(s in 1:steps) {
    place_holder <- base_mat
    place_holder[(fire1 == s)] <- 6 # fire front
    place_holder[(fire1 < s) & (fire1 > 0)] <- 7 # burned
    fire1_arr[, , s] <- place_holder
  }
  
  # make many-layers raster
  fire_list <- lapply(1:steps, function(s) rast_from_mat(fire1_arr[, , s], 
                                                         land_raster))
  fire_rast <- do.call("c", fire_list)
  names(fire_rast) <- 1:steps
  
  # get values to make df
  vv <- values(fire_rast)
  
  # show every k stpes:
  k <- 2
  keep <- seq(1, steps, by = k)
  if(!(max(steps) %in% keep)) keep <- c(keep, max(steps))
  
  coords <- xyFromCell(fire_rast, 1:ncell(fire_rast))
  rast_wide <- as.data.frame(cbind(coords, vv[, keep]))
  rast_long <- pivot_longer(rast_wide, all_of(which(colnames(rast_wide) %in% as.character(keep))),
                            names_to = "step", values_to = "val")
  
  
  rast_long$class <- factor(pix_levels[rast_long$val + 1], 
                            levels = pix_levels)
  rast_long$step <- as.numeric(rast_long$step)
  
  # single step plot
  map_colors <- c("darkolivegreen4", "darkgreen", "forestgreen", "#74c476",
                  "gray", "deepskyblue", 
                  "red", "black")
  for(i in 1:length(keep)) {
    print(paste("step", keep[i], "out of", steps))
    
    p1 <- ggplot(rast_long[rast_long$step == keep[i], ], 
                 aes(x = x, y = y, fill = class)) + 
      geom_tile() + 
      coord_fixed() + 
      scale_fill_manual(values = map_colors,
                        limits = levels(rast_long$class)) +
      theme_minimal() + 
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.title = element_blank(),
            legend.spacing.y = unit(0.5, 'mm')) +
      guides(fill = guide_legend(byrow = TRUE))
    # p1 
    
    nn <- str_pad(i, 3, pad = "0")
    
    ggsave(file.path("files", "animations", 
                     paste("fire_img_", f, "_", nn, ".jpeg", sep = "")), 
           plot = p1 + theme(legend.position = "none"), 
           width = 4, height = 5, units = "cm", dpi = 200)
  }
  
  # read images
  imgs <- list.files(file.path("files", "animations"), 
                     full.names = TRUE, pattern = paste("img", f, sep = "_"))
  img_list <- lapply(imgs, image_read)
  # join 
  img_joined <- image_join(img_list)
  # animate 
  img_animated <- image_animate(img_joined, fps = 10)
  # write
  image_write(image = img_animated,
              path = paste("files/animations/", "fire_animation_", f, 
                           ".gif", sep = ""))
# } # end looop over gifs


# plot with legend
rast_leg <- rast_long[rast_long$step == 1, ]
rast_leg$class[rast_leg$class == "frente del incendio"] <- "lenga"
p_labs1 <- ggplot(rast_leg, 
                  aes(x = x, y = y, fill = class)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_fill_manual(values = map_colors[1:6],
                    limits = levels(rast_long$class)[1:6]) +
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.spacing.y = unit(4, 'mm'),
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.margin = unit(rep(5, 4), "mm")) +
  guides(fill = guide_legend(byrow = TRUE))
# ggsave("files/animations/veg_map_1.jpeg", plot = p_labs1,
#        width = 14, height = 12, dpi = 300, units = "cm") 
  
p_labs2 <- ggplot(rast_leg, 
                  aes(x = x, y = y, fill = class)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_fill_manual(values = map_colors,
                    limits = levels(rast_long$class)) +
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.spacing.y = unit(4, 'mm'),
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.margin = unit(rep(5, 4), "mm")) +
  guides(fill = guide_legend(byrow = TRUE))
ggsave("files/animations/veg_map_2.jpeg", plot = p_labs2,
       width = 14, height = 12, dpi = 300, units = "cm") 

# The same with the real fire
tmp <- base_mat
tmp[land$burned_layer == 1] <- 7 # burned
tmp_rast <- rast_from_mat(tmp, land_raster)
tmp_df <- as.data.frame(cbind(coords, values(tmp_rast)))
colnames(tmp_df) <- c("x", "y", "val")
tmp_df$class <- factor(pix_levels[tmp_df$val + 1], levels = pix_levels)

p_labs3 <- ggplot(tmp_df, aes(x = x, y = y, fill = class)) + 
  geom_tile() + 
  coord_fixed() + 
  scale_fill_manual(values = map_colors[-7],
                    limits = levels(rast_long$class)[-7]) +
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.spacing.y = unit(4, 'mm'),
        plot.background = element_rect(fill = "white", colour = "white"),
        plot.margin = unit(rep(5, 4), "mm")) +
  guides(fill = guide_legend(byrow = TRUE))
p_labs3
# ggsave("files/animations/veg_map_3.jpeg", plot = p_labs3,
#        width = 14, height = 12, dpi = 300, units = "cm") 

# merge veg and fire plots
plots_veg_real <- egg::ggarrange(
  p_labs1 + theme(legend.position = "none"), 
  p_labs3,
  ncol = 2
)

ggsave("files/animations/veg_map_1.jpeg", plot = plots_veg_real,
     width = 22, height = 12, dpi = 300, units = "cm")

# gganimate ---------------------------------------------------------------

# gganimate option didn't work well with large fires, so I moved to 
# magick
  
  # many steps
  anim <- ggplot(rast_long, 
                 aes(x = x, y = y, fill = class, group = step)) + 
    geom_tile() + 
    coord_fixed() + 
    scale_fill_manual(values = map_colors,
                      limits = pix_levels)  +
    theme_minimal() + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "none") +
    transition_manual(step)
  
  # anim_save()
  # anim
  gif_name <- paste("fire_", f, ".gif", sep = "")
  
  gganimate::animate(
    anim,
    nframes = length(keep), 10, #
    fps = 5,
    width = 6, height = 7, units = "cm", res = 100,
    renderer = gifski_renderer(gif_name)
  )  
  