# The first part is only prior predictive checks, and the second is to explore
# discrepancy functions.
# the prior check ir repeated in the <prior predictive checks.R> script,
# and should be deleted later


# Packages and data -------------------------------------------------------

library(terra)
library(tidyverse)
library(Rcpp)

library(logitnorm)
library(LaplacesDemon)

sourceCpp("spread_functions.cpp")
source("spread_functions.R") # for rast_from_mat

# load Cholila fire data. The raster is used as a template for plotting.
land <- readRDS(file.path("data", "focal fires data", 
                          "landscapes_ig-known_non-steppe", "2015_50.rds"))
land_raster <- rast(file.path("data", "focal fires data", "raw data from GEE",
                              "fire_data_raw_2015_50.tif"))

# Comments ----------------------------------------------------------------

# Intercepts are the vegetation parameters, parameterized without a reference 
# category, and should not be restricted in sign;
# fwi, northing, windir, and slope, should have positive effects;
# elevation should have negative effect.

# fwi and elevation live in the standardized scale, [~ -5, ~ 5],
# while northing, wind, and slope are in [-1, 1].
# the priors for the former could be wider.

# the full-data model will include a random effect at the fire level, which is 
# centred at zero:

# linpred_f = fixed_effects + error_f
# error_f ~ Normal(mu_f, sigma)
# mu_f = 0 + beta_fwi * fwi_f

# where the global intercept for mu_f is  centred at zero for identifiability.
# In this way, fixed_effects represent the linear predictor for an average fire
# when fwi_f = 0.

# Constants ---------------------------------------------------------------

upper_limit <- 0.5
n_veg_types <- 6
n_terrain <- 4
vegetation_names <- c("shrubland", "subalpine", "wet",
                      "dry_a", "dry_b", "steppe")
terrain_names <- c("northing", "elev", "wind", "slope")

# Logistic curves for standardized variables ------------------------------

sd_z <- 5
r_z <- 0.15

par(mfrow = c(1, 2))
#normal
curve(upper_limit * plogis(abs(rnorm(1, 0, sd_z)) * x),
      main = paste("Normal(sd = ", sd_z, ")", sep = ""),
      from = -5, to = 5, col = rgb(0, 0, 0, 0.1),
      ylab = "burn prob",
      xlab = "[-1, 1] variable",
      ylim = c(0, upper_limit))
for(i in 1:1000) {
  curve(upper_limit * plogis(abs(rnorm(1, 0, sd_z)) * x), add = TRUE,
        col = rgb(0, 0, 0, 0.1))
}

# exponential
curve(upper_limit * plogis(rexp(1, r_z) * x),
      main = paste("exponential (r = ", r_z, ")", sep = ""),
      from = -5, to = 5, col = rgb(0, 0, 0, 0.1),
      ylab = "burn prob",
      xlab = "standardized variable",
      ylim = c(0, upper_limit))
for(i in 1:1000) {
  curve(upper_limit * plogis(rexp(1, r_z) * x), add = TRUE,
        col = rgb(0, 0, 0, 0.1))
}

par(mfrow = c(1, 1))

# densities
curve(dnorm(x, 0, sd_z) * 2, from = 0, to = 20, ylim = c(0, 0.3))
curve(dexp(x, r_z), col = 2, add = TRUE)
# curve(dlnorm(x, log(4), 2), col = 4, add = TRUE)   # hard to mimic with log-normal

# fit log-normal to exponential:
rsamps <- rexp(1e5, 0.25)
log_samps <- log(rsamps)
mm <- lm(log_samps ~ 1)
curve(dexp(x, 0.25), col = 2, from = 0, to = 20)
curve(dlnorm(x, log(coef(mm)), sigma(mm)), col = 4, add = TRUE) # feo feo el fit

# Logistic curves for [-1, 1] variables ---------------------------------

sd_01 <- 10
r_01 <- 0.05

par(mfrow = c(1, 2))

#normal
curve(upper_limit * plogis(abs(rnorm(1, 0, sd_01)) * x),
      main = paste("Normal(sd = ", sd_01, ")", sep = ""),
      from = -1, to = 1, col = rgb(0, 0, 0, 0.1),
      ylab = "burn prob",
      xlab = "[-1, 1] variable",
      ylim = c(0, upper_limit))
for(i in 1:1000) {
  curve(upper_limit * plogis(abs(rnorm(1, 0, sd_01)) * x), add = TRUE,
        col = rgb(0, 0, 0, 0.1))
}

# exponential
curve(upper_limit * plogis(rexp(1, r_01) * x),
      main = paste("exponential (r = ", r_01, ")", sep = ""),
      from = -1, to = 1, col = rgb(0, 0, 0, 0.1),
      ylab = "burn prob",
      xlab = "[-1, 1] variable",
      ylim = c(0, upper_limit))
for(i in 1:1000) {
  curve(upper_limit * plogis(rexp(1, r_01) * x), add = TRUE,
        col = rgb(0, 0, 0, 0.1))
}

par(mfrow = c(1, 1))


# densities
curve(dnorm(x, 0, sd_01) * 2, from = 0, to = 20, 
      ylim = c(0, 0.18))
curve(dexp(x, r_01), col = 2, add = TRUE)


# Intercepts distribution -------------------------------------------------

sd_int <- 10
mu_int <- 0#logit(0.2 / upper_limit)
curve(dlogitnorm(x, mu_int, sd_int), n = 1000)

# Function to simulate from the prior -------------------------------------

prior_sim <- function(mu_int = 0, sd_int = 20, r_01 = 0.05, r_z = 0.15,
                      sigma_sd = 3) {
  
  b_fwi <- rexp(1, r_z)
  sigma <- abs(rnorm(1, 0, sigma_sd))
  error <- rnorm(1, b_fwi * rnorm(1), sigma)
  
  intercepts <- rnorm(n_veg_types, 0, sd_int) + error
  names(intercepts) <- vegetation_names
  
  slopes <- c(
    "aspect" = rexp(1, r_01),                 # positive 
    "wind" = rexp(1, r_01),                   # positive
    "elevation" = (-1) * rexp(1, r_z),        # negative
    "slope" = rexp(1, r_01)                   # positive
  )
    
  return(c(intercepts, slopes))
}

prior_sim()


# Graphical prior check ---------------------------------------------------

# matrix to fill
hollow <- land$vegetation
hollow[land$vegetation == 99] <- -1
hollow[land$vegetation < 99] <- 0

pp <- prior_sim()

set.seed(1)
fire_prior <- simulate_fire_cpp(
  terrain = land$terrain, 
  vegetation = land$vegetation,
  ignition_cells = land$ig_rowcol,
  coef = pp,
  n_veg_types = n_veg_types,
  upper_limit = 0.5
)

hollow[land$vegetation < 99] <- 0
hollow[fire_prior == 1] <- 1
r <- rast_from_mat(hollow, land_raster[[1]])
plot(r, col = c("black", "forestgreen", "orange"), type = "classes",
     levels = c("non-burnable", "unburned", "burned"))
north(); sbar(xpd = T, xy = "bottom", col = "white")





# FROM HERE THIS IS EXPLORATION OF SIMILARITY METRICS ---------------------
# WARNING: it's still not updated to the latest data structure.


# Simulate many fires in Cholila landscape ---------------------------------

many <- 200
# 7 * 200 / 60 # 23 min at maximum
# 7 * 100 / 60

l <- lands[["2015_50"]]
dnames <- dimnames(l)$layers

## fire_sim <- matrix(NA, nrow = ncol(l) * nrow(l), ncol = many)

# fire_sim <- vector(mode = "list", length = many)
# rm(r); rm(fire_prior)
# gc()
# time_start <- Sys.time()
# for(i in 1:many) {
#   print(i)
#
#   pp <- prior_sim()
#
#   fire_i <- simulate_fire_compare(
#     landscape = l[, , 1:7],
#     burnable = l[, , "burnable"],
#     ignition_cells = attr(l, "ig_rowcol"),
#     coef = pp,#prior_sim(),
#     wind_layer = which(dnames == "wind") - 1,
#     elev_layer = which(dnames == "elev") - 1,
#     distances = distances,
#     upper_limit = upper_limit
#   )
#
#   #fire_sim[, i] <- as.numeric(fire_i)
#   fire_sim[[i]] <- list()
#   fire_sim[[i]]$fire <- fire_i
#
#   rm(fire_i) # important to avoid crashing because of memory usage
#   gc()
# }
# time_end <- Sys.time()
# (time_end - time_start) # Time difference of 10.537 mins to simulate 100
#                         # Time difference of 19.16281 mins to simulate 200

# saveRDS(fire_sim, file.path("data", "simulations", "fire_sim.rds"))
# save small subset for tests
# saveRDS(fire_sim[1:10], file.path("data", "simulations", "fire_sim_for_testing.rds"))
fire_sim <- readRDS(file.path("data", "simulations", "fire_sim.rds"))

# Compare discrepancies ---------------------------------------------------

example <- compare_fires(fire_sim[[1]], fire_sim[[2]])

nsim <- many
# create array to fill with discrepancies
disc_arr <- array(NA, dim = c(nsim, nsim, length(example) + 2),
                  dimnames = list(fires = 1:nsim,
                                  fires = 1:nsim,
                                  disc = c(names(example), "mean_size", "dif_size")))

# compute fire sizes
sizes_pix <- sapply(1:nsim, function(i) {
  sum(fire_sim[[i]]$counts_veg)
})
sizes <- sizes_pix / sum(as.numeric(l[, , "burnable"]), na.rm = T) * 100 # size as landscape area percentaje

for(c in 1:(nsim-1)) {  ## this loop also takes long
  print(c)
  for(r in (c+1):nsim) {
    # print(paste("col:", c, "/ row:", r))
    # r = 1; c = 2
    disc_arr[r, c, 1:length(example)] <- compare_fires(fire_sim[[c]], fire_sim[[r]])
    disc_arr[r, c, "mean_size"] <- mean(sizes[c], sizes[r])
    disc_arr[r, c, "dif_size"] <- abs(sizes[c] - sizes[r])
  }
}
saveRDS(disc_arr, file.path("data", "simulations", "fire_sim_disc_arr.rds"))
disc_arr <- readRDS(file.path("data", "simulations", "fire_sim_disc_arr.rds"))



# SEGUIR ACÁ, REHACER PLOTS



# clean array and turn into df
clean_arr <- function(x) na.omit(as.numeric(x))
discrep <- apply(disc_arr, 3, clean_arr) %>% as.data.frame()

# for(i in names(discrep)) {
#   hist(discrep[, i], main = i)
# }


# Many variables _________________________________________

# GGally::ggpairs(
#   data = discrep,
#   mapping = aes(alpha = 0.1)
# ) +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank())
#
# ggsave(file.path("files", "discrepancies_comparison_01_many.png"),
#        width = 23, height = 20, units = "cm")


# Overlaps pure _________________________________________

in_vars <- c("overlap_sp", "overlap_vd",
             "overlap_norm", "overlap_expquad", "overlap_quad")

p2 <- GGally::ggpairs(
  data = discrep[, in_vars],
  mapping = aes(alpha = 0.1)
) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(file.path("files", "discrepancies_comparison_02_pure_overlaps.png"),
       plot = p2,
       width = 23, height = 20, units = "cm")


# Overlap combinations _____________________________________

in_vars <- c("overlap_sp",
             "sp_norm_5050",
             "sp_norm_7525",
             "sp_expquad_5050",
             "sp_expquad_7525",
             "sp_quad_5050",
             "sp_quad_7525")

p3 <- GGally::ggpairs(
  data = discrep[, in_vars],
  mapping = aes(alpha = 0.05)
) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(file.path("files", "discrepancies_comparison_03_combinations.png"),
       plot = p3,
       width = 23, height = 20, units = "cm")

# Está bien probar estas combinaciones (y el overlap_sp).

# TAREA: LIMPIAR ESTE CÓDIGO, QUE SEA SOLO EXPLORAR DISCREPANCIAS. DEJAR EL
# PRIOR CHECK EN OTRO LADO.
