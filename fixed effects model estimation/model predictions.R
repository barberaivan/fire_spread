# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(tidyterra)
library(ggspatial)
library(viridis)
library(posterior)
library(coda)
library(FireSpread)
library(DHARMa)

library(Rcpp)
library(terra)
library(abind)         # manage arrays

library(foreach)       # parallelization
library(doMC)          # doRNG had problems with cpp functions

source(file.path("..", "FireSpread", "tests", "testthat", "R_spread_functions.R"))
# for rast_from_mat and a few constants


# Multicore settings -----------------------------------------------------

n_cores <- 15
registerDoMC(n_cores)

# Data to simulate fire ---------------------------------------------------

data_dir <- file.path("data", "focal fires data", "landscapes_ig-known")
filenames <- list.files(data_dir)

fire_ids <- sapply(filenames, function(x) {
  # x <- "fire_data_raw_CoVentana.tif"
  id <- strsplit(x, "[.]")[[1]][1]
  return(id)
}) %>% unname
n_fires <- length(fire_ids)

# dir to save output
target_dir <- file.path("files", "fixed efffects model predictions")

# constants for fire spread simulation
upper_limit <- 0.5
par_names <- c("intercept", "vfi", "tfi", "slope", "wind", "fwi")
n_coef <- length(par_names)

# fwi constants to standardize
fwi_values <- readRDS(file.path("data", "focal fires data", "fwi_values_ig-known.rds"))
fwi_mean <- mean(fwi_values)
fwi_sd <- sd(fwi_values)

# fire sizes
fire_size <- sapply(filenames, function(x) {
  f <- readRDS(file.path(data_dir, x))
  return(ncol(f$burned_ids))
})
names(fire_size) <- fire_ids

# Data --------------------------------------------------------------------

# Load samples (MCMC)

s <- readRDS(file.path("fixed effects model estimation", "draws_01_8chains.rds"))
# extract draws removing adaptation
draws_list <- lapply(1:length(s), function(i) {
  # remove adaptation
  d0 <- s[[i]]$samples[-(1:2000), ] # for 01
  # get a sample every 50 iter
  ii <- seq(1, nrow(d0), by = 50)
  return(mcmc(d0[ii, ]))
})

draws_mcmc <- mcmc.list(draws_list) # from coda
dl <- as_draws_list(draws_mcmc)     # from posterior
sm2 <- posterior::summarise_draws(dl)
print(sm2) # terrible, not converging with niter = 100000, adapt = 2000

# make matrix
smat <- do.call("rbind", draws_list)
d <- ncol(smat)
# transform parameters to original scale
smat[, 2:d] <- exp(smat[, 2:d])

# load samples from rejection sampling
smat <- readRDS(file.path("fixed effects model estimation", "vainilla_abc_samples.rds"))


# Functions ---------------------------------------------------------------

summ <- function(x) {
  qq <- quantile(x, prob = c(0.025, 0.25, 0.75, 0.975))
  names(qq) <- c("lower95", "lower50", "upper50", "upper95")
  return(c("mean" = mean(x), "median" = median(x), qq))
}

# function to simulate a fire and get the overlap and size
simulate_metrics <- function(particle, n_sim = 1, fire_data = NULL) {

  ## testo
  # particle <- particles_sim(N = 1)
  ## end testo
  metrics <- matrix(NA, n_sim, 2)
  colnames(metrics) <- c("overlap", "size")

  # simulate and compute metrics
  for(i in 1:n_sim) { # simulated fires by particle
    fire_sim <- simulate_fire_compare(
      landscape = fire_data$landscape,
      burnable = fire_data$burnable,
      ignition_cells = fire_data$ig_rowcol,
      coef = particle,
      upper_limit = upper_limit
    )

    metrics[i, "overlap"] <- overlap_spatial(
      fire_sim, fire_data[c("burned_layer", "burned_ids")]
    )

    metrics[i, "size"] <- ncol(fire_sim$burned_ids)
  }

  return(metrics)
}


# Emulate the loglik over a list of particles in parallel. Returns the metrics
# in rows.
simulate_metrics_parallel <- function(particles, # matrix with parameter vectors (rows)
                                      fire_data = NULL,
                                      n_sim = 1) {

  # define the fire's intercept as intercept + b_fwi * ((fwi[f] - mean) / sd)
  fwi_z <- (fire_data[["fwi"]]$fwi_expquad_day - fwi_mean) / fwi_sd

  new_intercept <- particles[, "intercept"] + particles[, "fwi"] * fwi_z

  particles2 <- particles[, -6]
  particles2[, "intercept"] <- new_intercept

  particles_list <- lapply(1:nrow(particles2),
                           function(x) particles2[x, ])

  # simulate in parallel
  result <- foreach(pp = particles_list) %dopar% {
    simulate_metrics(pp, fire_data = fire_data, n_sim = n_sim)
  }

  # inherit previous dimnames
  names_single <- colnames(result[[1]])

  # rbind list result
  res <- do.call("rbind", result) %>% as.data.frame()
  names(res) <- names_single

  # write intermediate output to disk
  saveRDS(res, file.path(target_dir,
                         paste("simulations_size-ov_", fire_data$fire_id, ".rds",
                               sep = "")))

  return(res)
}

# function to repeat the simulate_metrics_parallel, returning an array
# (3D: particles, metrics, fires) and writing the result to disk.
# this loads and removes every fire, and prints "simulating tal fire"
simulation_wave <- function(particles, n_sim = 1) {

  sims_list <- lapply(1:n_fires, function(f) {
    message(paste("simulating fire", fire_ids[f]))
    full_data <- readRDS(file.path(data_dir, filenames[f]))
    # subset data needed for spread (to be cloned across workers)
    spread_data <- full_data[c("landscape", "burnable", "ig_rowcol",
                               "burned_layer", "burned_ids",
                               "fwi", "fire_id")]

    res <- simulate_metrics_parallel(particles, spread_data, n_sim = n_sim)
    rm(full_data, spread_data)
    return(res)
  })

  # inherit previous dimnames
  dn_single <- dimnames(sims_list[[1]])

  # turn into array
  res <- abind::abind(sims_list, along = length(dn_single) + 1)

  # set dimnames
  dimnames(res) = c(
    dn_single,
    list(fire_id = fire_ids)
  )

  # write raw output
  saveRDS(res, file.path(target_dir,
                         paste("simulations_size-ov", "_all_fires.rds", sep = "")))

  # remove temporary results (by fire) from disk
  ff_all <- list.files(target_dir)
  ff_delete1 <- ff_all[grep("all_fires", ff_all, invert = TRUE)]
  ff_delete2 <- ff_delete1[grep("simulations_size-ov", ff_delete1)]
  lapply(ff_delete2, function(f) {
    unlink(file.path(target_dir, f))
  })

  gc()

  return(res)
}

# Partial predictions curves ----------------------------------------------

# all predictors can be centred at zero
ns <- 200
seq_list <- list(
  vfi = seq(-5, 5, length.out = ns),
  tfi = seq(-5, 5, length.out = ns),
  slope = seq(-0.8, 0.8, length.out = ns),
  wind = seq(-1, 1, length.out = ns),
  fwi = (seq(min(fwi_values), max(fwi_values), length.out = ns) - fwi_mean) / fwi_sd
)

ones <- rep(1, ns)

var_names <- names(seq_list)
nvar <- length(var_names)

preds <- do.call("rbind", lapply(1:nvar, function(v) {
  # v = 1
  b <- smat[, c(1, v+1)] %>% t
  X <- cbind(ones, seq_list[[v]])
  pp <- plogis(X %*% b) * 0.5
  psumm <- as.data.frame(apply(pp, 1, summ) %>% t)

  dd <- data.frame(
    x = seq_list[[v]],
    x_name = var_names[v]
  )

  d2 <- cbind(dd, psumm)
  return(d2)
}))

# rescale fwi to original scale
preds$x[preds$x_name == "fwi"] <- preds$x[preds$x_name == "fwi"] * fwi_sd + fwi_mean

var_labels <- c(
  "Inflamabilidad por\nvegetación",
  "Inflamabilidad por\ntopografía",
  "Pendiente",
  "Dirección del viento",
  "Fire Weather Index"
)
preds$x_name <- factor(preds$x_name, levels = c("vfi", "tfi", "slope", "wind", "fwi"),
                       labels = var_labels)

ggplot(preds, aes(x = x, y = median, ymin = lower50, ymax = upper50)) +
  geom_ribbon(data = preds, mapping = aes(x = x, ymin = lower95, ymax = upper95),
              color = NA, alpha = 0.3, inherit.aes = F) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  facet_wrap(vars(x_name), scales = "free_x", strip.position = "bottom") +
  ylab("Probabilidad de propagación") +
  theme(axis.title.x = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing = unit(4, "mm"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 11, vjust = 1),
        axis.title.y = element_text(size = 11))

# escalar fwi a los valores originales
ggsave(file.path("files", "images_rae_2023", "fixed_model_partial-predictions.png"),
       width = 17, height = 12, units = "cm")



# The same but plotting curves
preds2 <- do.call("rbind", lapply(1:nvar, function(v) {
  # v = 1
  b <- smat[, c(1, v+1)] %>% t
  X <- cbind(ones, seq_list[[v]])
  pp <- plogis(X %*% b) * 0.5
  psumm <- as.data.frame(apply(pp, 1, summ) %>% t)

  dd <- data.frame(
    x = seq_list[[v]],
    x_name = var_names[v]
  )

  colnames(pp) <- paste(rep("sim", ncol(pp)), 1:ncol(pp), sep = "_")
  d2 <- cbind(dd, psumm, pp)

  return(d2)
}))


# rescale fwi to original scale
preds2$x[preds2$x_name == "fwi"] <- preds2$x[preds2$x_name == "fwi"] * fwi_sd + fwi_mean

# rename predictors
preds2$x_name <- factor(preds2$x_name, levels = c("vfi", "tfi", "slope", "wind", "fwi"),
                       labels = var_labels)

# longanize simulations
sim_names <- grep("sim_", names(preds2))
sim_names <- sample(sim_names, 300)

preds_long <- pivot_longer(preds2, all_of(sim_names), names_to = "iter",
                           values_to = "p")

head(preds_long)

ggplot(preds, aes(x = x, y = median, ymin = lower95, ymax = upper95)) +
  geom_line(data = preds_long, mapping = aes(x = x, y = p, group = iter),
            alpha = 0.05) +
  geom_line(color = viridis(1, begin = 0.1), linewidth = 1) +
  facet_wrap(vars(x_name), scales = "free_x", strip.position = "bottom") +
  ylab("Probabilidad de propagación") +
  theme(axis.title.x = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing = unit(4, "mm"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 11, vjust = 1),
        axis.title.y = element_text(size = 11))

ggsave(file.path("files", "images_rae_2023", "fixed_model_partial-predictions-conditional.png"),
       width = 17, height = 12, units = "cm")

# Prior-posterior ---------------------------------------------------------

d <- ncol(smat)

# prior parameters
mu_int = 0
sd_int = 20
r_slope = 0.04
r_fi = 0.15
r_wind = 0.3
r_fwi = 0.15

# compute density for each parameter and compare with the prior
ns <- 2000
b_seq <- list(
  intercept = seq(-30, 30, length.out = ns),
  vfi = seq(0, 30, length.out = ns),
  tfi = seq(0, 20, length.out = ns),
  slope = seq(0, 150, length.out = ns),
  wind = seq(0, 20, length.out = ns),
  fwi = seq(0, 20, length.out = ns)
)

dprior <- list(
  intercept = dnorm(b_seq$intercept, mu_int, sd_int),
  vfi = dexp(b_seq$vfi, r_fi),
  tfi = dexp(b_seq$tfi, r_fi),
  slope = dexp(b_seq$slope, r_slope),
  wind = dexp(b_seq$wind, r_wind),
  fwi = dexp(b_seq$fwi, r_fwi)
)


dens <- do.call("rbind", lapply(1:d, function(v) {
  dd <- density(smat[, v], n = 2^10,
                from = min(b_seq[[v]]), to = max(b_seq[[v]]))
  dd2 <- data.frame(x = b_seq[[v]],
                    posterior = approx(x = dd$x, y = dd$y, xout = b_seq[[v]])$y,
                    previa = dprior[[v]],
                    ylow = rep(0, ns),
                    par = names(dprior)[v])
  d_long <- pivot_longer(dd2, c(2,3), names_to = "dens_type", values_to = "dens")
  return(d_long)
}))

dens$par <- factor(dens$par, levels = names(dprior),
                   labels = c("Intercept", var_labels))
dens$dens_type <- factor(dens$dens_type, levels = c("previa", "posterior"),
                         labels = c("Previa", "Posterior"))

ggplot(dens, aes(x = x, ymin = ylow, ymax = dens, fill = dens_type)) +
  geom_ribbon(alpha = 0.5) +
  facet_wrap(vars(par), scales = "free", strip.position = "bottom") +
  scale_fill_viridis(discrete = TRUE, end = 0.7, direction = -1) +
  ylab("Densidad") +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing = unit(4, "mm"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 11, vjust = 1),
        axis.title.y = element_text(size = 11))

ggsave(file.path("files", "images_rae_2023", "fixed_model_prior-posterior.png"),
       width = 17, height = 12, units = "cm")


# DHARMa residuals for fire size ------------------------------------------
800 * 20 # 16000 simulations by wave
# use 10000 samples from the posterior, which should take around 4 h to run.

nsim <- 8000
# output: matrix[nfires, nsim] with simulated fire size, to be compared with
# a vector of observed fire sizes.

set.seed(23143)
draws <- smat[sample(1:nrow(smat), size = nsim, replace = F), ]

# simulate all fires
# size_ov_sim <- simulation_wave(draws, n_sim = 1)
# saveRDS(size_ov_sim, file.path(target_dir, "size_ov_simulation_dharma.rds"))
size_ov_sim <- readRDS(file.path(target_dir, "size_ov_simulation_dharma.rds"))

size_sim <- size_ov_sim[, "size", ] %>% t %>% "*"(0.09) # turn into ha
size_obs <- fire_size * 0.09

# DHARMa checks
res_size <- createDHARMa(simulatedResponse = size_sim,
                         observedResponse = size_obs,
                         integerResponse = F,
                         fittedPredictedResponse = apply(size_sim, 1, median))

plotResiduals(res_size, form = log(size_obs), rank = F)
# parece haber sobredispersión? a los fuegos más grandes los simula más chicos,
# y a los fuegos más chicos los simula más grandes.

size_sim_median <- apply(size_sim, 1, median)
size_sim_mean <- apply(size_sim, 1, mean)

par(mfrow = c(2, 2))
plotResiduals(res_size, form = log(size_obs), rank = F)
plotResiduals(res_size, form = log(size_sim_median), rank = F)
plotResiduals(res_size, form = log(size_sim_mean), rank = F)
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
plotResiduals(res_size, form = size_obs, rank = F)
plotResiduals(res_size, form = size_sim_median, rank = F)
plotResiduals(res_size, form = size_sim_mean, rank = F)
par(mfrow = c(1, 1))

# bastante fiero
plot(size_obs ~ size_sim_median)
plot(log(size_obs) ~ log(size_sim_median),
     ylim = c(0, 14), xlim = c(0, 14))
abline(0, 1)

plot(log(size_obs) ~ log(size_sim_mean),
     ylim = c(0, 15), xlim = c(0, 15))
abline(0, 1)

# view densities
data_fit <- data.frame(
  res = res_size$scaledResiduals,
  fire_id = fire_ids,
  size = size_obs
)
data_fit$size_sim <- size_sim

data_fit <- data_fit[order(data_fit$res), ]
row.names(data_fit) <- 1:nrow(data_fit)

f <- 1
plot(density(data_fit$size_sim[f, ], from = min(data_fit$size_sim[f, ]),
             to = max(data_fit$size_sim[f, ])),
     main = data_fit$fire_id[f])
abline(v = data_fit$size[f], lty = 2, col = "red")

f <- 32
plot(density(data_fit$size_sim[f, ], from = min(data_fit$size_sim[f, ]),
             to = max(data_fit$size_sim[f, ])),
     main = data_fit$fire_id[f])
abline(v = data_fit$size[f], lty = 2, col = "red")

f <- 53
plot(density(data_fit$size_sim[f, ], from = min(data_fit$size_sim[f, ]),
             to = max(data_fit$size_sim[f, ])),
     main = data_fit$fire_id[f])
abline(v = data_fit$size[f], lty = 2, col = "red")


# Nice dharma plot
dharma1 <- data.frame(
  res = res_size$scaledResiduals[order(res_size$scaledResiduals)],
  qunif = ppoints(n_fires)
  # size_obs = size_obs,
  # size_sim_median = size_sim_median,
  # size_sim_mean = size_sim_mean
)

dharma2 <- data.frame(
  res = res_size$scaledResiduals,
  size_obs = size_obs,
  size_sim_median = size_sim_median,
  size_sim_mean = size_sim_mean
)

d1 <- ggplot(dharma1, aes(x = qunif, y = res)) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11)) +
  ylab("Cuantiles observados") +
  xlab("Cuantiles esperados")

d2a <- ggplot(dharma2, aes(x = size_sim_median, y = res)) +
  geom_vline(xintercept = range(dharma2$size_obs),
             linetype = "dashed", alpha = 0.5) +
  geom_smooth(method = "gam", method.args = list(family = mgcv::betar())) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11)) +
  ylab("Residuos") +
  xlab("Tamaño predicho (ha)") +
  scale_x_continuous(trans = "log10")

d2b <- ggplot(dharma2, aes(x = size_obs, y = res)) +
  geom_smooth(method = "gam", method.args = list(family = mgcv::betar())) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 11)) +
  ylab(NULL) +
  xlab("Tamaño observado (ha)") +
  scale_x_continuous(trans = "log10")

blank <- ggplot() + theme_void()#ggplot(dharma2, aes(x = size_obs, y = res)) + geom_blank()


dharma_plot <- egg::ggarrange(d1, blank, d2a, d2b, nrow = 2)
ggsave(file.path("files", "images_rae_2023", "fixed_model_dharma.png"), plot = dharma_plot,
       width = 17, height = 12, units = "cm")

# Simulated overlaps ------------------------------------------------------

ov_sim <- size_ov_sim[, "overlap", ] %>% t
ov_summ <- apply(ov_sim, 1, summ) %>% t %>% as.data.frame
ov_summ$fire_id <- rownames(ov_summ)
# hist(ov_summ$upper95, breaks = 10)

ov_summ2 <- left_join(ov_summ, data_fit[, c("fire_id", "res", "size")],
                      by = "fire_id")
ov_summ2 <- ov_summ2[order(ov_summ2$size), ]
ov_summ2$fire_id <- factor(ov_summ2$fire_id, levels = ov_summ2$fire_id)

ggplot(ov_summ2, aes(x = fire_id, y = mean, ymin = lower50, ymax = upper50,
                     color = res)) +
  geom_linerange(mapping = aes(x = fire_id, ymin = lower95, ymax = upper95),
                 alpha = 0.4, linewidth = 1) +
  geom_linerange(linewidth = 1, alpha = 0.5) +
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
        legend.position = "bottom",
        axis.title = element_text(size = 11),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_color_viridis(option = "B", end = 0.8,
                      name = "Proporción de simulaciones\ncon tamaño igual o menor\nal observado",
                      breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  # guides(color = gu)
  xlab("Incendios, de menor a mayor") +
  ylab("Solapamiento\n(observado, simulado)")

ggsave(file.path("files", "images_rae_2023", "fixed_model_fit_overlap_and_res.png"),
       width = 22, height = 14, units = "cm")




# Making maps for 3 fires -------------------------------------------------

# 1999_25j        # balcon: buen ajuste
# 2009_2007583421 # tromen: low dharma, large fire
# 2015_47S        # turbio 2015 S: high dharma, large fire

map_ids <- c("1999_25j", "2009_2007583421", "2015_47S")

raster_files <- list.files(file.path("data", "focal fires data", "raw data from GEE"))
raster_files <- sapply(1:3, function(i) {
  raster_files[grep(map_ids[i], raster_files)]
})

rasters <- lapply(raster_files, function(x) {
  rast(file.path("data", "focal fires data", "raw data from GEE", x))
})

names(rasters) <- c("balcon", "tromen", "turbio")

plot(rasters[[1]])
plot(rasters[[2]])
plot(rasters[[3]])



# load data to simulate fires
data_balcon <- readRDS(file.path(data_dir, "1999_25j.rds"))
data_turbio <- readRDS(file.path(data_dir, "2015_47S.rds"))
data_tromen <- readRDS(file.path(data_dir, "2009_2007583421.rds"))


# Function to simulate fires and record number of times each pixel was burned
count_burns <- function(fire_data, particles) {
  # testing
  # fire_data <- full_data1
  # n_by_core <- 200
  # particles <- smat[1:(1 * n_by_core), ]
  ####
  n_part <- nrow(particles)

  counts <- matrix(0, nrow = nrow(fire_data$burnable),
                   ncol = ncol(fire_data$burnable))

  for(i in 1:n_part) {
    # i = 1
    fsim <- simulate_fire(
      landscape = fire_data$landscape,
      burnable = fire_data$burnable,
      ignition_cells = fire_data$ig_rowcol,
      coef = particles[i, ],
      upper_limit = upper_limit
    )

    counts <- counts + fsim
  }

  attr(counts, "nsim") <- n_part
  return(counts)
}

simulate_burn_prob <- function(full_data, n_by_core = 500) {

  # testing
  # full_data <- full_data1
  # n_by_core <- 30
  #
  n_part <- n_cores * n_by_core
  particles <- smat[1:n_part, ]

  # define the fire's intercept as intercept + b_fwi * ((fwi[f] - mean) / sd)
  fwi_z <- (full_data[["fwi"]]$fwi_expquad_day - fwi_mean) / fwi_sd

  new_intercept <- particles[, "intercept"] + particles[, "fwi"] * fwi_z

  particles2 <- particles[, -6]
  particles2[, "intercept"] <- new_intercept

  # make a list of particles matrices for each core
  core_id <- rep(1:n_cores, each = n_by_core)
  particles_list <- lapply(1:n_cores, function(c) {
    particles2[core_id == c, ]
  })

  # subset data needed for spread (to be cloned across workers)
  spread_data <- full_data[c("landscape", "burnable", "ig_rowcol",
                             "burned_layer", "burned_ids",
                             "fwi", "fire_id")]

  # count burns in parallel
  result <- foreach(pp = particles_list) %dopar% {
    count_burns(pp, fire_data = spread_data)
  }

  res_arr <- simplify2array(result)
  probs <- apply(res_arr, 1:2, sum) / n_part
  return(probs)
}

# compute burn prob maps
prob_balcon <- simulate_burn_prob(data_balcon)
prob_tromen <- simulate_burn_prob(data_tromen)
prob_turbio <- simulate_burn_prob(data_turbio)

rastprobs <- list(
  balcon = rast_from_mat(prob_balcon, rasters[["balcon"]][[1]]),
  tromen = rast_from_mat(prob_tromen, rasters[["tromen"]][[1]]),
  turbio = rast_from_mat(prob_turbio, rasters[["turbio"]][[1]])
)
# saveRDS(rastprobs, file.path(target_dir, "burn_prob_rasters.rds"))

rastfull <- lapply(1:3, function(i) {
  r <- c(rasters[[i]], rastprobs[[i]])
  names(r)[nlyr(r)] <- "prob"
  return(r)
})
names(rastfull) <- names(rasters)


plot(rastfull[["balcon"]][[5:7]])
plot(rastfull[["tromen"]][[5:7]])
plot(rastfull[["turbio"]][[5:7]])

# plotear mejor lueguito.



# Maps with tidyterra -----------------------------------------------------

polys <- vect("data/patagonian_fires_spread.shp")
polys <- project(polys, "EPSG:5343")

points <- vect("data/ignition_points_checked.shp")

poli3 <- list(
  balcon = polys[polys$fire_id == "1999_25j"],
  tromen = polys[polys$fire_id == "2009_2007583421"],
  turbio = polys[polys$fire_id == "2015_47S"]
)

points3 <- list(
  balcon = points[points$Name == "1999_25j"],
  tromen = points[points$Name == "2009_2007583421"],
  turbio = points[points$Name == "2015_47S"]
)


## setting balcon non-burnable layer

# make layer with high_andean and lakes
elev_values <- values(rastfull[["balcon"]][["elev"]])
nb_rast <- rastfull[["balcon"]][["burnable"]]
nb_values <- values(nb_rast)
# make burnable the northern portion of Gutierrez (villa los coihues)
rc <- expand.grid(
  rows = 1:floor(nrow(nb_rast) * 0.14),
  cols = floor(ncol(nb_rast) * 0.5) : ncol(nb_rast)
)
ids_town <- cellFromRowCol(nb_rast, rc$rows, rc$cols)

nb_values[ids_town] <- 1
nb_balcon <- nb_rast
values(nb_balcon) <- nb_values


high_andean <- which(nb_values == 0 & elev_values > 1100)
lake <- which(nb_values == 0 & elev_values < 1000)
burnable <- which(nb_values == 1)

values(nb_balcon)[burnable] <- NA
values(nb_balcon)[high_andean] <- 1
values(nb_balcon)[lake] <- 2

bp_balcon <- rastfull[["balcon"]][["prob"]]
values(bp_balcon)[!is.na(values(nb_balcon))] <- NA

# make factor the category
nb_balcon <- as.factor(nb_balcon)
levels(nb_balcon) <- data.frame(
  value = 1:2, desc = c("Altoandino", "Lagos")
)

# colcol <- c("#EEE8CD", "#53868B") # bueno, clarito
# colcol <- c("#8B8B83", "#53868B")
colcol <- c("gray30", "#4A708B")  # bueno, oscurito

(plot_balcon <- ggplot(poli3$balcon) +

  # non-burnable
  scale_fill_manual(values = colcol,
                    na.translate = F,
                    na.value = "transparent",
                    name = "No quemable",
                    guide = guide_legend(order = 2),
                    drop = TRUE) +
  scale_alpha_manual(na.value = 0) +
  geom_spatraster(data = nb_balcon) +

  # burn probability
  ggnewscale::new_scale_fill() +
  scale_fill_viridis(na.value = "transparent", option = "F",
                     direction = 1, begin = 0, end = 0.9,
                     guide = guide_colourbar(order = 1),
                     name = "Probabilidad de quemarse") +
  geom_spatraster(data = bp_balcon) +

  # fire
  geom_spatvector(fill = NA, size = 1, color = "white", linewidth = 0.4) +

  # ignition point
  geom_spatvector(data = points3$balcon,
                  fill = "white", size = 2, color = "black",
                  shape = 21, stroke = 1) +

  # map scale
  annotation_scale(height = unit(1, "mm"),
                   text_col = "white", location = "tr", style = "ticks",
                   line_col = "white") +

  theme(legend.title = element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
        legend.position = "right") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Balcón del Gutiérrez (1999)"))


# non burnable for turbio

# make layer with high_andean and lakes
elev_values <- values(rastfull[["turbio"]][["elev"]])
nb_turbio <- rastfull[["turbio"]][["burnable"]]
nb_values <- values(nb_turbio)

high_andean <- which(nb_values == 0)# & elev_values >= 600)
lake <- which(nb_values == 0 & elev_values < 600)
burnable <- which(nb_values == 1)

values(nb_turbio)[burnable] <- NA
values(nb_turbio)[high_andean] <- 1
values(nb_turbio)[lake] <- 2

bp_turbio <- rastfull[["turbio"]][["prob"]]
values(bp_turbio)[!is.na(values(nb_turbio))] <- NA

# turn fake lakes into high-andean to the west
rc <- expand.grid(row = 1:nrow(nb_turbio), col = 1:(ncol(nb_turbio) * 0.25))
cc <- cellFromRowCol(nb_turbio, row = rc$row, col = rc$col)
target_values <- values(nb_turbio)[cc]
change_these <- cc[which(target_values == 2)]
values(nb_turbio)[change_these] <- 1

# make factor the category
nb_turbio <- as.factor(nb_turbio)
levels(nb_turbio) <- data.frame(
  value = 1:2, desc = c("Altoandino", "Lagos")
)
# plot(nb_turbio)
## plot turbio
(plot_turbio <- ggplot(poli3$turbio) +

    # non-burnable
    scale_fill_manual(values = colcol,
                      na.translate = F,
                      na.value = "transparent",
                      name = "No quemable",
                      drop = TRUE,
                      guide = guide_legend(order = 2)) +
    geom_spatraster(data = nb_turbio) +

    # burn probability
    ggnewscale::new_scale_fill() +
    scale_fill_viridis(na.value = "transparent", option = "F",
                       limits = c(0, 1), discrete = F,
                       direction = 1, begin = 0, end = 0.9,
                       guide = guide_colourbar(order = 1),
                       name = "Probabilidad de quemarse") +
    geom_spatraster(data = bp_turbio) +

    # fire
    geom_spatvector(fill = NA, size = 1, color = "white", linewidth = 0.4) +

    # ignition point
    geom_spatvector(data = points3$turbio,
                    fill = "white", size = 2, color = "black",
                    shape = 21, stroke = 1) +

    # map scale
    annotation_scale(height = unit(1, "mm"),
                     text_col = "white", location = "tr", style = "ticks",
                     line_col = "white") +

    theme(legend.title = element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
          legend.position = "right") +
    # guides(fill = guide(order = 1)) +

    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +

    ggtitle("Río Turbio (2015)"))



# non burnable for tromen

# make layer with high_andean and lakes
elev_values <- values(rastfull[["tromen"]][["elev"]])
nb_tromen <- rastfull[["tromen"]][["burnable"]]
nb_values <- values(nb_tromen)

high_andean <- which(nb_values == 0)# & elev_values >= 600)
lake <- which(nb_values == 0 & elev_values < 1100)
burnable <- which(nb_values == 1)

values(nb_tromen)[burnable] <- NA
values(nb_tromen)[high_andean] <- 1
values(nb_tromen)[lake] <- 2

# correct some lakes in the north
rc <- expand.grid(row = 1:(nrow(nb_tromen) * 0.28), col = 1:(ncol(nb_tromen) * 0.5))
cc <- cellFromRowCol(nb_tromen, row = rc$row, col = rc$col)
target_values <- values(nb_tromen)[cc]
change_these <- cc[which(target_values == 1)]
values(nb_tromen)[change_these] <- 2

# and now, unlake the portion of high-andean
rc <- expand.grid(row = 1:(nrow(nb_tromen) * 0.1), col = 1:(ncol(nb_tromen) * 0.1))
cc <- cellFromRowCol(nb_tromen, row = rc$row, col = rc$col)
target_values <- values(nb_tromen)[cc]
change_these <- cc[which(target_values == 2)]
values(nb_tromen)[change_these] <- 1

# make factor the category
nb_tromen <- as.factor(nb_tromen)
levels(nb_tromen) <- data.frame(
  value = 1:2, desc = c("Altoandino", "Lagos")
)
# plot(nb_tromen)

bp_tromen <- rastfull[["tromen"]][["prob"]]
values(bp_tromen)[!is.na(values(nb_tromen))] <- NA

## plot tromen
(plot_tromen <- ggplot(poli3$tromen) +

    # non-burnable
    scale_fill_manual(values = colcol,
                      na.translate = F,
                      na.value = "transparent",
                      name = "No quemable",
                      guide = guide_legend(order = 2),
                      drop = TRUE) +
    scale_alpha_manual(na.value = 0) +
    geom_spatraster(data = nb_tromen) +

    # burn probability
    ggnewscale::new_scale_fill() +
    scale_fill_viridis(na.value = "transparent", option = "F",
                       direction = 1, begin = 0, end = 0.9,
                       guide = guide_colourbar(order = 1),
                       name = "Probabilidad de quemarse") +
    geom_spatraster(data = bp_tromen) +

    # fire
    geom_spatvector(fill = NA, size = 1, color = "white", linewidth = 0.4) +

    # ignition point
    geom_spatvector(data = points3$tromen,
                    fill = "white", size = 2, color = "black",
                    shape = 21, stroke = 1) +

    # map scale
    annotation_scale(height = unit(1, "mm"),
                     text_col = "white", location = "tr", style = "ticks",
                     line_col = "white") +

    theme(legend.title = element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
          legend.position = "right") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("Tromen (2009)"))



# save maps

ggsave(file.path("files", "images_rae_2023", "fixed_model_maps_balcon.png"),
       plot = plot_balcon, width = 18, height = 12, units = "cm")

ggsave(file.path("files", "images_rae_2023", "fixed_model_maps_turbio.png"),
       plot = plot_turbio, width = 18, height = 12, units = "cm")

ggsave(file.path("files", "images_rae_2023", "fixed_model_maps_tromen.png"),
       plot = plot_tromen, width = 18, height = 12, units = "cm")



# Map with predictors (Balcon) --------------------------------------------

## correcting balcon non-burnable layer

# make layer with high_andean and lakes
elev_values <- values(rasters[["balcon"]][["elev"]])
veg_rast <- rasters[["balcon"]][["veg"]]
plot(veg_rast)
veg_values <- values(veg_rast)

# make burnable (shrubland) the northern portion of Gutierrez (villa los coihues)
rc <- expand.grid(
  rows = 1:floor(nrow(veg_rast) * 0.14),
  cols = floor(ncol(veg_rast) * 0.5) : ncol(veg_rast)
)
ids_town0 <- cellFromRowCol(veg_rast, rc$rows, rc$cols)
ids_town <- ids_town0[which(veg_values[ids_town0] == 1)]
veg_values[ids_town] <- 5
values(veg_rast) <- veg_values
# plot(veg_rast)

# separate non burnable into lake and high-andean
high_andean <- which(veg_values == 1 & elev_values > 1100)
lake <- which(veg_values == 1 & elev_values < 1100)

values(veg_rast)[high_andean] <- 1
values(veg_rast)[lake] <- 0

# make factor the category
veg_rast <- as.factor(veg_rast)
levels(veg_rast) <- data.frame(
  value = 0:5,
  desc = c("Lago", "Altoandino", "Bosque de lenga", "Bosque de coihue",
           "Bosque de ciprés", "Matorral")
)

plot(veg_rast)

### Plots de predictoras

# rangos de coordenadas
x_range <- c(1542390, 1552710)
y_range <- c(5433660, 5444640)
x_seq <- seq(x_range[1], x_range[2], length.out = 4)
y_seq <- seq(y_range[1], y_range[2], length.out = 4)


barplot(rep(1, 4), col = viridis(4, begin = 0.5, end = 1))
k <- 6
barplot(rep(1, k), col = hcl.colors(k, palette = "Greens"))
# vegetation colors
vegcolor <- c(#"#4A708B", "gray30", # usados en mapas de prob
              hcl.colors(1, palette = "Blues 3"), # Lakes
              hcl.colors(1, palette = "Grays"), # Altoandino

              # viridis(4, begin = 0.5, end = 1))
              hcl.colors(k, palette = "Greens")[1:4])

# vegetación
(map_veg <- ggplot() +
  geom_spatraster(data = veg_rast) +
  # scale_fill_viridis(discrete = TRUE, name = "Tipo de vegetación",
  #                    option = "B") +

  # green colors
  scale_fill_manual(values = vegcolor) +

  # map scale
  # annotation_scale(height = unit(1, "mm"),
  #                  text_col = "black", location = "tr", style = "ticks",
  #                  line_col = "black",
  #                  text_cex = 0.6) +

  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.text = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 11),
        legend.text = element_text(size = 8)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Vegetación"))

# altitud
(map_alt <- ggplot() +
    geom_spatraster(data = rasters[["balcon"]][["elev"]]) +
    scale_fill_viridis(discrete = F, name = "Altitud (msnm)",
                       option = "D") +
    # map scale
    # annotation_scale(height = unit(1, "mm"),
    #                  text_col = "white", location = "tr", style = "ticks",
    #                  line_col = "white",
    #                  text_cex = 0.6) +

    theme(legend.title = element_blank(),#element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
          legend.position = "right",
          axis.text = element_text(size = 6),
          axis.text.x = element_text(angle = 50, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 11),
          legend.text = element_text(size = 8)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("Altitud (msnm)"))

# fire
# make factor the category
fire_rast <- veg_rast
values(fire_rast)[values(fire_rast) > 1] <- 2
values(fire_rast)[values(fire_rast) <= 1] <- 1
values(fire_rast)[values(rasters[["balcon"]][["burned"]]) == 1] <- 3

levels(fire_rast) <- data.frame(
  value = 1:3,
  desc = c("No quemable", "Quemable", "Quemado")
)

(map_fire <- ggplot() +
    geom_spatraster(data = fire_rast) +
    scale_fill_manual(
      values = c(
        viridis(1, option = "A"),
        viridis(1, option = "A", begin = 1),
        viridis(1, option = "A", begin = 0.6)
      )) +

    # scale_fill_viridis(option = "A", discrete = TRUE, begin = 0.6,
    #                    na.value = "transparent",
    #                    na.translate = F) +


    # map scale
    annotation_scale(height = unit(1, "mm"),
                     text_col = "black", location = "tr", style = "ticks",
                     line_col = "black",
                     text_cex = 0.6) +

    # ignition point
    geom_spatvector(data = points3$balcon,
                    fill = "white", size = 1.2, color = "black",
                    shape = 21, stroke = 0.6) +

    theme(legend.title = element_blank(),#element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
          legend.position = "right",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          # axis.text = element_text(size = 6),
          # axis.text.x = element_text(angle = 50, hjust = 1),
          plot.title = element_text(size = 11),
          legend.text = element_text(size = 8)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("Incendio"))


# merge maps

map_predictors <- egg::ggarrange(map_veg, map_alt, map_fire,
                                 nrow = 2)

ggsave(file.path("files", "images_rae_2023", "fixed_model_maps_predictors.png"),
       plot = map_predictors, width = 17, height = 13, units = "cm")


# Map overlap examples (balcon) -------------------------------------------

fire_obs <- fire_rast
fire_sim <- fire_obs

# make burned burnable
values(fire_sim)[values(fire_obs) == 3] <- 2

# simulate fire 1
set.seed(23143)
bb <- smat[sample(1:nrow(smat), size = 1), ]
fwi_z <- (data_balcon[["fwi"]]$fwi_expquad_day - fwi_mean) / fwi_sd
new_intercept <- bb["intercept"] + bb["fwi"] * fwi_z
bb2 <- bb[-6]
bb2["intercept"] <- new_intercept

balcon_sim_raw <- simulate_fire_compare(
  landscape = data_balcon$landscape,
  burnable = data_balcon$burnable,
  ignition_cells = data_balcon$ig_rowcol,
  coef = bb2,
  upper_limit = upper_limit
)

simrast0 <- rast_from_mat(balcon_sim_raw$burned_layer, fire_sim)
values(fire_sim)[values(simrast0) == 1] <- 3
plot(fire_sim)

ov1 <- overlap_spatial(
  balcon_sim_raw, data_balcon[c("burned_layer", "burned_ids")]
)

fire_sim <- as.factor(fire_sim)
levels(fire_sim) <- data.frame(
  value = 1:3,
  desc = c("No quemable", "Quemable", "Quemado")
)


### hacemos más capas para overlapear.

# every category will be a layer, so all legends are equivalent

# burnable
# burnable_layer <- fire_sim
# values(burnable_layer)[values(burnable_layer) != 2] <- NA
# burnable_layer <- as.factor(burnable_layer)
# levels(burnable_layer) <- data.frame(
#   value = 2,
#   desc = "Quemable"
# )
burnable_layer <- fire_sim
values(burnable_layer)[values(burnable_layer) > 1] <- 2
values(burnable_layer)[values(burnable_layer) == 1] <- NA
burnable_layer <- as.factor(burnable_layer)
levels(burnable_layer) <- data.frame(
  value = 2,
  desc = "Quemable"
)


# no quemable
nonburnable_layer <- fire_sim
values(nonburnable_layer)[values(nonburnable_layer) != 1] <- NA
nonburnable_layer <- as.factor(nonburnable_layer)
levels(nonburnable_layer) <- data.frame(
  value = 1,
  desc = "No quemable"
)
# plot(nonburnable_layer)

# incendio simulado solo (NA en lo otro)
fire_sim_alone <- fire_sim
values(fire_sim_alone)[values(fire_sim_alone) < 3] <- NA
# plot(fire_sim_alone)
levels(fire_sim_alone) <- data.frame(
  value = 3,
  desc = c("Incendio simulado")
)

# observado solo
fire_obs_alone <- fire_obs
values(fire_obs_alone)[values(fire_obs_alone) < 3] <- NA
# plot(fire_sim_alone)
levels(fire_obs_alone) <- data.frame(
  value = 3,
  desc = c("Incendio real")
)
# plot(fire_obs_alone)

(map_fire_sim1 <- ggplot() +
    geom_spatraster(data = nonburnable_layer) +
    scale_fill_manual(values = "black",
                      na.value = "transparent", na.translate = F,
                      guide = guide_legend(order = 1)) +

    ggnewscale::new_scale_fill() +
    geom_spatraster(data = burnable_layer) +
    scale_fill_viridis(option = "A", discrete = TRUE, begin = 1,
                       na.value = "transparent",
                       na.translate = F,
                       guide = guide_legend(order = 2)) +

    ggnewscale::new_scale_fill() +
    geom_spatraster(data = fire_sim_alone, alpha = 1) +
    scale_fill_viridis(option = "A", discrete = TRUE, begin = 0.2,
                       na.value = "transparent",
                       na.translate = F,
                       guide = guide_legend(order = 3)) +

    # ignition point
    geom_spatvector(data = points3$balcon,
                    fill = "white", size = 1.2, color = "black",
                    shape = 21, stroke = 0.8) +

    theme(legend.title = element_blank(),#element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
          legend.position = "right",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 11),
          legend.text = element_text(size = 8)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("Incendio simulado"))


(map_fire_union <- ggplot() +
  geom_spatraster(data = nonburnable_layer) +
  scale_fill_manual(values = "black",
                    na.value = "transparent", na.translate = F,
                    guide = guide_legend(order = 1)) +

  ggnewscale::new_scale_fill() +
  geom_spatraster(data = burnable_layer) +
  scale_fill_viridis(option = "A", discrete = TRUE, begin = 1,
                     na.value = "transparent",
                     na.translate = F,
                     guide = guide_legend(order = 2)) +

  ggnewscale::new_scale_fill() +
  geom_spatraster(data = fire_sim_alone, alpha = 1) +
  scale_fill_viridis(option = "A", discrete = TRUE, begin = 0.2,
                     na.value = "transparent",
                     na.translate = F,
                     guide = guide_legend(order = 3)) +

  ggnewscale::new_scale_fill() +
  geom_spatraster(data = fire_obs_alone, alpha = 0.7) +
  scale_fill_viridis(option = "A", discrete = TRUE, begin = 0.6,
                     na.value = "transparent",
                     na.translate = F,
                     guide = guide_legend(order = 4)) +

  # ignition point
  geom_spatvector(data = points3$balcon,
                  fill = "white", size = 1.2, color = "black",
                  shape = 21, stroke = 0.8) +

    theme(legend.title = element_blank(),#element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
          legend.position = "right",
          axis.text.y = element_text(size = 6),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(size = 11),
          legend.text = element_text(size = 8)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("Superposición"))


# transparent observed fire
(map_fire_transparent <- ggplot() +
    geom_spatraster(data = nonburnable_layer) +
    scale_fill_manual(values = "black",
                      na.value = "transparent", na.translate = F,
                      guide = guide_legend(order = 1)) +

    ggnewscale::new_scale_fill() +
    geom_spatraster(data = burnable_layer) +
    scale_fill_viridis(option = "A", discrete = TRUE, begin = 1,
                       na.value = "transparent",
                       na.translate = F,
                       guide = guide_legend(order = 2)) +

    ggnewscale::new_scale_fill() +
    geom_spatraster(data = fire_obs_alone, alpha = 0.7) +
    scale_fill_viridis(option = "A", discrete = TRUE, begin = 0.6,
                       na.value = "transparent",
                       na.translate = F,
                       guide = guide_legend(order = 4)) +

    # ignition point
    geom_spatvector(data = points3$balcon,
                    fill = "white", size = 1.2, color = "black",
                    shape = 21, stroke = 0.8) +

    # map scale
    annotation_scale(height = unit(1, "mm"),
                     text_col = "black", location = "tr", style = "ticks",
                     line_col = "black",
                     text_cex = 0.6) +

    theme(legend.title = element_blank(),
          legend.position = "right",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 11),
          legend.text = element_text(size = 8)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("Incendio real"))


# Intersection
# Make NA the burned area if the other fire is NA
fire_obs_int <- fire_obs_alone
values(fire_obs_int)[is.na(values(fire_sim_alone))] <- NA

fire_sim_int <- fire_sim_alone
values(fire_sim_int)[is.na(values(fire_obs_alone))] <- NA

levels(fire_obs_int) <- data.frame(
  value = 3,
  desc = c("Incendio real")
)
levels(fire_sim_int) <- data.frame(
  value = 3,
  desc = c("Incendio simulado")
)


(map_fire_int <- ggplot() +
    geom_spatraster(data = nonburnable_layer) +
    scale_fill_manual(values = "black",
                      na.value = "transparent", na.translate = F,
                      guide = guide_legend(order = 1)) +

    ggnewscale::new_scale_fill() +
    geom_spatraster(data = burnable_layer) +
    scale_fill_viridis(option = "A", discrete = TRUE, begin = 1,
                       na.value = "transparent",
                       na.translate = F,
                       guide = guide_legend(order = 2)) +

    ggnewscale::new_scale_fill() +
    geom_spatraster(data = fire_sim_int, alpha = 1) +
    scale_fill_viridis(option = "A", discrete = TRUE, begin = 0.2,
                       na.value = "transparent",
                       na.translate = F,
                       guide = guide_legend(order = 3)) +

    ggnewscale::new_scale_fill() +
    geom_spatraster(data = fire_obs_int, alpha = 0.7) +
    scale_fill_viridis(option = "A", discrete = TRUE, begin = 0.6,
                       na.value = "transparent",
                       na.translate = F,
                       guide = guide_legend(order = 4)) +

    # ignition point
    geom_spatvector(data = points3$balcon,
                    fill = "white", size = 1.2, color = "black",
                    shape = 21, stroke = 0.8) +

    theme(legend.title = element_blank(),#element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
          legend.position = "right",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 50, hjust = 1, size = 6),
          plot.title = element_text(size = 11),
          legend.text = element_text(size = 8),
          legend.spacing.y = unit(0.1, 'mm'),
          legend.margin = margin(),)+
          # legend.spacing = unit(0.1, "mm")) +

    guides(fill = guide_legend(byrow = TRUE)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("Intersección"))

maps_fires <- egg::ggarrange(map_fire_transparent + theme(legend.position = "none"),
                             map_fire_sim1 + theme(legend.position = "none"),
                             map_fire_union + theme(legend.position = "none"),
                             map_fire_int,
                             nrow = 2)

ggsave(file.path("files", "images_rae_2023", "fixed_model_maps_overlap.png"),
       plot = maps_fires, width = 14, height = 12, units = "cm")
ov1 # 0.5153793


# Spread animation by hand ------------------------------------------------

set.seed(23143)
bb <- smat[sample(1:nrow(smat), size = 1), ]
fwi_z <- (data_balcon[["fwi"]]$fwi_expquad_day - fwi_mean) / fwi_sd
new_intercept <- bb["intercept"] + bb["fwi"] * fwi_z
bb2 <- bb[-6]
bb2["intercept"] <- new_intercept

balcon_sim_anim <- simulate_fire_animate(
  landscape = data_balcon$landscape,
  burnable = data_balcon$burnable,
  ignition_cells = data_balcon$ig_rowcol,
  coef = bb2,
  upper_limit = upper_limit
)

simrast_anim <- rast_from_mat(balcon_sim_anim, fire_sim)
values(simrast_anim)[values(simrast_anim) < 1] <- NA
plot(simrast_anim)
range(simrast0)

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


for(s in 1:n_frames) {

# s = 1
advance <- ggplot() +
  geom_spatraster(data = nonburnable_layer) +
  scale_fill_manual(values = "black",
                    na.value = "transparent", na.translate = F,
                    guide = guide_legend(order = 1)) +

  ggnewscale::new_scale_fill() +
  geom_spatraster(data = burnable_layer) +
  scale_fill_viridis(option = "A", discrete = TRUE, begin = 1,
                     na.value = "transparent",
                     na.translate = F,
                     guide = guide_legend(order = 2)) +

  ggnewscale::new_scale_fill() +
  geom_spatraster(data = sim_steps[[s]][["back"]], alpha = 1) +
  scale_fill_viridis(option = "A", discrete = TRUE, begin = 0.2,
                     na.value = "transparent",
                     na.translate = F,
                     guide = guide_legend(order = 3)) +

  ggnewscale::new_scale_fill() +
  geom_spatraster(data = sim_steps[[s]][["front"]], alpha = 1) +
  scale_fill_viridis(option = "A", discrete = TRUE, begin = 0.6,
                     na.value = "transparent",
                     na.translate = F,
                     guide = guide_legend(order = 4)) +

  # ignition point
  geom_spatvector(data = points3$balcon,
                  fill = "white", size = 1.2, color = "black",
                  shape = 21, stroke = 0.8) +

  annotation_scale(height = unit(1, "mm"),
                   text_col = "black", location = "tr", style = "ticks",
                   line_col = "black") +

  theme(legend.title = element_blank(),#element_text(size = 10, margin = margin(r = 4, unit = 'mm')),
        legend.position = "right",
        axis.text = element_text(size = 6),
        axis.text.x = element_text(angle = 50, hjust = 1),
        plot.title = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.spacing.y = unit(0.1, 'mm'),
        legend.margin = margin())+

  guides(fill = guide_legend(byrow = TRUE)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Incendio simulado")


ggsave(file.path("files", "images_rae_2023",
                 paste("fixed_model_maps_spread_anim_",
                       stringr::str_pad(s, 2, pad = "0"),
                       ".png", sep = "")),
       plot = advance, width = 12, height = 9, units = "cm")

}

