# Visualize separate likelihoods by fire to detect flat dimensions.
library(tidyverse); theme_set(theme_bw())
library(viridis)
library(mgcv)
library(deeptime)

# dir to save output
# target_dir <- file.path("files", "steps_optimization_not-scaled")
target_dir <- file.path("files", "steps_optimization")

# load file with constants to standardize
ndvi_params <- readRDS(file.path("data", "NDVI_regional_data",
                                 "ndvi_optim_and_proportion.rds"))

slope_sd <- ndvi_params$slope_term_sd

par_names <- c("forest", "shrubland", "grassland",
               "ndvi", "north", "elev", "slope", "wind",
               "steps")

par_fixed <- par_names[-length(par_names)]
n_par <- length(par_fixed)

ext_alpha <- 50
ext_beta <- 30

params_lower <- c("forest" = -ext_alpha,
                  "shrubland" = -ext_alpha,
                  "grassland" = -ext_alpha,
                  "ndvi" = -ext_beta,
                  "north" = 0,
                  "elev" = -ext_beta,
                  "slope" = 0,
                  "wind" = 0,
                  "steps" = 5)

params_upper <- c("forest" = ext_alpha,
                  "shrubland" = ext_alpha,
                  "grassland" = ext_alpha,
                  "ndvi" = 0,
                  "north" = ext_beta,
                  "elev" = 0,
                  "slope" = ext_beta / slope_sd, # because it is not standardized
                  "wind" = ext_beta,
                  "steps" = NA)


# Load simulations --------------------------------------------------------

sim_files <- list.files(target_dir, pattern = "_simulations.rds")
fire_names <- sapply(sim_files, function(x) {
  strsplit(x, "-steps_optimization_simulations.rds")[[1]][1]
}) %>% unname
sims <- lapply(sim_files, function(x) readRDS(file.path(target_dir, x)))
names(sims) <- fire_names
n_fires <- length(fire_names)

sims_df <- do.call("rbind", lapply(1:length(sims), function(i) {
  x <- sims[[i]]
  df <- data.frame(fire_id = fire_names[i],
                   overlap = x$overlap,
                   wave = x$wave)
  df2 <- cbind(df, as.data.frame(x$par_values))
  return(df2)
}))

sims_df$fire_id <- factor(sims_df$fire_id)


# Explore separate and joint overlap functions -----------------------------

models_list <- vector("list", n_par)
plots_list <- vector("list", n_par)
names(models_list) <- names(plots_list) <- par_fixed

for(p in 1:n_par) {
  # p = 7
  param <- par_names[p]
  print(param)

  # Get maximum overlaps in a fine sequence
  cut_points <- seq(params_lower[p], params_upper[p], length.out = 201)
  sims_df$xcat <- cut(sims_df[, param], breaks = cut_points,
                      include_lowest = T)

  dfsub <- sims_df[, c("overlap", param, "xcat", "fire_id")]
  names(dfsub)[2] <- "parvalues"
  dfagg <- aggregate(cbind(overlap, parvalues) ~ xcat + fire_id, dfsub, max)

  # Fit unidimensional GAM by fire
  nk <- 20
  kk <- seq(params_lower[p], params_upper[p], length.out = nk)

  mm <- bam(overlap ~
              fire_id +
              s(parvalues, by = fire_id, k = nk, bs = "cr"),
            data = dfagg, method = "fREML", knots = list(parvalues = kk),
            discrete = T, nthreads = 8)
  models_list[[p]] <- mm
  mm <- models_list[[p]]

  # data to predict
  ns <- 150
  dpred <- expand.grid(parvalues = seq(params_lower[p], params_upper[p],
                                       length.out = ns),
                       fire_id = levels(sims_df$fire_id))

  dpred$overlap <- predict(mm, dpred, se.fit = F)
  dpred$prediction <- "Overlap by fire"

  # widenize to compute log-sum

  ovmat <- matrix(dpred$overlap, ns, n_fires)
  dpred_wide <- dpred[1:ns, ]
  dpred_wide$fire_id <- "All"
  dpred_wide$prediction <- "Joint overlap"
  dpred_wide$overlap_logsum <- rowSums(log(ovmat), na.rm = T)
  op <- apply(ovmat, 1, prod)
  dpred_wide$overlap_prod <- op / sum(op)
  dpred_wide$overlap_mean <- rowMeans(ovmat)

  p1 <-
  ggplot(dpred, aes(parvalues, overlap, color = fire_id, group = fire_id)) +
    geom_line() +
    scale_color_viridis(discrete = T, option = "C", end = 0.9) +
    theme(legend.position = "none",
          panel.grid.minor = element_blank()) +
    ylim(0, 1) +
    scale_x_continuous(limits = c(params_lower[p], params_upper[p]),
                       breaks = seq(params_lower[p], params_upper[p],
                                    length.out = 5)) +
    ylab("Separate overlaps") +
    xlab(param)

  p2 <-
  ggplot(dpred_wide, aes(parvalues, overlap_logsum)) +
    geom_line() +
    ylab("Overlap log-sum") +
    xlab(param) +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(limits = c(params_lower[p], params_upper[p]),
                       breaks = seq(params_lower[p], params_upper[p],
                                    length.out = 5))

  p3 <-
  ggplot(dpred_wide, aes(parvalues, overlap_prod)) +
    geom_line() +
    theme(legend.position = "none",
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank()) +
    ylab("Overlaps product\n(normalized)") +
    xlab(param) +
    ylim(0, max(dpred_wide$overlap_prod)) +
    scale_x_continuous(limits = c(params_lower[p], params_upper[p]),
                       breaks = seq(params_lower[p], params_upper[p],
                                    length.out = 5))

  p4 <-
    ggplot(dpred_wide, aes(parvalues, overlap_mean)) +
    geom_line() +
    theme(legend.position = "none",
          axis.ticks.y = element_blank(),
          panel.grid.minor = element_blank()) +
    ylab("Overlaps mean") +
    xlab(param) +
    ylim(0, 1) +
    scale_x_continuous(limits = c(params_lower[p], params_upper[p]),
                       breaks = seq(params_lower[p], params_upper[p],
                                    length.out = 5))
  # plot(p4)

  plots_list[[p]] <- list(p1, p2, p3, p4)
}

# merge plots
kp <- 4
all_plots <- vector("list", n_par * kp)

idpar <- rep(1:n_par, each = kp)
idplot <- rep(1:kp, n_par)

for(i in 1:(n_par * kp)) {
  all_plots[[i]] <- plots_list[[idpar[i]]][[idplot[i]]]
}

x11()
pmerge_all <- ggarrange2(plots = all_plots, ncol = kp)
dev.off()

ggsave("files/steps_optimization/overlap_max_marginals_trial2.png",
       plot = pmerge_all,
       width = 30, height = 38, units = "cm")



# How many simulations to run? Will run again! ----------------------------

sims_df
head(sims_df)

sims_ord <- sims_df[order(sims_df$fire_id, sims_df$wave), ]
head(sims_ord)

n_waves <- max(sims_ord$wave)

ovmax <- matrix(NA, n_waves, n_fires)
ovmean <- matrix(NA, n_waves, n_fires)
steps_op <- matrix(NA, n_waves, n_fires)

for(f in 1:n_fires) {
  for(w in 1:n_waves) {
    # f = 4; w = 3
    rows <- which(sims_ord$fire_id == fire_names[f] &
                  sims_ord$wave <= w)
    ovmax[w, f] <- max(sims_ord$overlap[rows], na.rm = T)
    ovmean[w, f] <- mean(sims_ord$overlap[rows], na.rm = T)

    # get best step
    id_max_temp <- which.max(sims_ord$overlap[rows])
    id_max <- rows[id_max_temp]
    steps_op[w, f] <- sims_ord$steps[id_max]
  }
}

wave_summ <- data.frame(
  wave = rep(1:n_waves, n_fires),
  fire_id = rep(fire_names, each = n_waves),
  overlap_max = as.vector(ovmax),
  overlap_mean = as.vector(ovmean),
  steps_optim = round(as.vector(steps_op))
)

ggplot(wave_summ) +
  geom_line(aes(wave, overlap_mean), color = "gray") +
  geom_line(aes(wave, overlap_max), color = "black") +
  facet_wrap(vars(fire_id), ncol = 6) +
  theme(panel.grid.minor = element_blank())

ggplot(wave_summ) +
  geom_vline(xintercept = 60, color = "red") +
  geom_line(aes(wave, steps_optim), color = "black") +
  facet_wrap(vars(fire_id), ncol = 6, scales = "free_y") +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

# 74 is OK.


# NDVI and slope still at the limits --------------------------------------

# Correlations NDVI-intercepts??

# Subset best particles, as max - 0.1
thres <- 0.05

sims_best <- do.call("rbind", lapply(fire_names, function(f) {
  d <- sims_df[sims_df$fire_id == f, ]
  dbest <- d[d$overlap >= max(d$overlap - thres), ]
  return(dbest)
}))

## ggplot

ndvi_for <-
  ggplot(sims_best, aes(forest, ndvi)) +
  geom_point(alpha = 0.1) +
  facet_wrap(vars(fire_id), ncol = 7) +
  theme(panel.grid.minor = element_blank())
ggsave("files/steps_optimization/zz_plots/corr-ndvi-forest-all.png",
       plot = ndvi_for, width = 25, height = 32, units = "cm")

ndvi_shrub <-
  ggplot(sims_best, aes(shrubland, ndvi)) +
  geom_point(alpha = 0.1) +
  facet_wrap(vars(fire_id), ncol = 7) +
  theme(panel.grid.minor = element_blank())
ggsave("files/steps_optimization/zz_plots/corr-ndvi-shrubland-all.png",
       plot = ndvi_shrub, width = 25, height = 32, units = "cm")

ndvi_grass <-
  ggplot(sims_best, aes(grassland, ndvi)) +
  geom_point(alpha = 0.1) +
  facet_wrap(vars(fire_id), ncol = 7) +
  theme(panel.grid.minor = element_blank())
ggsave("files/steps_optimization/zz_plots/corr-ndvi-grassland-all.png",
       plot = ndvi_grass, width = 25, height = 32, units = "cm")

# slope-elev
slope_elev <-
  ggplot(sims_best, aes(elev, slope)) +
  geom_point(alpha = 0.1) +
  facet_wrap(vars(fire_id), ncol = 7) +
  theme(panel.grid.minor = element_blank())
ggsave("files/steps_optimization/zz_plots/corr-slope-elev-all.png",
       plot = slope_elev, width = 25, height = 32, units = "cm")



