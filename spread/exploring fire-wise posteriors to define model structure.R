# To define the final model, check which parameters vary across fires and
# whether that variation can be explained by FWI.

# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(bayesplot)
library(bayestestR) # hdi
library(posterior)
library(deeptime)
theme_set(theme_bw())

# Functions ---------------------------------------------------------------

summarise <- function(x) {
  q <- quantile(x, probs = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975),
                method = 8)
  names(q) <- c("eti_lower_95", "eti_lower_90", "eti_lower_80", "median",
                "eti_upper_80", "eti_upper_90", "eti_upper_95")

  hdi_95 <- hdi(x, ci = 0.95)
  hdi_90 <- hdi(x, ci = 0.90)
  hdi_80 <- hdi(x, ci = 0.80)

  hdis <- c(
    "hdi_lower_95" = hdi_95$CI_low, "hdi_lower_90" = hdi_90$CI_low,
    "hdi_lower_80" = hdi_80$CI_low,
    "hdi_upper_80" = hdi_80$CI_high, "hdi_upper_90" = hdi_90$CI_high,
    "hdi_upper_95" = hdi_95$CI_high
  )

  res <- c("mean" = mean(x), hdis, q)

  return(res)
}

unconstrain <- function(x, support) {
  xun <- x

  names_logit <- colnames(x)[colnames(x) != "steps"]
  for(j in names_logit) {
    xun[, j] <- qlogis((x[, j] - support[1, j]) / (support[2, j] - support[1, j]))
  }

  xun[, "steps"] <- log(x[, "steps"])

  return(xun)
}

# climatic data ------------------------------------------------------------

# fwi_data <- read.csv(file.path("data", "patagonian_fires_data_with_fwi.csv"))
fwi_data_fort <- read.csv(file.path("data", "climatic_data_by_fire_fwi-fortnight-cumulative.csv"))
fwi_data_day <- read.csv(file.path("data", "climatic_data_by_fire_fwi-day-cumulative.csv"))

fwi_data <- left_join(
  fwi_data_day[, c("fire_id", "fwi_focal_day", "fwi_expquad_day", "fwi_exp_day")],
  fwi_data_fort[, c("fire_id", "fwi_focal_fortnight", "fwi_expquad_fortnight",
                    "fwi_exp_fortnight", "area_ha")],
  by = "fire_id"
)

wind_data <- read.csv(file.path("data", "climatic_data_by_fire_FWI-wind_corrected.csv"))

clim_data <- left_join(fwi_data,
                       wind_data[, c("fire_id", "speed_mean")],
                       by = "fire_id")

clim_data$fire_id2 <- clim_data$fire_id # just to use synonyms


# A few constants ---------------------------------------------------------

par_names <- c("intercept", "vfi", "tfi", "slope", "wind", "steps")
n_coef <- length(par_names)
# load file with constants to standardize
ndvi_params <- readRDS(file.path("data", "flammability indices",
                                 "ndvi_optim_and_proportion.rds"))
slope_sd <- ndvi_params$slope_term_sd

ext_alpha <- 50
ext_beta <- 30

params_lower <- c(-ext_alpha, rep(0, n_coef-2), 5)
params_upper <- c(ext_alpha, rep(ext_beta, n_coef-2), NA)
names(params_lower) <- names(params_upper) <- par_names
params_upper["slope"] <- ext_beta / slope_sd

support <- rbind(params_lower, params_upper)
colnames(support) <- names(params_lower) <- names(params_upper) <- par_names
support_width <- apply(support, 2, diff)

# Load posteriors ---------------------------------------------------------

# dir to load files
target_dir <- file.path("files", "posterior_samples_stage1")

samples_files <- list.files(target_dir, pattern = "-samples.rds")

# import in single df
samples_df <- do.call("rbind", lapply(samples_files, function(s) {
  # s <- samples_files[30]
  rrr <- readRDS(file.path(target_dir, s))

  # resample to get 10000
  ids <- sample(1:nrow(rrr$samples), 10000, prob = rrr$samples$prob,
                replace = T)

  r <- cbind(
    data.frame(fire_id = rrr$fire_id),
    as.data.frame(rrr$samples$par_values[ids, ])
  )

  return(r)
}))

# add steps at log scale
samples_df$steps_log <- log(samples_df$steps)
par_names_all <- c(par_names, "steps_log")

# longanize to summarize
samples_long <- pivot_longer(samples_df,
                             all_of(which(names(samples_df) %in% par_names_all)),
                             names_to = "par_name", values_to = "par_value")

posteriors_summ <- aggregate(par_value ~ par_name + fire_id, samples_long, summarise)

posteriors_summ <- cbind(posteriors_summ[, c("par_name", "fire_id")],
                         as.data.frame(posteriors_summ$par_value))

# merge with FWI data

# a few fires were divided for spread, but have the same data
posteriors_summ$fire_id2 <- posteriors_summ$fire_id
posteriors_summ$fire_id2[grep("2011_19", posteriors_summ$fire_id)] <- "2011_19"
posteriors_summ$fire_id2[grep("2015_47", posteriors_summ$fire_id)] <- "2015_47"

posteriors_summ <- left_join(posteriors_summ,
                             clim_data[, !(colnames(clim_data) %in% "fire_id")],
                             by = "fire_id2")

# tidy factors
posteriors_summ$par_name <- factor(posteriors_summ$par_name,
                                   levels = par_names_all)

fires_unique <- posteriors_summ[!duplicated(posteriors_summ[, c("fire_id", "area_ha")]), ]
fires_unique <- fires_unique[order(fires_unique$area_ha), ]

# View(fires_unique)
posteriors_summ$fire_id <- factor(posteriors_summ$fire_id,
                                  levels = fires_unique$fire_id)

# add log area
posteriors_summ$area_ha_log <- log10(posteriors_summ$area_ha)


# Remove flat posteriors --------------------------------------------------

pars_check <- par_names[-length(par_names)]

posteriors_summ$too_wide <- FALSE

for(i in 1:nrow(posteriors_summ)) {
  local_par <- posteriors_summ$par_name[i]
  if(local_par %in% pars_check) {
    ic_width <- posteriors_summ$hdi_upper_95[i] - posteriors_summ$hdi_lower_95[i]
    posteriors_summ$too_wide[i] <- ic_width >= (0.90 * support_width[local_par])
  }
}

# Plots -------------------------------------------------------------------

# ~ fire
ggplot(posteriors_summ[!posteriors_summ$too_wide, ]) +
  geom_linerange(aes(fire_id, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  geom_linerange(aes(fire_id, ymin = hdi_lower_80, ymax = hdi_upper_80), alpha = 0.7) +
  geom_point(aes(fire_id, y = mean), alpha = 0.7) +
  facet_wrap(vars(par_name), ncol = 2, scales = "free_y") +
  theme(axis.text.x = element_blank())

# ~ area
ggplot(posteriors_summ) +
  geom_smooth(aes(area_ha_log, y = mean), method = "lm") +
  geom_linerange(aes(area_ha_log, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  geom_linerange(aes(area_ha_log, ymin = hdi_lower_80, ymax = hdi_upper_80), alpha = 0.7) +
  geom_point(aes(area_ha_log, y = mean), alpha = 0.7) +
  facet_wrap(vars(par_name), ncol = 3, scales = "free_y") +
  xlab("Fire size (log-ha)") +
  ylab("Estimates") +
  theme(strip.text = element_text(vjust = -0.5),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_blank())


# ~ fwi fortnight
ggplot(posteriors_summ) +#[!posteriors_summ$too_wide, ]) +
  geom_smooth(aes(fwi_expquad_fortnight, y = mean), method = "lm") +
  geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  # geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_80, ymax = hdi_upper_80), alpha = 0.7) +
  geom_point(aes(fwi_expquad_fortnight, y = mean), alpha = 0.7) +
  facet_wrap(vars(par_name), ncol = 3, scales = "free_y") +
  xlab("FWI (expquad, fortnight)") +
  ylab("Estimates") +
  theme(strip.text = element_text(vjust = -0.5),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_blank())

# quadratic
ggplot(posteriors_summ[!posteriors_summ$too_wide, ]) +
  geom_smooth(aes(fwi_expquad_fortnight, y = mean), method = "lm",
              formula = "y ~ x + I(x ^ 2)") +
  geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_80, ymax = hdi_upper_80), alpha = 0.7) +
  geom_point(aes(fwi_expquad_fortnight, y = mean), alpha = 0.7) +
  facet_wrap(vars(par_name), ncol = 3, scales = "free_y") +
  xlab("FWI (expquad, fortnight)") +
  ylab("Estimates") +
  theme(strip.text = element_text(vjust = -0.5),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_blank())

# GAM
ggplot(posteriors_summ[!posteriors_summ$too_wide, ]) +
  geom_smooth(aes(fwi_expquad_fortnight, y = mean), method = "gam") +
  geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_80, ymax = hdi_upper_80), alpha = 0.7) +
  geom_point(aes(fwi_expquad_fortnight, y = mean), alpha = 0.7) +
  facet_wrap(vars(par_name), ncol = 3, scales = "free_y") +
  xlab("FWI (expquad, fortnight)") +
  ylab("Estimates") +
  theme(strip.text = element_text(vjust = -0.5),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_blank())

# log(x)
ggplot(posteriors_summ[!posteriors_summ$too_wide, ]) +
  geom_smooth(aes(fwi_expquad_fortnight, y = mean), method = "lm",
              formula = "y ~ log(x)") +
  geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_80, ymax = hdi_upper_80), alpha = 0.7) +
  geom_point(aes(fwi_expquad_fortnight, y = mean), alpha = 0.7) +
  facet_wrap(vars(par_name), ncol = 3, scales = "free_y") +
  xlab("FWI (expquad, fortnight)") +
  ylab("Estimates") +
  theme(strip.text = element_text(vjust = -0.5),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_blank())


# ~ fwi day
ggplot(posteriors_summ[!posteriors_summ$too_wide, ]) +
  geom_smooth(aes(fwi_expquad_day, y = mean), method = "lm",
              formula = "y ~ x + I(x ^ 2)") +
  geom_linerange(aes(fwi_expquad_day, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  geom_linerange(aes(fwi_expquad_day, ymin = hdi_lower_80, ymax = hdi_upper_80), alpha = 0.7) +
  geom_point(aes(fwi_expquad_day, y = mean), alpha = 0.7) +
  facet_wrap(vars(par_name), ncol = 3, scales = "free_y") +
  xlab("FWI (expquad, day)") +
  ylab("Estimates") +
  theme(strip.text = element_text(vjust = -0.5),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_blank())

ggplot(posteriors_summ[!posteriors_summ$too_wide, ]) +
  geom_smooth(aes(fwi_expquad_day, y = mean), method = "lm") +
  geom_linerange(aes(fwi_expquad_day, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  geom_linerange(aes(fwi_expquad_day, ymin = hdi_lower_80, ymax = hdi_upper_80), alpha = 0.7) +
  geom_point(aes(fwi_expquad_day, y = mean), alpha = 0.7) +
  facet_wrap(vars(par_name), ncol = 3, scales = "free_y") +
  xlab("FWI (expquad, day)") +
  ylab("Estimates") +
  theme(strip.text = element_text(vjust = -0.5),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_blank())


# ~ windspeed
ggplot(posteriors_summ) +
  geom_smooth(aes(speed_mean, y = mean), method = "lm") +
  geom_linerange(aes(speed_mean, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  geom_linerange(aes(speed_mean, ymin = hdi_lower_80, ymax = hdi_upper_80), alpha = 0.7) +
  geom_point(aes(speed_mean, y = mean), alpha = 0.7) +
  facet_wrap(vars(par_name), ncol = 3, scales = "free_y") +
  xlab("Wind speed") +
  ylab("Estimates") +
  theme(strip.text = element_text(vjust = -0.5),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_blank())
# clear effect of windspeed. the low data-speed is compensated with a
# large wind effect.

# Intercepts ~ FWI

# ~ fwi fortnight
ggplot(posteriors_summ[!posteriors_summ$too_wide &
                         posteriors_summ$par_name %in% par_names[1:5], ]) +
  geom_smooth(aes(fwi_expquad_fortnight, y = mean,
                  group = par_name, color = par_name, fill = par_name),
              method = "lm",
              linewidth = 0.7, alpha = 0.2) +
  # geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  # geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_80, ymax = hdi_upper_80), alpha = 0.7) +
  geom_point(aes(fwi_expquad_fortnight, y = mean,
                 group = par_name, color = par_name),
             alpha = 0.7) +
  # facet_wrap(vars(par_name), ncol = 3, scales = "free_y") +
  xlab("FWI (expquad, fortnight)") +
  ylab("Estimates") +
  theme(strip.text = element_text(vjust = -0.5),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_blank())


# Slope
ggplot(posteriors_summ[!posteriors_summ$too_wide &
                         posteriors_summ$par_name == "slope", ]) +
  geom_smooth(aes(fwi_expquad_fortnight, y = mean,
                  group = par_name, color = par_name, fill = par_name),
              method = "lm",
              linewidth = 0.7, alpha = 0.2) +
  geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  geom_point(aes(fwi_expquad_fortnight, y = mean,
                 group = par_name, color = par_name),
             alpha = 0.7) +
  # facet_wrap(vars(par_name), ncol = 3, scales = "free_y") +
  xlab("FWI (expquad, fortnight)") +
  ylab("Estimates") +
  theme(strip.text = element_text(vjust = -0.5),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_blank())


# Compare FWI effect logged or not

steps_data <- posteriors_summ[posteriors_summ$par_name == "steps_log", ]

lm1 <- lm(mean ~ fwi_expquad_fortnight, data = steps_data)
lm2 <- lm(mean ~ log(fwi_expquad_fortnight), data = steps_data)
summary(lm1)
summary(lm2) # Much better the linear effect.



# The same at unconstrained scale -----------------------------------------

# import in single df
samples_df_unc <- do.call("rbind", lapply(samples_files, function(s) {
  # s <- samples_files[30]
  rrr <- readRDS(file.path(target_dir, s))

  # resample to get 10000
  ids <- sample(1:nrow(rrr$samples), 10000, prob = rrr$samples$prob,
                replace = T)
  xx <- unconstrain(rrr$samples$par_values[ids, ], support = support)

  # remove non-finite
  keep_id <- apply(xx, 1, function(x) all(is.finite(x)))

  r <- cbind(
    data.frame(fire_id = rrr$fire_id),
    as.data.frame(xx[keep_id, ])
  )

  return(r)
}))

# longanize to summarize
samples_unc_long <- pivot_longer(samples_df_unc,
                             all_of(which(names(samples_df) %in% par_names)),
                             names_to = "par_name", values_to = "par_value")

posteriors_summ_unc <- aggregate(par_value ~ par_name + fire_id,
                                 samples_unc_long, summarise)

posteriors_summ_unc <- cbind(posteriors_summ_unc[, c("par_name", "fire_id")],
                             as.data.frame(posteriors_summ_unc$par_value))

# merge with FWI data

# a few fires were divided for spread, but have the same data
posteriors_summ_unc$fire_id2 <- posteriors_summ_unc$fire_id
posteriors_summ_unc$fire_id2[grep("2011_19", posteriors_summ_unc$fire_id)] <- "2011_19"
posteriors_summ_unc$fire_id2[grep("2015_47", posteriors_summ_unc$fire_id)] <- "2015_47"

posteriors_summ_unc <- left_join(posteriors_summ_unc,
                             clim_data[, !(colnames(clim_data) %in% "fire_id")],
                             by = "fire_id2")

# tidy factors
posteriors_summ_unc$par_name <- factor(posteriors_summ_unc$par_name,
                                   levels = par_names_all)

fires_unique <- posteriors_summ_unc[!duplicated(posteriors_summ_unc[, c("fire_id", "area_ha")]), ]
fires_unique <- fires_unique[order(fires_unique$area_ha), ]

# View(fires_unique)
posteriors_summ_unc$fire_id <- factor(posteriors_summ_unc$fire_id,
                                  levels = fires_unique$fire_id)


# ~ fwi fortnight
ggplot(posteriors_summ_unc) +#[!posteriors_summ$too_wide, ]) +
  geom_smooth(aes(fwi_expquad_fortnight, y = mean), method = "lm") +
  geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  # geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_80, ymax = hdi_upper_80), alpha = 0.7) +
  geom_point(aes(fwi_expquad_fortnight, y = mean), alpha = 0.7) +
  facet_wrap(vars(par_name), ncol = 3, scales = "free_y") +
  xlab("FWI (expquad, fortnight)") +
  ylab("Estimates") +
  theme(strip.text = element_text(vjust = -0.5),
        strip.background = element_rect(fill = "white", color = "white"),
        panel.grid.minor = element_blank())


# Posterior correlation within fires ---------------------------------------

fire_id <- unique(samples_df$fire_id)
nfires <- length(fire_id)

combs <- combn(par_names, 2) %>% t %>% as.data.frame()
combs$cor <- NA
ncombs <- nrow(combs)

cor_data <- do.call("rbind", lapply(fire_id, function(ff) {
  cc <- cov2cor(cov(samples_df[samples_df$fire_id == ff, par_names]))
  for(i in 1:ncombs) {
    combs$cor[i] <- cc[combs$V2[i], combs$V1[i]]
  }
  combs$fire_id <- ff
  return(combs)
}))

cor_data$V1 <- factor(cor_data$V1, levels = par_names)
cor_data$V2 <- factor(cor_data$V2, levels = par_names[-1])

ggplot(cor_data, aes(cor)) +
  geom_histogram(binwidth = 0.5, breaks = seq(-1, 1, by = 1 / 3),
                 color = "black", fill = "black", linewidth = 0.25,
                 alpha = 0.6) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  facet_grid(V2 ~ V1) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "white"),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10)) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  xlab("Correlation coefficient") +
  ylab("Number of fires")


# Parameters correlation among fires --------------------------------------


# Hay que hacer esto a mano, porque para cada par de parámetros quiero quitar
# algunos datos, pero no la fila completa.

combs <- combn(par_names, 2) %>% t %>% as.data.frame()
ncombs <- nrow(combs)

cordata <- posteriors_summ_unc#[!posteriors_summ$too_wide, ]




# fwi model
mm <- lm(mean ~ par_name * fwi_expquad_fortnight, data = cordata)
cordata$mean_res <- cordata$mean - fitted(mm)

# Create data list with paired variables

dl <- vector("list", ncombs)
plots_list_marg <- vector("list", ncombs)
plots_list_cond <- vector("list", ncombs)

for(cc in 1:ncombs) {
  # cc = 1
  var_x <- combs[cc, "V1"]
  var_y <- combs[cc, "V2"]

  dsub <- cordata[cordata$par_name %in% c(var_x, var_y),
                  c("par_name", "fire_id", "mean", "mean_res")]

  dx <- dsub[dsub$par_name == var_x, ]
  colnames(dx) <- c("par_name_x", "fire_id", "mean_x", "mean_res_x")

  dy <- dsub[dsub$par_name == var_y, ]
  colnames(dy) <- c("par_name_y", "fire_id", "mean_y", "mean_res_y")

  dboth <- inner_join(dx, dy, by = "fire_id")

  dl[[cc]] <- dboth

  # make plots

  plots_list_marg[[cc]] <- ggplot(dboth, aes(mean_x, mean_y)) +
    geom_smooth(method = "lm", linewidth = 0.5) +
    geom_point(alpha = 0.7) +
    theme(panel.grid.minor = element_blank()) +
    xlab(var_x) + ylab(var_y)

  plots_list_cond[[cc]] <- ggplot(dboth, aes(mean_res_x, mean_res_y)) +
    geom_smooth(method = "lm", linewidth = 0.5) +
    geom_point(alpha = 0.7) +
    theme(panel.grid.minor = element_blank()) +
    xlab(var_x) + ylab(var_y)
}


lay <- matrix(
  c(
    1, 2, 3, 4, 5,
    NA, 6, 7, 8, 9,
    NA, NA, 10, 11, 12,
    NA, NA, NA, 13, 14,
    NA, NA, NA, NA, 15
  ),
  ncol = 5
)

plot_cor_marg <- ggarrange2(plots = plots_list_marg,
                            layout = lay)
ggsave("hierarchical model/pairs_plot_marginal.png", plot = plot_cor_marg,
       width = 20, height = 20, units = "cm")


plot_cor_cond <- ggarrange2(plots = plots_list_cond,
                            layout = lay)
ggsave("hierarchical model/pairs_plot_conditional.png", plot = plot_cor_cond,
       width = 20, height = 20, units = "cm")
