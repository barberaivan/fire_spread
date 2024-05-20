# To define the final model, check which parameters vary across fires and
# whether that variation can be explained by FWI.

# Packages ----------------------------------------------------------------

library(tidyverse)
library(viridis)
library(bayesplot)
library(bayestestR) # hdi
library(posterior)
theme_set(theme_bw())

# Functions ---------------------------------------------------------------

# the same as tidy samples to be used with the result from
# rejection_sample_parallel
tidy_samples_ind <- function(samples_list) {

  nc <- length(samples_list)
  rr <- nrow(samples_list[[1]])
  cc <- ncol(samples_list[[1]])

  # turn into array of samples
  arr0 <- array(NA, dim = c(rr, cc, nc))
  for(c in 1:nc) {
    arr0[, , c] <- samples_list[[c]]
  }
  arr <- aperm(arr0, c(1, 3, 2))
  draws_arr <- as_draws_array(arr)

  draws_arr2 <- rename_variables(draws_arr,
                                 "intercept" = "...1",
                                 "vfi" = "...2",
                                 "tfi" = "...3",
                                 "slope" = "...4",
                                 "wind" = "...5",
                                 "steps" = "...6")

  return(draws_arr2)
}

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


# climatic data ------------------------------------------------------------

# fwi_data <- read.csv(file.path("data", "patagonian_fires_data_with_fwi.csv"))
fwi_data <- read.csv(file.path("data", "climatic_data_by_fire_fwi-fortnight-cumulative.csv"))
fwi_data <- fwi_data[, c("fire_id", "fwi_focal_fortnight", "fwi_expquad_fortnight",
                         "fwi_exp_fortnight", "area_ha")]

wind_data <- read.csv(file.path("data", "climatic_data_by_fire_FWI-wind_corrected.csv"))

clim_data <- left_join(fwi_data,
                       wind_data[, c("fire_id", "speed_mean")],
                       by = "fire_id")

clim_data$fire_id2 <- clim_data$fire_id # just to use synonyms

# Load posteriors ---------------------------------------------------------

# dir to load files
target_dir <- file.path("files", "pseudolikelihoods")

par_names <- c("intercept", "vfi", "tfi", "slope", "wind", "steps")

results_names <- list.files(target_dir)
samples_names <- results_names[grep("posterior_samples", results_names)]

# import in single df
samples_df <- do.call("rbind", lapply(samples_names, function(s) {
  # s <- samples_names[30]
  rrr <- readRDS(file.path(target_dir, s))
  fire_name <- strsplit(s, split = "-posterior_samples.rds")[[1]]

  if(is.data.frame(rrr)) {
    rrr$fire_id <- fire_name
    return(rrr)
  } else {
    rrr <- do.call("rbind", rrr)
    rrr <- as.data.frame(rrr)
    rrr$fire_id <- fire_name
    return(rrr)
  }
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



# Plots -------------------------------------------------------------------

# ~ fire
ggplot(posteriors_summ) +
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


# ~ fwi
ggplot(posteriors_summ) +
  geom_smooth(aes(fwi_expquad_fortnight, y = mean), method = "lm") +
  geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_95, ymax = hdi_upper_95), alpha = 0.7) +
  geom_linerange(aes(fwi_expquad_fortnight, ymin = hdi_lower_80, ymax = hdi_upper_80), alpha = 0.7) +
  geom_point(aes(fwi_expquad_fortnight, y = mean), alpha = 0.7) +
  facet_wrap(vars(par_name), ncol = 3, scales = "free_y") +
  xlab("FWI (expquad, fortnight)") +
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


# fwi ~ wind?
plot(fwi_expquad_fortnight ~ speed_mean, clim_data, xlim = c(0, 16))

# wind as a function of fwi and windspeed
mm <- lm(mean ~ speed_mean + fwi_expquad_fortnight,
         data = posteriors_summ[posteriors_summ$par_name == "wind", ])
summary(mm)

mm2 <- lm(mean ~ fwi_expquad_fortnight,
         data = posteriors_summ[posteriors_summ$par_name == "wind", ])
summary(mm2)

mm3 <- lm(mean ~ speed_mean,
          data = posteriors_summ[posteriors_summ$par_name == "wind", ])
summary(mm3)

# el efecto del viento disminuye con la vel media, como que lo compensa un poco.
# Pero no debe ser una compensación 1:1. chequear esto igual.
# Si! En el univariado es -0.9842.

wind_term = bwind * wind_pixel
bwind = (alpha + bb * wind_landscape)
windspeed_pixel ~ wind_landscape # si fuera casi 1:1

wind_term = (alpha + bb * wind_landscape) * wind_pixel
wind_term = alpha * wind_pixel + bb * wind_landscape * wind_pixel

# seguir pensando esto.


# wind issues -------------------------------------------------------------

dwind <- posteriors_summ[posteriors_summ$par_name == "wind", ]

# turn wind betas and means into factors relative to the mean.

dwind$wind_b_factor <- dwind$mean / mean(dwind$mean)
dwind$wind_s_factor <- dwind$speed_mean / mean(dwind$speed_mean)

wf1 <- lm(wind_b_factor ~ wind_s_factor, dwind)
summary(wf1)
plot(wind_b_factor ~ wind_s_factor, dwind)
abline(coef(wf1))
mean(dwind$wind_b_factor)
mean(dwind$wind_s_factor)
# tiene el mismo R2 que la reg de antes.
# es como si compensara a la mitad... o algo así.

ggplot(dwind, aes(wind_s_factor, wind_b_factor)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_abline(intercept = 2, slope = -1, linetype = "dashed", linewidth = 0.3, color = "red") +
  geom_smooth(method = "lm") +
  geom_point() +
  coord_fixed()
# Compensa, pero no perfectamente. Siempre hay más variación en los betas
# que en el windspeed. La slope es ~ 0.5, no 1.

ggplot(dwind, aes(speed_mean, mean)) +
  geom_hline(yintercept = mean(dwind$mean), linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = mean(dwind$speed_mean), linetype = "dashed", linewidth = 0.3) +
  geom_smooth(method = "lm") +
  geom_point() +
  coord_fixed()
# claro, acá la slope es 1, pero pero los betas varían el doble.

# si las escalo
dwind$beta_scaled <- as.numeric(scale(dwind$mean))
dwind$speed_scaled <- as.numeric(scale(dwind$speed_mean))

ggplot(dwind, aes(speed_scaled, beta_scaled)) +
  geom_hline(yintercept = mean(dwind$beta_scaled), linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = mean(dwind$speed_scaled), linetype = "dashed", linewidth = 0.3) +
  geom_abline(intercept = 0, slope = -1, linetype = "dashed", linewidth = 0.3, color = "red") +
  geom_smooth(method = "lm") +
  geom_point() +
  coord_fixed()

hist(dwind$wind_b_factor / dwind$wind_s_factor)
mean(dwind$wind_b_factor / dwind$wind_s_factor)

# pairwise quotients
qbeta <- outer(dwind$mean, dwind$mean, FUN = "/")
qspeed <- outer(dwind$speed_mean, dwind$speed_mean, FUN = "/")

qbeta <- qbeta[lower.tri(qbeta)]
qspeed <- qspeed[lower.tri(qspeed)]

plot(qbeta ~ qspeed)
curve(1 / x, add = T, col = 2)

qdf <- data.frame(beta = qbeta,
                  speed = qspeed)

ggplot(qdf, aes(speed, beta)) +
  geom_point(alpha = 0.5) +
  # geom_smooth(method = "gam", method.args = list(family = Gamma(link = "log"))) +
  geom_smooth(method = "gam", method.args = list(family = Gamma(link = "log")),
              formula = y ~ s(x, k = 50),
              color = "red", se = F) +
  # geom_smooth(method = "glm", method.args = list(family = Gamma(link = "log")),
  #             color = "pink", se = F) +
  geom_function(fun = function(x) 1 / x, colour = "blue", n = 201, linewidth = 1)

summary(qdf)

# tareas ------------------------------------------------------------------

# explorar las correlaciones. el wind parece variar mucho, pero en general
# estaba bastante correlated con el intercept.

rr <- abs(rnorm(500))
qrr <- outer(rr, rr, FUN = "/")
qrr <- qrr[lower.tri(qrr)]
summary(qrr)

# Posterior correlation within fires ---------------------------------------

fire_id <- unique(samples_df$fire_id)
nfires <- length(fire_id)

combs <- combn(par_names, 2) %>% t %>% as.data.frame()
combs$cor <- NA

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
  xlab("Correlation coefficient") +
  ylab("Number of fires")


# Parameters correlation among fires --------------------------------------

post_means <- aggregate(as.matrix(samples_df[, par_names_all]) ~ fire_id,
                        samples_df, mean)

post_means$fire_id2 <- post_means$fire_id
post_means$fire_id2[grep("2011_19", post_means$fire_id)] <- "2011_19"
post_means$fire_id2[grep("2015_47", post_means$fire_id)] <- "2015_47"

post_means <- left_join(post_means,
                        clim_data[, !(colnames(clim_data) %in% "fire_id")],
                        by = "fire_id2")

npar <- length(par_names)
par_names_log <- par_names
par_names_log[npar] <- par_names_all[npar+1]

# Raw correlations
GGally::ggpairs(post_means[, par_names_log])

# Correlations removing FWI effect

fff <- post_means[, par_names_log]
for(i in par_names_log) {
  # i = 1
  y <- post_means[, i, drop = T]
  x <- post_means$fwi_expquad_fortnight
  mm <- lm(y ~ x)
  fff[, i] <- y - fitted(mm)
}

GGally::ggpairs(fff)
