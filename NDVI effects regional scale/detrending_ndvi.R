# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_classic())
library(terra)
library(mgcv)
library(DHARMa)
library(mgcViz)

# Data --------------------------------------------------------------------

# vegetation names
dveg <- readxl::read_excel("/home/ivan/Insync/Mapa vegetación WWF - Lara et al. 1999/clases de vegetacion y equivalencias.xlsx",
                           sheet = "Sheet2")

# ndvi time series
v1 <- vect(file.path("..", "data", "NDVI_regional_data",
                     "ndvi_ts_detrend.shp"))

# images by summer used to compute ndvi
nn <- read.csv(file.path("..", "data", "NDVI_regional_data",
                         "ndvi_images-by-summer.csv"))


# Tidy data ----------------------------------------------------------------

names(dveg) <- c("vegnum1", "vegfac1", "vegnum2", "vegfac2")
dveg$vegfac1 <- factor(dveg$vegfac1, levels = dveg$vegfac1)
dveg$vegfac2 <- factor(dveg$vegfac2, levels = unique(dveg$vegfac2))

d1 <- as.data.frame(v1)
d1 <- left_join(d1, dveg, by = "vegnum1")

# remove non-burnable and burned
d2 <- d1[d1$burned == 0 & d1$vegnum2 < 6, grep("b_", colnames(d1))]
d2 <- d2[, order(colnames(d2))]

ndvi_summ <- data.frame(
  year = 1998:2022,
  mean = apply(as.matrix(d2), 2, mean),
  sd = apply(as.matrix(d2), 2, sd),
  n = nn$n
)

str(d2)

# Exploratory plots -------------------------------------------------------

ggplot(ndvi_summ, aes(year, mean)) +
  geom_point() + geom_smooth() +
  theme(panel.grid.major = element_line())

ggplot(ndvi_summ, aes(n, mean)) +
  geom_point() + geom_smooth() +
  ylim(0.6, 0.8)

ggplot(ndvi_summ, aes(year, n)) +
  geom_point() + geom_smooth()

ggplot(ndvi_summ, aes(year, sd)) +
  geom_point() + geom_smooth() +
  ylim(0.15, 0.20)

ggplot(ndvi_summ, aes(n, sd)) +
  geom_point() + geom_smooth()


# sd within years ranges between
round(range(ndvi_summ$sd), 4) # 0.0878 0.1099
# between years, is
round(sd(ndvi_summ$mean), 4) # 0.0199
# spatial variability is way much larger than temporal.
mean(ndvi_summ$sd) / sd(ndvi_summ$mean) # 4.851748
# spatial variation in NDVI IN FORESTS is 4.85 times the variation in time.
# Correct it?


# Detrend model ------------------------------------------------------------

# scale NDVI to [0, 1]
dmat <- as.matrix(d2)
dmat_01 <- (dmat + 1) / 2
ny <- ncol(dmat)

diff_mat_logit <- qlogis(dmat_01[, -ny]) - qlogis(dmat_01[, ny])
colnames(diff_mat_logit) <- 1998:(2022-1)
diff_mat_logit <- as.data.frame(diff_mat_logit)
diff_mat_logit$ndvi01_22 <- dmat_01[, ny]
diff_mat_logit$row_id <- 1:nrow(diff_mat_logit)

difflong <- pivot_longer(
  diff_mat_logit,
  all_of(which(colnames(diff_mat_logit) %in% as.character(1998:(2022-1)))),
  names_to = "year", values_to = "diff_logit"
)
difflong$year <- as.numeric(difflong$year)

# fit model to the logit-difference
md <- bam(diff_logit ~ te(year, ndvi01_22, k = c(6, 6), bs = "cr"),#s(year, k = 10, bs = "cr"),#
          data = difflong)
saveRDS(md, file.path("..", "data", "NDVI_regional_data", "ndvi_detrender_model.rds"))

# res <- simulateResiduals(md, n = 500)
# plot(res)
# plotResiduals(res, form = difflong$ndvi01_22)
# plotResiduals(res, form = difflong$year)

difflong$diff_fit <- fitted(md)
# widenize
diffwide <- pivot_wider(difflong[, c("diff_fit", "row_id", "year")],
                        names_from = "year", values_from = "diff_fit")
diff_fit_mat <- as.matrix(diffwide[, -1])

# detrend data
mat_logit <- qlogis(dmat_01[, -ny])
mat_logit_dt <- mat_logit - diff_fit_mat
mat_01_dt <- plogis(mat_logit_dt)
ndvi_dt <- mat_01_dt * 2 - 1
hist(ndvi_dt)
hist(dmat)

plot(density(ndvi_dt %>% as.vector, to = 1), ylim = c(0, 6))
lines(density(dmat[, -ny] %>% as.vector, to = 1), col = 2)
lines(density(dmat[, ny], to = 1), col = 3)

ndvi_ok <- cbind(ndvi_dt, dmat[, ny])

ndvi_summ_dt <- data.frame(
  year = 1998:2022,
  mean = apply(ndvi_ok, 2, mean),
  sd = apply(ndvi_ok, 2, sd),
  n = nn$n,
  dataset = "detrended"
)

ndvi_summ$dataset = "original"

ndvi_summ_both <- rbind(ndvi_summ, ndvi_summ_dt)

ggplot(ndvi_summ_both, aes(year, mean)) +
  geom_point() + geom_smooth() +
  facet_wrap(vars(dataset)) +
  theme(panel.grid.major = element_line())

ggplot(ndvi_summ_both, aes(year, sd)) +
  geom_point() + geom_smooth() +
  facet_wrap(vars(dataset)) +
  theme(panel.grid.major = element_line())

mviz <- mgcViz::getViz(md)
plot(sm(mviz, 1))

# detrend ndvi and fit the bayesian model later.

# get the year with the most average NDVI (ignoring 2022), so the non-year
# pixels are assigned that value
mm <- ndvi_summ$mean[-ny] %>% mean
year_mean <- ndvi_summ$year[which.min(abs(ndvi_summ$mean - mm))]
# 2003