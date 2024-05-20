# To reduce the dimensionality of predictors in the spread model we build
# a vegetation flammability index using a logistic regression
# that predicts burned vs unburned pixels as a function of vegetation.

# Based on Barberá et al. 2024, we found that topographic variables are somewhat
# correlated, and including an elevation term and the directional slope effect
# would be enough. Hence, we include those variables in the vegetation model to
# remove their effects only. The vegetation terms are quadratic functions of the
# NDVI by vegetation type, with an intercept by vegetation type.

# We use the same data as in Barberá et al. 2024, with a very similar but reduced
# model.

# ndvi_max is ndvi_max_max, the maximum of all summer maximums.
# ndvi_mean is ndvi_mean_max, the maximum of all summer means.

############## above this, old comment

# Ideas:
# probarel efecto del NDVI a secas sobre la burn prob, a ver qué onda.
# Ignorar el agua pero considerar lo no quemable.

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_minimal())
library(terra)
library(rstan)
library(posterior)
library(mgcv)
library(DHARMa)
library(HDInterval)

# Functions ---------------------------------------------------------------

# to compute mean and 95 % CI from samples.
mean_ci <- function(x, name = "p_") {
  qs <- quantile(x, probs = c(0.025, 0.975), method = 8)
  tmp <- c(mean(x), qs) %>% unname
  names(tmp) <- paste(name, c("mean", "lower", "upper"), sep = "")
  return(tmp)
}

mean_hdi <- function(x, name = "p_", prob = 0.95) {
  hh <- hdi(x, prob)
  tmp <- c(mean(x), hh["lower"], hh["upper"])
  names(tmp) <- paste(name, c("mean", "lower", "upper"), sep = "")
  return(tmp)
}

# Data --------------------------------------------------------------------

v0 <- vect(file.path("..", "data", "flammability indexes",
                     "flammability_indexes_points.shp"))

dveg <- readxl::read_excel("/home/ivan/Insync/Mapa vegetación WWF - Lara et al. 1999/clases de vegetacion y equivalencias.xlsx",
                           sheet = "Sheet2")

names(dveg) <- c("vegnum1", "vegfac1", "vegnum2", "vegfac2")
dveg$vegfac1 <- factor(dveg$vegfac1, levels = dveg$vegfac1)
dveg$vegfac2 <- factor(dveg$vegfac2, levels = unique(dveg$vegfac2))

names(v0)[names(v0) == "veg"] <- "vegnum1"
names(v0)[names(v0) == "elev"] <- "elevation"


# data to chech NDVI trend through time
v1 <- vect(file.path("..", "data", "flammability indexes",
                     "flammability_indexes_points_ndvi_ts_check.shp"))

# images by summer used to compute ndvi
nn <- read.csv(file.path("..", "data", "flammability indexes",
                         "ndvi_images-by-summer.csv"))


# Tidy data ---------------------------------------------------------------

# remove water and ice and snow
data <- as.data.frame(v0)
data <- left_join(data, dveg, by = "vegnum1")
agua <- c(9, 11)
data <- data[!(data$vegnum1 %in% agua), ]

# compute slope-weighted northing
# (see code <northing importance function.R>)
data$northing <- cos(data$aspect * pi / 180) * plogis(-5 + 0.35 * data$slope)


d1 <- as.data.frame(v1)
d1 <- left_join(d1, dveg, by = "vegnum1")
# remove non-burnable and burned
d2 <- d1[d1$burned == 0 & d1$vegnum2 %in% c(1, 2), grep("b_", colnames(d1))]
d2 <- d2[, order(colnames(d2))]

ndvi_summ <- data.frame(
  year = 1998:2022,
  mean = apply(as.matrix(d2), 2, mean),
  sd = apply(as.matrix(d2), 2, sd),
  n = nn$n
)

ggplot(ndvi_summ, aes(year, mean)) +
  geom_point() + geom_smooth() +
  ylim(0.7, 0.85)

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

# NDVI should be harmonized across sensors.
# Or translate to the 2022 values using a GAM. Model the logit-diff with 2022
# as a function of year and 2022-ndvi. Then, add this logit-diff and take the
# inv_logit.


# Exploratory -------------------------------------------------------------

ggplot(data, aes(ndvi_max, ndvi_dyn)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam",
              formula = y ~ s(x, k = 10, bs = "cr"),
              method.args = list(family = mgcv::betar())) +
  geom_abline(slope = 1, intercept = 0) +
  coord_fixed()
# quite similar
