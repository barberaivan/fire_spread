library(tidyverse)
library(randtoolbox)   # sobol sequences
library(mgcv)          # Fit Gaussian Processes - difícil de definir
library(GauPro)        # Fit Gaussian Processes - muy complicado y sin modelado para la media
library(mlegp)         # Fit Gaussian Processes - sin modelado para la media (aparentemente)
library(brms)          # Fit Gaussian Processes - tarda demasiado, y peor con muchas partículas
library(DiceKriging)   # Fit Gaussian Processes - Parece muy bueno! permite definir el ruido de las simulaciones estocásticas.

library(Rcpp)
sourceCpp("spread_functions.cpp")


# A few constants ---------------------------------------------------------

d <- 9 # likelihood dimension

# Exploring the sobol sequence ----------------------------------------

m1 <- sobol(n = 1000, dim = d, init = T) # use init = F to grow the seqence.
m2 <- sobol(n = 1000, dim = d, init = F)
m3 <- sobol(n = 1000, dim = d, init = F)

# plot(m1[, 1], m1[, 2])
# points(m2[, 1], m2[, 2], col = "red")
# points(m3[, 1], m3[, 2], col = "blue")

# If the first sample (m1) is initialized, all the following will always be the
# same. i.e., the sequence is fixed. init = F is needed so the following sequences
# are not the same as the first.


# simulate prior function -------------------------------------------------

prior_sim <- function(sd_int = 2, mu_int = 0, r_01 = 0.1, r_z = 0.25) {
  betas <- c(
    "intercept" = rnorm(1, mu_int, sd_int),   # shrubland logit (reference class)
    "subalpine" = rnorm(1, 0, sd_int),        # veg coefficients
    "wet" = rnorm(1, 0, sd_int),
    "dry" = rnorm(1, 0, sd_int),
    "fwi" = rexp(1, r_z),                     # positive
    "aspect" = rexp(1, r_01),                 # positive (northing)
    "wind" = rexp(1, r_01),                   # positive
    "elevation" = (-1) * rexp(1, r_z),        # negative
    "slope" = rexp(1, r_01)                   # positive
  )

  return(betas)
}

prior_sim()

# function to compute prior quantiles from the percentiles, which will be
# obtained as sobols sequences in [0, 1] ^ d
prior_q <- function(p, mu_int = 0, sd_int = 2, r_01 = 0.1, r_z = 0.25) {
  q <- matrix(NA, nrow(p), ncol(p))
  colnames(q) <- c("intercept", "subalpine", "wet", "dry",
                   "fwi", "aspect", "wind", "elevation", "slope")

  q[, 1] <- qnorm(p[, 1], mean = mu_int, sd = sd_int)

  for(i in 2:4) {
    q[, i] <- qnorm(p[, i], mean = 0, sd = sd_int)
  }

  names_01 <- which(colnames(q) %in% c("aspect", "wind", "slope"))
  for(i in names_01) { # [0, 1] predictors
    q[, i] <- qexp(p[, i], rate = r_01)
  }

  names_z <- which(colnames(q) %in% c("fwi", "elevation"))
  for(i in names_z) { # standardized predictors
    q[, i] <- qexp(p[, i], rate = r_z)
  }

  return(q)
}


n = 500
p <- sobol(n, d)
w0 <- prior_q(p)

pairs(w0, pch = 19, col = rgb(0, 0, 0, 0.05))



# gp trials ----------------------------------------------------------------

y <- data.frame(y = rnorm(n))

data <- cbind(as.data.frame(w0), y)

mod <- gam(y ~ s(intercept, k = 6, bs = "gp") +
               s(subalpine, k = 6, bs = "gp") +
               s(wet, k = 6, bs = "gp") +
               s(dry, k = 6, bs = "gp") +
               s(fwi, k = 6, bs = "gp") +
               s(aspect, k = 6, bs = "gp") +
               s(wind, k = 6, bs = "gp") +
               s(elevation, k = 6, bs = "gp") +
               s(slope, k = 6, bs = "gp"),
           data = data, method = "REML")

# it's difficult to fit a gp because all interactions should be considered.
summary(mod)

# It's easier to use a GP directly. We need one that can fit a nugget term, to
# take into account unexplained variability.
# In addition, raw results could be used to account for heteroscedasticity and
# model it. However, that would be more complex.
# Try a first wave and compare the results between 2 gp:
#   including or not a quadratic term for the mean.
# Perhaps the best package by now is GauPro. GPfit is only for deterministic
# simulators. (it would be useful to parameterize FARSITE!)



# GP trials with GauPro ---------------------------------------------------

diamonds_subset <- diamonds[sample(1:nrow(diamonds), 60), ]
dm <- gpkm(price ~ carat + cut + color + clarity + depth,
           diamonds_subset,
           kernel = "gauss")
summary(dm)

predict(dm, diamonds_subset)
dm$predict(diamonds_subset)

names(dm)
dm$plotLOO()
dm$fit()

?GauPro_kernel_model


curve((10) * x + (-10) * x ^ 2)

# with my data
colnames(p) <- colnames(w0)
y <- rnorm(nrow(p),
           (10) * p[, 1] + (-10) * p[, 1] ^ 2,
           sd = 1)

mod1 <- gpkm(p, y, kernel = "gauss")
summary(mod1)

y_hat <- mod1$predict(p, se.fit = TRUE)
head(y_hat)

table_plot <- cbind(
  data.frame(y = y, y_hat = y_hat),
  p
)

ggplot(table_plot, aes(x = intercept, y = y)) +
  geom_point()

ggplot(table_plot, aes(x = intercept, y = y_hat)) +
  geom_point()



# 1 dimension, to tell the difference between se and s2
gp <- GauPro(p[, 1, drop = F], y)

plot(p[, 1], y)
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x)+2*gp$predict(x, se=T)$se, add=T, col=4)
curve(gp$predict(x)-2*gp$predict(x, se=T)$se, add=T, col=4)
curve(gp$predict(x)+2*gp$predict(x, se=T)$s2, add=T, col=3)
curve(gp$predict(x)-2*gp$predict(x, se=T)$s2, add=T, col=3)

gp$s2_hat


all.equal(sqrt(ppp$s2), ppp$se)


N <- 900
x <- c(runif(N / 3, 0, 0.1),
       runif(N / 3, 0.45, 0.55),
       runif(N / 3, 0.9, 1))
y <- rnorm(N,
           (10) * x + (-10) * x ^ 2,
           sd = 1 + 2 * x)


tr <- trend_LM$new(D=1)
gp <- GauPro(x, y, trend = tr, parallel = TRUE, parallel_cores = 10)

gp$cool1Dplot()
# it is not heteroscedastic

gp$plot()

gp$cool_plot_1d()
gp$Z

xseq <- seq(0, 1, length.out = 150)
ppp <- gp$predict(xseq, se.fit= T)
head(ppp)

hist(ppp$s2) # very very tight both plots.
hist(ppp$se)

# why does it vary?
gp$s2_hat # doesn't make sense...
# the se from predict does make sense

gp$param.est
gp$corr_func()
gp$nug

gp$summary()

?trend_LM

# GauPro seems incredible, but I can't fit a quadratic linear model to the trend,
# and the documentations is terrible.


# GP with brms ------------------------------------------------------------

n = 500
p <- sobol(n, d)
w0 <- prior_q(p)

d0 <- as.data.frame(w0)
head(d0)

d0$y <- rnorm(nrow(w0),
              (10) * d0$intercept + (-10) * d0$intercept ^ 2,
              sd = 1)

# fit gp
# m1 <- brm(y ~ intercept + I(intercept ^ 2) + gp(intercept), data = d0)
# stancode(m1)
# es lentísimo, pero MUUUUUY!!



# with funGP --------------------------------------------------------------
# no support for mean modelling



# DiceKriging -------------------------------------------------------------

n = 20
p <- sobol(n, d)
w0 <- prior_q(p)

colnames(p) <- colnames(w0)
d0 <- as.data.frame(p)
d0$y <- rnorm(nrow(d0),
              (10) * d0$intercept + (-10) * d0$intercept ^ 2,
              sd = 1)
d0$noise <- exp((-3) * d0$intercept + (3) * d0$intercept ^ 2)

# plot(noise ~ intercept, d0)

des <- d0[, "intercept", drop = F]

m2 <- km(y ~ intercept + I(intercept ^ 2), design = des, response = d0$y,
         nugget.estim = TRUE)

m4 <- km(y ~ intercept + I(intercept ^ 2), design = des, response = d0$y,
         noise.var = d0$noise)

# m3 <- km(~1, design = des, response = d0$y,
#          nugget.estim = TRUE)

mod <- m4
dpred <- data.frame(intercept = seq(-2, 2, length.out = 100))
preds <- predict(mod, dpred, type = "UK")
dpred$trend <- as.numeric(preds$trend)
dpred$mean <- preds$mean
dpred$lower95 <- preds$lower95
dpred$upper95 <- preds$upper95

plot(mean ~ trend, dpred); abline(0, 1) # mean and trend are the same here.

plot(mean ~ intercept, data = dpred, type = "l", ylim = c(-10, 10))
points(y ~ intercept, data = d0)
lines(lower95 ~ intercept, dpred, col = "red")
lines(upper95 ~ intercept, dpred, col = "red")

# Esto funciona!! permite que las evaluaciones se salgan de rango :)

# deberia usar modelos del tipo m4, y parece que son rapidos.
