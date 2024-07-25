
# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(sn)
library(rstan)
library(extraDistr)  # half t density (dht)
library(DHARMa)
library(mnormt) # for pd.solve
library(trialr) # rlkjcorr
library(truncnorm)  # check truncated models
library(bayestestR)      # highest density intervals
library(randtoolbox)   # sobol sequences
library(adaptMCMC)
library(mgcv)

# Data --------------------------------------------------------------------

d <- readRDS(file.path("files", "pseudolikelihood_estimation", "smc_waves_2008_3.rds"))
bounds <- readRDS("files/pseudolikelihood_estimation/bounds_2008_3.rds")
par_names <- colnames(d$par_values)
n_par <- length(par_names)

d$par_values_log <- d$par_values
for(i in 2:ncol(d$par_values)) {
  d$par_values_log[, i] <- log(d$par_values[, i])
}

d <- cbind(d, as.data.frame(d$par_values))

# Functions ---------------------------------------------------------------

# equal-tailed credible interval and mean
etimean <- function(x, ci = 0.95, name = "mu") {
  out <- 1 - ci
  q <- quantile(x, probs = c(out / 2, 1 - out / 2), method = 8)
  result <- c(q[1], mean(x), q[2])
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

# Functions to summaryze
hdi_lower <- function(x, ci = 0.9) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = ci)[["CI_low"]])
}

hdi_upper <- function(x, ci = 0.9) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = ci)[["CI_high"]])
}

hdmean <- function(x, ci = 0.9, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, mean(x), ci$CI_high)
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

# plot data
plot_simulation <- function(data, log = F, alpha = 0.5) {
  x <- data[, "par_values"]
  nn <- colnames(x)
  d <- ncol(x)
  if(log) x <- data[, "par_values_log"]
  par(mfrow = c(ceiling(d / 2), 2))
  for(i in 1:ncol(x)) {
    plot(
      data$overlap ~ x[, i], ylab = "overlap", xlab = nn[i],
      col = rgb(0, 0, 0, alpha), pch = 19
    )
  }
  par(mfrow = c(1, 1))
}

# function to make new data varying only one predictor.
# the mle, if provided, must be named.
# Data is used to take the limits, and it must have a column (matrix) named
# "par_values".
make_newdata <- function(varying = "intercept",
                         data = NULL,
                         mle = NULL, lower = NULL, upper = NULL, ci = 0.98) {

  # limits for prediction
  if(is.null(lower)) {
    if(is.null(ci)) {
      lower <- min(data$par_values[, varying])
    } else {
      lower <- hdi_lower(data$par_values[, varying], ci = ci)
    }
  }

  if(is.null(upper)) {
    if(is.null(ci)) {
      upper <- max(data$par_values[, varying])
    } else {
      upper <- hdi_upper(data$par_values[, varying], ci = ci)
    }
  }

  # list with values for each predictor
  values_list <- lapply(par_names, function(v) {
    if(v != varying) {
      res <- mle[v]
    } else {
      res <- seq(lower, upper, length.out = 150)
    }
    return(res)
  })
  names(values_list) <- par_names

  # grid
  new_data <- expand.grid(values_list)

  # add columns indicating which is the varying predictor and which are its
  # values, useful to plot later
  new_data$varying_var <- varying
  new_data$varying_val <- new_data[, varying]

  return(new_data)
}


# transform parameters vector from and to log scale
params_log <- function(x) {
  x_log <- x
  x_log[2:n_par] <- log(x[2:n_par])
  return(x_log)
}

params_unlog <- function(x_log) {
  x <- x_log
  x[2:n_par] <- exp(x_log[2:n_par])
  return(x)
}

params_unlog_matrix <- function(x_log) {
  x <- x_log
  for(i in 2:n_par) x[, i] <- exp(x_log[, i])
  return(x)
}

# in_bounds: evaluate whether the particles are in_bounds, returning a binary
# vector to fit in_bounds model.
# similarity_bounds argument as provided by the get_bounds function.
in_bounds <- function(data, similarity_bounds,
                      prop = 0.85, var_factor = 2,
                      edge_prop_upper = 0.05,
                      edge_abs_upper = NULL) {

  # # # test
  # data <- particles_data_new
  # similarity_bounds <- bb
  # var_factor = 10
  # prop = 0.85
  # var_factor = 100
  # edge_prop_upper = 0.05
  # # # end testo

  size_diff_max <- similarity_bounds["largest", "size_diff"]
  size_diff_min <- similarity_bounds["smallest", "size_diff"]

  size_diff_lower <- prop * size_diff_min
  size_diff_upper <- prop * size_diff_max

  # get variance in overlap
  low_var <- c(
    similarity_bounds["largest", "var"], # overlap variance, original scale
    similarity_bounds["smallest", "var"]
  ) %>% max

  lowest_var <- low_var * var_factor

  edge_upper <- ifelse(
    is.null(edge_abs_upper),
    edge_prop_upper * similarity_bounds["largest", "edge"],
    edge_abs_upper
  )

  too_large <- ((data$size_diff >= size_diff_upper) &
                  (data$var <= lowest_var)) |
    (data$edge > edge_upper)
  too_small <- (data$size_diff <= size_diff_lower) &
    (data$var <= lowest_var)

  keep <- as.numeric(!(too_large | too_small))

  return(keep)
}

# function to make new data varying only one predictor.
# the mle, if provided, must be named.
# Data is used to take the limits, and it must have a column (matrix) named
# "par_values".
make_newdata <- function(varying = "intercept",
                         data = NULL,
                         mle = NULL, lower = NULL, upper = NULL, ci = 0.98) {

  # limits for prediction
  if(is.null(lower)) {
    if(is.null(ci)) {
      lower <- min(data$par_values[, varying])
    } else {
      lower <- hdi_lower(data$par_values[, varying], ci = ci)
    }
  }

  if(is.null(upper)) {
    if(is.null(ci)) {
      upper <- max(data$par_values[, varying])
    } else {
      upper <- hdi_upper(data$par_values[, varying], ci = ci)
    }
  }

  # list with values for each predictor
  values_list <- lapply(par_names, function(v) {
    if(v != varying) {
      res <- mle[v]
    } else {
      res <- seq(lower, upper, length.out = 150)
    }
    return(res)
  })
  names(values_list) <- par_names

  # grid
  new_data <- expand.grid(values_list)

  # add columns indicating which is the varying predictor and which are its
  # values, useful to plot later
  new_data$varying_var <- varying
  new_data$varying_val <- new_data[, varying]

  return(new_data)
}

# Turn bounded data into zero ----------------------------------------------

# turn data out of bounds into zeroes. Fit a GAM with that dataset,
# and assume overlap zero out of the range of data.

d$inside <- in_bounds(d, bounds)

d$overlap_trans <- d$overlap
d$overlap_trans[d$inside == 0] <- 0.0001

k_side <- 20
bs <- "cr"

m1 <- gam(
  overlap_trans ~
    s(intercept, k = k_side, bs = bs, m = 1) +
    s(vfi, k = k_side, bs = bs, m = 1) +
    s(tfi, k = k_side, bs = bs, m = 1) +
    s(slope, k = k_side, bs = bs, m = 1) +
    s(wind, k = k_side, bs = bs, m = 1) +
    s(steps, k = k_side, bs = bs, m = 1),
  data = d, family = betar(), method = "REML"
)

plot(fitted(m1) ~ d$overlap)
points(fitted(m1)[d$inside == 1] ~ d$overlap[d$inside == 1],
       col = 2, pch = 19)
points(fitted(m1)[d$inside == 0] ~ d$overlap[d$inside == 0],
       col = 1, pch = 19)
mean(d$inside)

summary(m1)

# partial predictions
id_max <- which.max(fitted(m1))
start <- d$par_values[id_max, ]

lowers <- apply(d$par_values, 2, min)
uppers <- apply(d$par_values, 2, max)

fn <- function(x) {
  if(any(x < lowers) | any(x > uppers)) return(0)
  m <- as.data.frame(matrix(x, nrow = 1))
  names(m) <- par_names
  p <- predict(m1, m, type = "response")
  return(p)
}
opt <- optim(start, fn = fn, control = list(fnscale = -1))
opt

new_data <- do.call("rbind", lapply(par_names, function(v) {
  make_newdata(varying = v, data = d, mle = opt$par, ci = NULL)
}))

# evaluate bases at new_data
pp <- predict(m1, new_data, se.fit = TRUE)
pred1 <- cbind(
  new_data,
  mu_mean = plogis(pp$fit),
  mu_lower = plogis(pp$fit - qnorm(0.975) * pp$se.fit),
  mu_upper = plogis(pp$fit + qnorm(0.975) * pp$se.fit)
)

# plot
d_long <- cbind(overlap = d$overlap, as.data.frame(d$par_values))
d_long <- pivot_longer(d_long, all_of(which(names(d_long) %in% par_names)),
                       values_to = "varying_val", names_to = "varying_var")

pred1$varying_var <- factor(pred1$varying_var, levels = par_names)
d_long$varying_var <- factor(d_long$varying_var, levels = par_names)

ggplot(pred1, aes(x = varying_val, y = mu_mean, ymin = mu_lower,
                  ymax = mu_upper)) +
  geom_point(data = d_long, mapping = aes(x = varying_val, y = overlap),
             inherit.aes = F, alpha = 0.05) +
  geom_ribbon(color = NA, alpha = 0.3) +
  geom_line() +
  facet_wrap(vars(varying_var), ncol = 2, scales = "free_x",
             strip.position = "bottom") +
  ylab("overlap") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text.x = element_text(margin = margin(t = 0), vjust = 1),
        axis.title.x = element_blank(),
        panel.spacing.y = unit(5, "mm"))


# fitted values
d$fitted <- fitted(m1)

for(par in par_names) {
  par(mfrow = c(1, 2))
  plot(d$overlap_trans ~ d$par_values[, par], pch = 19, col = rgb(0, 0, 0, 0.1),
       main = "data", xlab = par, ylab = "overlap",
       ylim = range(c(d$overlap, d$fitted)))
  plot(d$fitted ~ d$par_values[, par], pch = 19, col = rgb(0, 0, 0, 0.1),
       main = "predictions", xlab = par, ylab = "overlap",
       ylim = range(c(d$overlap, d$fitted)))
  par(mfrow = c(1, 1))
}


# prediction over square space
# make sobol data to see the effect of knot location
perc_data <- sobol(nrow(d), dim = n_par, seed = 123, init = TRUE)
data_sobol <- sapply(1:n_par, function(i) {
  qunif(perc_data[, i], lowers[i], uppers[i])
})
colnames(data_sobol) <- par_names
data_sobol <- as.data.frame(data_sobol)
data_sobol$fitted <- predict(m1, data_sobol, type = "response")

for(par in par_names) {
  par(mfrow = c(1, 2))
  plot(d$overlap_trans ~ d$par_values[, par], pch = 19, col = rgb(0, 0, 0, 0.1),
       main = "data", xlab = par, ylab = "overlap",
       ylim = range(c(d$overlap, d$fitted)))
  plot(data_sobol$fitted ~ data_sobol[, par], pch = 19, col = rgb(0, 0, 0, 0.1),
       main = "predictions\nin square", xlab = par, ylab = "overlap",
       ylim = range(c(d$overlap, d$fitted)))
  par(mfrow = c(1, 1))
}
# OK, it doesn't take huge values outside the data range.


# posterior distribution.

# function to sample at raw scale, returning -Inf if any parameter
# is out of bounds
log_like <- function(x) {
  if(any(x < lowers) | any(x > uppers)) return(-Inf)
  # evaluate the bases at x
  XX <- as.data.frame(matrix(x, nrow = 1))
  names(XX) <- par_names
  ld <- predict(m1, XX, type = "response") %>% log
  return(ld)
}

# MCMC at log scale
sampling_iters <- 20000
adapt_iters <- 2000
mcmc1 <- MCMC(log_like,
              n = sampling_iters + adapt_iters,
              adapt = adapt_iters,
              scale = as.numeric(apply(d$par_values, 2, sd) / 2),
              init = start, acc.rate = 0.234)

draws1 <- mcmc1$samples
colnames(draws1) <- par_names

par(mfrow = c(3, 2))
for(par in par_names) {

  if(par == "intercept") {
    den <- density(draws1[-(1:adapt_iters), par])
  } else {
    den <- density(draws1[-(1:adapt_iters), par], from = 0)
  }

  den$y <- den$y * (max(d$overlap) / max(den$y))
  plot(d$overlap ~ d$par_values[, par], main = par,
       col = rgb(0, 0, 0, 0.05))
  lines(den, main = par, lwd = 2, col = "red")

}
par(mfrow = c(1, 1))


# Am I sampling wrong? MCMC problems?
# resample a discrete space according to likelihood.
large_grid_unit <- sobol(1e6, dim = n_par, seed = 123, init = TRUE)
large_grid <- sapply(1:n_par, function(i) {
  qunif(large_grid_unit[, i], lowers[i], uppers[i])
})
colnames(large_grid) <- par_names
large_grid <- as.data.frame(large_grid)
large_grid$fitted <- predict(m1, large_grid, type = "response")

# resample
sample_ids <- sample(1:nrow(large_grid), size = 1e5, prob = large_grid$fitted,
                     replace = T)


par(mfrow = c(3, 2))
for(par in par_names) {

  if(par == "intercept") {
    den <- density(large_grid[sample_ids, par])
  } else {
    den <- density(large_grid[sample_ids, par], from = 0)
  }

  den$y <- den$y * (max(d$overlap) / max(den$y))
  plot(d$overlap ~ d$par_values[, par], main = par,
       col = rgb(0, 0, 0, 0.05))
  lines(den, main = par, lwd = 2, col = "red")

}
par(mfrow = c(1, 1))



### IDEA:
# to avoid the problem of space-filling,
# fit a convoluted GP but only considering pair-wise interactions. this way
# it would be like a spline, but with constrained shape.

# pair-wise interactions:
(n_interact <- ncol(combn(1:n_par, 2)))
side_dim <- 3
side_dim ^ 2
(n_interact * side_dim ^ 2)
n_par * 10 # 60 marginal basis

# 6 marginal sigmas
# 1 interaction sigma (for all pairs)

# se podría ajustar con lm() y ya? sin restricciones?
# raro.

