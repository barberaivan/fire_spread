# Code to test how to fit a GP to data, and how to sample from it.
# To constrain the sampling space, a density if fitted to data above a threshold,
# and where the density is low, the GP is not evaluated, the log-like is set to
# zero.

# Fit to all data points, or only to good ones?
# overlap scale, or log?

# Packages ----------------------------------------------------------------

library(tidyverse); theme_set(theme_bw())
library(viridis)
library(GauPro)
library(DHARMa)
library(bayestestR)    # highest density intervals
library(randtoolbox)   # sobol sequences
library(kdevine)       # fit multivariate distributions to data
library(adaptMCMC)
library(posterior)     # to manage samples

# Data --------------------------------------------------------------------

d <- readRDS(file.path("files", "pseudolikelihood_estimation",
                       "smc_waves_2008_3.rds"))
similarity_bounds <- readRDS(file.path("files", "pseudolikelihood_estimation",
                                       "bounds_2008_3.rds"))

par_names <- colnames(d$par_values)
n_par <- length(par_names)

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

# find MLE in the fitted GP
gp_optimize <- function(model, data) {

  fitted_ov <- model$pred(data$par_values, se.fit = F)
  id_max <- which.max(fitted_ov)
  start <- data$par_values[id_max, ]

  # bounds
  lowers <- apply(data$par_values, 2, min)
  uppers <- apply(data$par_values, 2, max)

  op <- optim(start, fn = function(x) model$pred(matrix(x, nrow = 1), se.fit = F),
              control = list(fnscale = -1, maxit = 1e5),
              method = "L-BFGS-B", lower = lowers, upper = uppers)

  return(op)
}

# function to make new data varying only one predictor.
# the mle, if provided, must be named.
# Data is used to take the limits, and it must have a column (matrix) named
# "par_values".
make_newdata <- function(varying = "intercept",
                         data = NULL, mle = NULL, ci = 0.98) {

  # limits for prediction
  if(is.null(ci)) {
    lower <- min(data$par_values[, varying])
    upper <- max(data$par_values[, varying])
  } else {
    lower <- hdi_lower(data$par_values[, varying], ci = ci)
    upper <- hdi_upper(data$par_values[, varying], ci = ci)
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

# partial predictor function
partial_predictions <- function(data, model, mle, ci = NULL) {

  ## TEST
  # varying = "all"
  # data = data_pred
  # model = model; mle = mle
  ###

  new_data <- do.call("rbind", lapply(par_names, function(v) {
    make_newdata(varying = v, data = d, mle = mle, ci = ci)
  }))

  pred <- model$pred(new_data[, par_names], se.fit = TRUE)

  new_data$mle <- pred$mean
  new_data$upper <- pred$mean + qnorm(0.975) * pred$se
  new_data$lower <- pred$mean - qnorm(0.975) * pred$se

  return(new_data)
}

# Function to plot the GP.
gp_partials <- function(model, data, ci = NULL) { # parameter to vary in the 1d plot

  data$in_bounds <- factor(as.character(data$in_bounds), levels = c("0", "1"))

  # compute partial predictions
  op <- gp_optimize(model, data)
  preds <- partial_predictions(data = data, model = model,
                               mle = op$par, ci = ci)

  # longanize data to plot
  data_ <- cbind(overlap = data$overlap, in_bounds = data$in_bounds,
                 as.data.frame(data$par_values))
  data_long <- pivot_longer(data_, all_of(which(names(data_) %in% par_names)),
                            values_to = "varying_val", names_to = "varying_var")

  data_long$varying_var <- factor(data_long$varying_var, levels = par_names)
  preds$varying_var <- factor(preds$varying_var, levels = par_names)

  p <-
    ggplot() +

    # data
    geom_point(data = data_long,
               mapping = aes(x = varying_val, y = overlap,
                             color = in_bounds, shape = in_bounds),
               size = 2, alpha = 0.5) +
    scale_shape_manual(values = c(16, 17)) +
    scale_color_viridis(end = 0.7, discrete = TRUE) +

    # model predictions
    ggnewscale::new_scale_color() +
    scale_color_viridis(end = 0.7, discrete = TRUE, option = "A") +
    ggnewscale::new_scale_fill() +
    scale_fill_viridis(end = 0.7, discrete = TRUE, option = "A") +

    geom_ribbon(data = preds, mapping = aes(x = varying_val, y = mle,
                                            ymin = lower, ymax = upper),
                alpha = 0.3, color = NA) +

    geom_line(data = preds, mapping = aes(x = varying_val, y = mle)) +

    facet_wrap(vars(varying_var), scales = "free_x", ncol = 2,
               strip.position = "bottom") +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right",
          strip.placement = "outside",
          strip.background = element_rect(fill = "white", color = "white"),
          strip.text.x = element_text(margin = margin(t = 0), vjust = 1),
          axis.title.x = element_blank(),
          panel.spacing.y = unit(5, "mm")) +
    ylab("overlap")

  print(p)
  return(list(plot = p, data_points = data_long, data_preds = preds))
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

# gp fit ------------------------------------------------------------------

d$in_bounds <- in_bounds(d, similarity_bounds)

# gp <- gpkm(
#   X = d$par_values,
#   Z = d$overlap,
#   parallel = F, useC = TRUE, nug.max = 100,
#   kernel = "matern52"
# ) # 5 min para 1800 obs. No es tanto.
# saveRDS(
#   gp,
#   file.path("files", "pseudolikelihood_estimation", "gp-2008_3.rds")
# )
gp <- readRDS(file.path("files", "pseudolikelihood_estimation", "gp-2008_3.rds"))
p1 <- gp_partials(gp, d)
p1$plot

gp_high <- gpkm(
  X = d$par_values[d$overlap >= quantile(d$overlap, 0.65), ],
  Z = d$overlap[d$overlap >= quantile(d$overlap, 0.65)],
  parallel = F, useC = TRUE, nug.max = 100,
  kernel = "matern52"
)# 19 s para 630 obs.
p2 <- gp_partials(gp_high, d)
p2$plot # ajusta re mal!!


# check fitted values
gp_preds <- gp$pred(d$par_values, se.fit = T)
d$overlap_fitted <- gp_preds$mean

plot(overlap_fitted ~ overlap, data = d, pch = 19, col = rgb(0, 0, 0, 0.1))
abline(0, 1, col = "red")

# DHARMa
y_sim <- sapply(1:2000, function(i) {
  rnorm(nrow(d), gp_preds$mean, gp_preds$se)
})

res <- createDHARMa(y_sim, d$overlap, fittedPredictedResponse = gp_preds$mean)
plot(res) # muy mal, pero es esperable, por el supuesto de normalidad.
for(p in par_names) {
  plotResiduals(res, form = d$par_values[, p], rank = F, xlab = p,
                quantreg = F)
}


# fit density to data > threshold.
threshold <- 0.1
d$overlap_high <- d$overlap
d$overlap_high[d$overlap < 0.1] <- 0

nsim <- 5000
rows_boot <- sample(1:nrow(d), nsim, prob = d$overlap, replace = T)

den_fun <- kdevine(d$par_values[rows_boot, ], xmin = c(-Inf, rep(0, n_par - 1)))
den_fun <- kdevine(d$par_values, xmin = c(-Inf, rep(0, n_par - 1)))

# evaluate density estimate
fitted_density <- dkdevine(d$par_values, den_fun)
plot(fitted_density ~ d$overlap)
plot(fitted_density)


draw <- rkdevine(3000, den_fun) # re tarda en muestrearlo.

for(i in 1:n_par) {
  plot(density(d$par_values[, i]), main = par_names[i])
  lines(density(draw[, i]), col = "red")
}

# the functions vine() and vinecop() from the 'rvinecopulib'
# package as replacements for kdevine() and kdevinecop().


# MCMC trials -------------------------------------------------------------

# can we sample from the GP?
lowers <- apply(d$par_values, 2, min)
uppers <- apply(d$par_values, 2, max)
threshold <- 0.1

# function to sample at raw scale, returning -Inf if any parameter
# is out of bounds
log_like <- function(x) {
  if(any(x < lowers) | any(x > uppers)) return(-Inf)
  like <- gp$pred(matrix(x, nrow = 1), se.fit = F)
  if(like < threshold) return(-Inf)
  return(log(like))
}

# MCMC at log scale
sampling_iters <- 50000
adapt_iters <- 5000

make_init <- function() {
  row <- sample(which(d$overlap > median(d$overlap)), 1)
  return(d$par_values[row, ])
}

nc <- 8
mcmc1 <- MCMC.parallel(log_like,
                       n = sampling_iters + adapt_iters,
                       adapt = adapt_iters,
                       n.cpu = nc,
                       n.chain = nc,
                       scale = as.numeric(apply(d$par_values, 2, sd) / 2),
                       init = make_init(), acc.rate = 0.234,
                       packages = c("GauPro"))

# extract draws removing adaptation
draws_list <- lapply(1:nc, function(i) {
  draws <- mcmc1[[i]]$samples[-(1:adapt_iters), ]
  draws_thinned <- draws[seq(1, nrow(draws), by = 10), ]
  return(mcmc(draws_thinned))
})

draws_mcmc <- mcmc.list(draws_list) # from coda
dl <- as_draws_list(draws_mcmc)     # from posterior
sm2 <- posterior::summarise_draws(dl)
print(sm2)

draws_dfdf <- as_draws_df(dl)
plot(dl)
plot(draws_mcmc)
pairs(draws_dfdf[seq(1, nrow(draws_dfdf), by = 5), par_names],
      col = rgb(0, 0, 0, 0.05))

str(draws_dfdf)

# plot data and marginal densities
par(mfrow = c(3, 2))
for(p in par_names) {
  # p = "slope"
  if(p != "intercept") {
    ddd <- density(draws_dfdf[, p, drop = T], from = 0, adjust = 0.75)
  } else {
    ddd <- density(draws_dfdf[, p, drop = T], adjust = 0.75)
  }
  plot(d$overlap ~ d$par_values[, p], ylab = "overlap", xlab = p,
       main = NULL, col = rgb(0, 0, 0, 0.2), pch = 19)
  ddd$y <- (ddd$y / max(ddd$y)) * max(d$overlap)
  lines(ddd, col = "red", lty = 1.5)
}
par(mfrow = c(1, 1))


# plot data and fitted overlap (GP)
par(mfrow = c(3, 2))
for(p in par_names) {
  plot(d$overlap ~ d$par_values[, p], ylab = "overlap", xlab = p,
       main = NULL, col = rgb(0, 0, 0, 0.2), pch = 19)
  points(d$overlap_fitted ~ d$par_values[, p], col = "red")
}
par(mfrow = c(1, 1))




# Cambia la posterior si simulamos del GP en el muestreo? -----------------

log_like_sim <- function(x) {
  if(any(x < lowers) | any(x > uppers)) return(-Inf)
  pred <- gp$pred(matrix(x, nrow = 1), se.fit = T)
  like <- rnorm(1, pred$mean, pred$se)
  if(like < threshold) return(-Inf)
  return(log(like))
}

mcmc2 <- MCMC.parallel(log_like_sim,
                       n = sampling_iters + adapt_iters,
                       adapt = adapt_iters,
                       n.cpu = nc,
                       n.chain = nc,
                       scale = as.numeric(apply(d$par_values, 2, sd) / 2),
                       init = make_init(), acc.rate = 0.234,
                       packages = c("GauPro"))

# extract draws removing adaptation
draws_list2 <- lapply(1:nc, function(i) {
  draws <- mcmc2[[i]]$samples[-(1:adapt_iters), ]
  draws_thinned <- draws[seq(1, nrow(draws), by = 10), ]
  return(mcmc(draws_thinned))
})

draws_mcmc2 <- mcmc.list(draws_list2) # from coda
dl2 <- as_draws_list(draws_mcmc2)     # from posterior
sm2 <- posterior::summarise_draws(dl2)
# le remil cuesta muestrear así

draws_dfdf2 <- as_draws_df(dl2)


# compute quantiles
pp <- ppoints(100)
q_det <- apply(as.matrix(draws_dfdf[, par_names]), 2, quantile, probs = pp)
q_sim <- apply(as.matrix(draws_dfdf2[, par_names]), 2, quantile, probs = pp)

for(p in par_names) {
  plot(q_det[, p] ~ q_sim[, p], xlab = "stochastic", ylab = "deterministic",
       pch = 19, main = p)
  abline(0, 1, col = "red", lty = 2)
}
# no es comparable, porque muestreó para la mierda.


# probarlo unidimensionalmente, fijando los demás params en la media.



# MCMC accepting only particles above a threshold -------------------------
# like the classica ABC, but only allowing high overlap.

# evaluate high observed and fitted overlap to define a high threshold
quantile(d$overlap, probs = seq(0.9, 1, by = 0.02))
quantile(d$overlap_fitted, probs = seq(0.9, 1, by = 0.02))
range(gp_preds$se[d$overlap_fitted > 0.5])

# 3 sigmas below from the max_fitted

threshold_high <- gp_preds$mean[which.max(gp_preds$mean)] - 10 * gp_preds$se[which.max(gp_preds$mean)]
# near the observed 98 % of data.

sum(gp_preds$mean > threshold_high) # only 26 observations above thres.

# function to sample at raw scale, returning -Inf if any parameter
# is out of bounds
overlap_80 <- quantile(d$overlap, prob = 0.8)
lowers_80 <- apply(d$par_values[d$overlap > overlap_80, ], 2, min)
uppers_80 <- apply(d$par_values[d$overlap > overlap_80, ], 2, max)

log_like_high <- function(x) {
  if(any(x < lowers_80) | any(x > uppers_80)) return(-Inf)
  like <- gp$pred(matrix(x, nrow = 1), se.fit = F)
  res <- ifelse(like < threshold_high, -Inf, 0)
  return(res)
}

# MCMC at log scale
sampling_iters <- 50000
adapt_iters <- 5000

make_init_high <- function() {
  row <- sample(which(d$overlap > threshold_high), 1)
  return(d$par_values[row, ])
}

nc <- 8
mcmc3 <- MCMC.parallel(log_like_high,
                       n = sampling_iters + adapt_iters,
                       adapt = adapt_iters,
                       n.cpu = nc,
                       n.chain = nc,
                       scale = as.numeric(apply(d$par_values[d$overlap >= threshold_high, ], 2, sd) / 2),
                       init = make_init_high(), acc.rate = 0.234,
                       packages = c("GauPro"))

# extract draws removing adaptation
draws_list <- lapply(1:nc, function(i) {
  draws <- mcmc3[[i]]$samples[-(1:adapt_iters), ]
  draws_thinned <- draws[seq(1, nrow(draws), by = 10), ]
  return(mcmc(draws_thinned))
})

draws_mcmc3 <- mcmc.list(draws_list) # from coda
dl3 <- as_draws_list(draws_mcmc3)     # from posterior
sm3 <- posterior::summarise_draws(dl3)
print(sm3)

draws_dfdf3 <- as_draws_df(dl3)
# plot(draws_mcmc3)

# compare overlap sampling and high quality sampling
par(mfrow = c(3, 2))
for(par in par_names) {
  # densities
  if(par == "intercept") {
    d1 <- density(draws_dfdf[, par, drop = T])
    d3 <- density(draws_dfdf3[, par, drop = T])
  } else {
    d1 <- density(draws_dfdf[, par, drop = T], from = 0)
    d3 <- density(draws_dfdf3[, par, drop = T], from = 0)
  }

  xx <- range(c(d1$x, d3$x))
  yy <- range(c(d1$y, d3$y))

  plot(d1, type = "l", main = NA, xlab = par, ylab = "density",
       xlim = xx, ylim = yy)
  lines(d3, col = "red")
}
par(mfrow = c(1, 1))



# Get posterior density with kdevine --------------------------------------

draws_dfdf3 %>% nrow()
n_few <- 4000
ids_few <- seq(1, nrow(draws_dfdf3), by = nrow(draws_dfdf3) / n_few)
draws_dfdf3_few <- draws_dfdf3[ids_few, par_names]
pairs(draws_dfdf3_few, col = rgb(0, 0, 0, 0.1), pch = 19)

den_fun3 <- kdevine(as.matrix(draws_dfdf3_few),
                    xmin = c(-Inf, rep(0, n_par - 1)))

fitted_density <- dkdevine(d$par_values, den_fun3)
plot(fitted_density ~ d$overlap)
plot(log(fitted_density) ~ d$overlap)
# llega a -Inf llegado el punto,
dlog <- log(fitted_density)
min(dlog[is.finite(dlog)])
sum(is.infinite(dlog)) / nrow(d) # 57 % es -Inf
sum(fitted_density == 0) / nrow(d) # 57 % es 0

apply(d$par_values[fitted_density > 0, ], 2, range)
apply(d$par_values[d$overlap_fitted > threshold_high, ], 2, range)

# La density suaviza un poquito, pero la extensión de rango que hace es muy
# pequeña. O sea que en la práctica vamos a tener soporte compacto.


# Y una MVN? tenderá a cero?
library(MGMM)          # fit MVN distribution

log_params <- function(x) {
  xlog <- x
  for(i in 2:n_par) xlog[, i] <- log(x[, i])
  return(xlog)
}

exp_params <- function(x) {
  xexp <- x
  for(i in 2:n_par) xexp[, i] <- exp(x[, i])
  return(xexp)
}

d$par_values_log <- log_params(d$par_values)
draws_few_log <- as.matrix(log_params(draws_dfdf3_few))
mvn_fitted <- FitGMM(draws_few_log)

fitted_den_mvn <- apply(
  d$par_values_log, 1,
  FUN = function(x) mgcv::dmvn(x, mu = mvn_fitted@Mean,
                               V = mvn_fitted@Covariance * 0.25 ^ 2)
)
range(fitted_den_mvn)
# Ok, it doesn't underflow

plot(fitted_den_mvn ~ d$overlap, col = rgb(1, 0, 0, 0.1), pch = 19)
points(dlog ~ d$overlap, col = rgb(0, 0, 0, 0.8), pch = 19)
# It works. So, whenever we find a zero density, we use the normal MVN

plot(fitted_den_mvn[is.finite(dlog)] ~ dlog[is.finite(dlog)],
     col = rgb(1, 0, 0, 0.1), pch = 19)
abline(0, 1)
# Las MVN son más altas. mejor encogerle la varianza.
# con un shrink factor de 0.25 ^ 2 anda bien.


# Mejor aún: para disminuir el nro de aproximaciones, usar el GP para evaluar
# la likelihood. Y cuando el GP tire cero, ahí usar la MVN. Pero acá la MVN
# debería estar escalada para que el salto no sea tan bruto. O sea, tien que
# tener la misma densidad máxima que el GP.
# Y el GP debería tener como piso el cero.


# max density mvn:
mgcv::dmvn(mvn_fitted@Mean, mu = mvn_fitted@Mean,
           V = mvn_fitted@Covariance * 0.25 ^ 2)
# esto debería ser seteado en
log(max(op))
gp_op <- gp_optimize(gp, d)
log(gp_op$value)

# igual no estoy seguro...
# la mínima logden del GP será
log(threshold_high)
# luego de eso habrá un salto re grande hacia la mvn density.
# quizás sólo sea cuestión de escalar bien esa MVN, sin "desnormalizarla"..

# Ajustar una skew-MVN a la puntita? --------------------------------------

# primero hacerle MCMC al GP y luego estimar skew-mvn.

library(sn)

f2 <- sn::selm.fit(
  x = matrix(rep(1, nrow(draws_few_log)), ncol = 1),
  y = draws_few_log,
  family = "SN"
)
f2$param$dp
d2 <- sn::makeSECdistr(dp = f2$param$dp, family = "SN", compNames = par_names)
plot(d2)
r2 <- sn::rmsn(n = 4e4, dp = f2$param$dp)
r2exp <- exp_params(r2)

# y comparar también con las muestras del kdevine
r3 <- rkdevine(2000, den_fun3) # re tarda en muestrearlo.
r3log <- log_params(r3)

par(mfrow = c(3, 2))
for(i in 1:n_par) {

  d1 <- density(draws_few_log[, i]) # posterior samples
  d2 <- density(r2[, i])  # skew mvn
  d3 <- density(r3log[, i]) # kdevine

  xx <- range(d1$x, d2$x, d3$x)
  yy <- range(d1$y, d2$y, d3$y)

  plot(d1, xlab = par_names[i], ylim = yy, xlim = xx, main = NA)
  lines(d2, col = "red")
  lines(d3, col = "blue")
}
par(mfrow = c(1, 1))

# El MVN anda muy mal, y el kdevine, muy bien.


# Likelihoods de una poisson con datos muy serparados. se hace compacto? ----

a1 <- 1
a2 <- 1e10
curve(dpois(a1, x), to = a2 * 1.3, n = 500)
curve(dpois(a2, x), add = TRUE, col = "red", n = 500)
curve(dpois(a1, x) * dpois(a2, x), to = a2 * 1.3, n = 500)

# en log
curve(dpois(a1, x, log = T), to = a2 * 1.3, n = 500)
curve(dpois(a2, x, log = T), add = T, col = "red", n = 500)
curve(dpois(a1, x, log = T) + dpois(a2, x, log = T), to = a2 * 1.3, n = 500)

# Nunca llega a hacer underflow en escala log, aunque en escala exp, sí.