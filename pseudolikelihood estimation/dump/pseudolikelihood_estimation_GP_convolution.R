# Try to fit a GP by convolution, at the logit scale. The intercept will be
# fixed at zero, so when there is no data the GP tends to zero. Fit at the
# overlap scale.

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

# Data --------------------------------------------------------------------

d <- readRDS(file.path("files", "pseudolikelihood_estimation", "smc_waves_2008_3.rds"))
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

# unnormalized log-density of a multivariate normal distribution, computed from
# the location vector xi and vcov matrix Omega.
multi_normal_lupdf <- function(X, xi, Omega, log = T) {
  N <- nrow(X); K <- ncol(X)
  log_den <- numeric(N)
  P <- pd.solve(Omega, silent = TRUE, log.det = F) # same computation as in {sn}

  for(n in 1:N) {
    x_centred = X[n, ] - xi
    log_den[n] = -0.5 * t(x_centred) %*% P %*% x_centred;
  }

  if(log) return (log_den) else return(exp(log_den))
}

# unnormalized log-density of a multivariate skew-normal distribution, defined by
# location vector xi, precision matrix P, and slant vector alpha.
multi_skew_normal_lupdf <- function(X, xi, Omega, alpha) {

  N <- nrow(X); K <- ncol(X)
  log_den_skew <- numeric(N)
  sigma <- sqrt(diag(Omega))
  P <- pd.solve(Omega, silent = TRUE, log.det = F) # same computation as in {sn}
  alpha <- matrix(alpha, nrow = 1)
  shifted <- numeric(N)

  for(n in 1:N) {
    x_centred <- X[n, ] - xi
    # normal log-density
    log_den <- -0.5 * t(x_centred) %*% P %*% x_centred;
    # normal cumulative probability for shifted value
    shifted[n] <- alpha %*% (x_centred / sigma)
    log_cum <- pnorm(shifted[n], log.p = T)
    # skew normal log-density
    log_den_skew[n] <- log(2) + log_den + log_cum
  }

  return (list(mu = log_den_skew, shifted = shifted))
}


multi_skew_normal_debug <- function(X, xi, Omega, alpha) {
  N <- nrow(X); K <- ncol(X)
  log_den_skew <- numeric(N)
  sigma <- sqrt(diag(Omega))
  P <- pd.solve(Omega, silent = TRUE, log.det = F) # same computation as in {sn}
  alpha <- matrix(alpha, nrow = 1)

  log_den <- numeric(N); log_cum <- numeric(N)
  shifted <- numeric(N)

  for(n in 1:N) {
    x_centred <- X[n, ] - xi
    # normal log-density
    log_den[n] <- -0.5 * t(x_centred) %*% P %*% x_centred;
    # normal cumulative probability for shifted value
    shifted[n] <- alpha %*% (x_centred / sigma)
    log_cum[n] <- pnorm(shifted[n], log.p = T)
    # skew normal log-density
    log_den_skew[n] <- log(2) + log_den[n] + log_cum[n]
  }

  res <- cbind("log_den" = log_den, "log_cum" = log_cum, "shifted" = shifted,
               "log_sn" = log_den_skew)

  return (res)
}


# make correlation matrix from unconstrainde partial correlation parameters,
# matrix dimension and bound on partial correlation coefficients. Although
# tanh transforms from (-Inf, Inf) to (-1, 1), for extreme values the
# resulting matrix might be non-positive-definite.
# This is how Stan defines correlation matrices to meet all constraints.
# y is the vector of unconstrained (-Inf, Inf) partial correlation parameters,
# and d is te matrix dimension. cor_bound is the maximum absolute partial
# correlation allowed.
make_corr <- function(y, d, cor_bound = 0.95) {
  # # TEST
  # d <- 4
  # y <- rnorm(cor_num(d), sd = 1e5)
  # cor_bound <- 0.95
  # # END TEST

  # upper triangular matrix of Canonical Partial Correlations
  z <- matrix(0, d, d)
  z[upper.tri(z)] <- tanh(y) * cor_bound

  # compute omega (upper choleski factor for corr matrix)
  w <- matrix(0, d, d)

  w[1, ] <- c(1, z[1, 2:d])
  for(i in 2:d) {
    for(j in i:d) {
      if(i == j) w[i, j] <- prod(sqrt(1 - z[1:(i-1), j] ^ 2))
      if(i < j) w[i, j] <- z[i, j] * prod(sqrt(1 - z[1:(i-1), j] ^ 2))
    }
  }

  Rho <- t(w) %*% w

  ## Check
  # is_positive_definite(Rho); isSymmetric(Rho)
  # range(Rho[lower.tri(Rho)])
  ## End check

  return(Rho)
}

# function to compute the number of correlation parameters as a function of the
# multivariate distribution dimension
cor_num <- function(d) ncol(combn(d, 2))

# inverse function: detect dimension of the matrix from the length of the 1D
# vector of correlation coefficients
cor_dim <- function(l) 1 / 2 * (sqrt(8 * l + 1) + 1)
# using wolframalpha: https://www.wolframalpha.com/widgets/view.jsp?id=c86d8aea1b6e9c6a9503a2cecea55b13


# function to get unconstrained xi (position) from the real xi
# (used to provide starting values from a selm fit)
unconstrain_xi <- function(xi, xi_lower, xi_upper, prior_xi_logit_sd) {
  if(any(xi < xi_lower) | any(xi > xi_upper)) stop("xi out of range")
  p <- (xi - xi_lower) / (xi_upper - xi_lower)
  xi_raw <- qlogis(p) / prior_xi_logit_sd
  return(xi_raw)
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


normalize <- function(x) x / sum(x)


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

# Prepare data ------------------------------------------------------------

# work at log scale, because steps has a heavily skew distribution,
# and a log-normal would be better.

d$par_values_log <- d$par_values
for(i in 2:ncol(d$par_values)) {
  d$par_values_log[, i] <- log(d$par_values[, i])
}

plot_simulation(d)

plot(d$overlap ~ d$par_values[, "intercept"])
low <- 0.2 # threshold
y_lower <- log(low)
data <- d[d$overlap_log > y_lower, ]
nrow(data)

# reasonable values for sigma
x_ranges <- sapply(1:ncol(data$par_values), function(c) {
  diff(range(data$par_values_log[, c]))
})

sds <- apply(data$par_values_log, 2, sd)

# data out
X = data$par_values_log
y = data$overlap_log
N = nrow(data)
K = ncol(X)

# parameters will be heavily constrained for the surface to be fitted.
sdata <- list(
  # data
  X = X, y = y, N = N, K = K,

  # lower truncation for the response
  L = y_lower,

  # force b0 to be near the optimum
  prior_b0_mean = max(y),
  prior_b0_sd = diff(range(y)), # sd(y) * 4

  # constrain xi (~mu) within the 90 % HDI of the observations
  xi_lower = apply(X, 2, hdi_lower),
  xi_upper = apply(X, 2, hdi_upper),
  prior_xi_logit_sd = 1.5,

  # prior for sigma (student_t)
  prior_sigma_sd = x_ranges * 2,
  prior_sigma_nu = 3,

  # prior for slant (alpha) (half-normal)
  prior_alpha_sd = 3,

  # residual standard deviation (student_t)
  prior_tau_sd = sd(y) * 2,
  prior_tau_nu = 3,

  # prior eta parameter for correlation matrix
  prior_corr_eta = 1
)

# compile stan model
smodel <- stan_model(
  file.path("pseudolikelihood estimation",
            "pseudolikelihood_estimation_smvn_01.stan"), # truncated, skew, mvn
  verbose = F
)



# Simulations: GP 1D ------------------------------------------------------

# make real function (2 gaussian bases)
gaussian_base <- function(x = seq(-3, 3, length.out = 100),
                          top = 1, sd = 1, mu = 0) {
  exp(-((x - mu) / sd) ^ 2) * top
}
curve(gaussian_base(x), from = -3, to = 3)

# function to create real data
my_fun <- function(x) {
  gaussian_base(x) + gaussian_base(x, top = 1.5, mu = 2.5, sd = 1)
}
curve(my_fun, from = -5, to = 5)

# synthetic data:
N <- 500
x <- runif(N, -2, 4.5)
tau_noise <- 0.05
mu <- my_fun(x)
y <- rnorm(N, mu, tau_noise)
points(y ~ x)

# Bases:
k <- 5
sigma <- (diff(range(x)) / k) * 2
xk <- seq(min(x), max(x), length.out = k) # k for knots

deltas_norm <- outer(x, xk, "-") / sigma
w <- exp(-deltas_norm ^ 2)
range(w)
w_norm <- apply(w, 1, normalize) %>% t

# weighted average at knots?
# here we use a small sigma, to make the avg very local
deltas_norm_knots <- outer(x, xk, "-") / (sigma / 3)
w_knots <- exp(-deltas_norm_knots ^ 2)
w_norm_knots <- apply(w_knots, 2, normalize) %>% t
str(w_norm_knots)
knots_mean <- as.numeric(w_norm_knots %*% y)

curve(my_fun, from = -20, to = 20, lwd = 2, ylim = c(-10, 10))
abline(v = xk, lty = 2)
points(y ~ x, pch = 19, col = rgb(0, 0, 0, 0.1))
points(knots_mean ~ xk, col = "red", pch = 19)

# linear model to unnormalized weights
m1 <- lm(y ~ w - 1)
m1 <- lm(y ~ w_norm - 1)
xseq <- seq(min(x) - 10, max(x) + 6, length.out = 200)
deltas_norm_pred <- outer(xseq, xk, "-") / sigma
w_pred <- exp(-deltas_norm_pred ^ 2)
w_pred <- apply(w_pred, 1, normalize) %>% t
mu_pred <- w_pred %*% coef(m1)
lines(mu_pred ~ xseq, col = "blue", lwd = 2)

coef(m1); knots_mean

# If all coeffs are positive, is the smooth positive?

k <- 15
sigma <- (diff(range(x)) / k) * 1
xk <- seq(min(x), max(x), length.out = k) # k for knots

xseq <- seq(min(x) - 5, max(x) + 5, length.out = 200)
deltas_norm_pred <- outer(xseq, xk, "-") / sigma
w_pred <- exp(-deltas_norm_pred ^ 2)
# w_pred <- apply(w_pred, 1, normalize) %>% t
mu_pred <- w_pred %*% abs(rnorm(k, sd = 0.01))#coef(m1)
plot(mu_pred ~ xseq, col = "blue", lwd = 2, type = "l")
# abline(v = xk, lty = 2)

# yes, positive coeffs make a positive function.



# Try lm with real data ---------------------------------------------------

plot_simulation(data)

n_par <- ncol(data$par_values)
par_names <- colnames(data$par_values)

lowers <- apply(data$par_values, 2, hdi_lower, ci = 0.98)
uppers <- apply(data$par_values, 2, hdi_upper, ci = 0.98)
widths <- uppers - lowers

n_knots <- 50
perc <- sobol(n_knots, dim = n_par, seed = 123, init = TRUE)
knots <- sapply(1:n_par, function(i) {
  qunif(perc[, i], lowers[i], uppers[i])
})
colnames(knots) <- par_names

# define kernel
# sigma is a fraction of the distance between knots
sigmas <- (widths / 2)
V <- diag(sigmas ^ 2)

# weights for all knots
weights_un <- sapply(1:n_knots, function(k) {
  multi_normal_lupdf(data$par_values, knots[k, ], V, log = F)
})
str(weights_un)
hist(colSums(weights_un), breaks = 10)
hist(rowSums(weights_un), breaks = 10)

m1 <- lm(data$overlap ~ weights_un - 1)
sm1 <- summary(m1)
sm1$r.squared
hist(coef(m1))

plot(data$overlap ~ fitted(m1)); abline(0, 1, col = "red")
cc <- cov2cor(vcov(m1))
hist(cc[lower.tri(cc)])
cc[lower.tri(cc)] %>% range # no high correlation


# with large sigmas the model fits perfectly, but this might make it hard
# to go towards zero if the particles are out of range.

# predict over the whole data set to check:
weights_un_pred <- sapply(1:n_knots, function(k) {
  multi_normal_lupdf(d$par_values, knots[k, ], V, log = F)
})

predd <- weights_un_pred %*% coef(m1)
plot(predd ~ d$overlap); abline(0, 1, col = "red")


# bayesian model, prior check  ---------------------------------------------

# fit it with stan so we can constrain the sign and shrink the coefficients.

# separate data in two groups
plot(d$overlap ~ d$par_values[, "intercept"])
abline(h = 0.1, col = 2, lwd = 2)

threshold <- 0.2
data_in <- d[d$overlap > threshold, ]
data_out <- d[d$overlap <= threshold, ]

# prior check.
lowers <- apply(data_in$par_values, 2, hdi_lower, ci = 0.98)
uppers <- apply(data_in$par_values, 2, hdi_upper, ci = 0.98)

minims <- apply(d$par_values, 2, min)
maxims <- apply(d$par_values, 2, max)


# define kernel
# sigma is a fraction of the distance between knots
sigmas <- (widths / 2)
V <- diag(sigmas ^ 2)
# knots
n_knots <- 400

# (prueba)
# lowers <- apply(data$par_values, 2, min)
# uppers <- apply(data$par_values, 2, max)

perc <- sobol(n_knots, dim = n_par, seed = 123, init = TRUE)
knots <- sapply(1:n_par, function(i) {
  qunif(perc[, i], lowers[i], uppers[i])
})
colnames(knots) <- par_names

# make sobol data to see the effect of knot location
perc_data <- sobol(nrow(d), dim = n_par, seed = 123, init = TRUE)
data_sobol <- sapply(1:n_par, function(i) {
  qunif(perc_data[, i], minims[i], maxims[i])
})
colnames(data_sobol) <- par_names

weights_sobol <- sapply(1:n_knots, function(k) {
  multi_normal_lupdf(data_sobol, knots[k, ], V, log = F)
})

# weights for all data points (includes out of range)
d_out <- d[d$overlap <= 0.2, ]
weights_pred_out <- sapply(1:n_knots, function(k) {
  multi_normal_lupdf(data_out$par_values, knots[k, ], V, log = F)
})
weights_pred_in <- sapply(1:n_knots, function(k) {
  multi_normal_lupdf(data_in$par_values, knots[k, ], V, log = F)
})

ccc <- abs(rnorm(n_knots, sd = 0.05))
pred_out <- weights_pred_out %*% ccc
pred_in <- weights_pred_in %*% ccc

ccc <- rep(1, n_knots)
pred_sobol <- weights_sobol %*% ccc
vv <- "intercept"
plot(pred_sobol ~ data_sobol[, vv], xlab = vv,
     col = rgb(0, 0, 1, 0.3), pch = 19)
rug(knots[, vv])

yy <- range(c(pred_out, pred_in))
par(mfrow = c(1, 2))
plot(pred_out ~ data_out$par_values[, vv], ylim = yy,
     col = rgb(0, 0, 1, 0.3), pch = 19)
plot(pred_in ~ data_in$par_values[, vv], ylim = yy,
       col = rgb(1, 0, 0, 0.3), pch = 19)
par(mfrow = c(1, 1))


sum(pred_out > threshold) / length(pred_out) # mucho error
sum(pred_in <= threshold) / length(pred_in)  # poco error.

# plot(pred_in ~ data_in$overlap)

# Con sigma muy peque (widths / 10), muchos puntos quedan siempre
# lejos de los knots, por lo que son siempre cero.
# Esos puntos tienen siempre predicción cero, no importa si sus
# coef asociados son altos o no.
# Esto implica que sí o sí hay que dejar a sigma relativamente grande
# (widths / 1:3) para que no haya puntos a los que les pasa eso.
# Y algo interesante es que esto no se resuelve agregando más knots.

# El número de knots tiene un efecto similar al sigma. Más knots hacen que,
# a un mismo sigma, mu pueda llegar más arriba (porque se suman más términos.)

# los puntitos que son siempre cero suelen ser los que están bastante
# fuera de rango. Es muy obvio mirando el parámetro steps.

# plot_simulation(d)



# GP by convolution in stan ---------------------------------------

# 100 knots, previa jerárquica en los betas, sigmas = widths / 2.

smodel1 <- stan_model(
  file.path("pseudolikelihood estimation",
            "pseudolikelihood_estimation_GP_convolution_01.stan"),
  verbose = F
)

sdata1 <- list(
  N = nrow(data_in), K = n_knots,
  y = data_in$overlap, X = weights_pred_in,
  y_lower = threshold,
  prior_sigma_b_sd = 0.04
)

m1 <- sampling(
  smodel1, data = sdata1,
  chains = 4, cores = 4, iter = 2000
)
# 97.151 / 60 = 1.62 min

sm1 <- summary(m1, pars = c("sigma_b", "sigma"))[[1]]
print(sm1)
mu_hat <- as.matrix(m1, "mu") %>% t
mu_mean <- rowMeans(mu_hat)

plot(mu_mean ~ data_in$overlap); abline(0, 1, col = 2)

b_hat <- as.matrix(m1, "b") %>% t
str(b_hat)

data_in$overlap_pred <- rowMeans(weights_pred_in %*% b_hat)
data_out$overlap_pred <- rowMeans(weights_pred_out %*% b_hat)

yy <- range(c(data_in$overlap, data_in$overlap_pred,
              data_out$overlap, data_out$overlap_pred))

vv <- "wind"
par(mfrow = c(1, 2))
# predicted
plot(data_out$overlap_pred ~ data_out$par_values[, vv], ylim = yy,
     col = rgb(0, 0, 1, 0.3), pch = 19,
     main = "predictions",
     xlab = vv, ylab = "overlap")
points(data_in$overlap_pred ~ data_in$par_values[, vv], ylim = yy,
       col = rgb(1, 0, 0, 0.3), pch = 19)

# observed
plot(data_out$overlap ~ data_out$par_values[, vv], ylim = yy,
     col = rgb(0, 0, 1, 0.3), pch = 19,
     main = "observations",
     xlab = vv, ylab = "overlap")
points(data_in$overlap ~ data_in$par_values[, vv], ylim = yy,
       col = rgb(1, 0, 0, 0.3), pch = 19)
par(mfrow = c(1, 1))


# DHARMa residuals
sigma_hat <- as.matrix(m1, "sigma")
y_sim <- sapply(1:length(sigma_hat), function(i) {
  rtruncnorm(nrow(data_in), a = sdata1$y_lower, b = Inf, mean = mu_hat[, i],
             sd = sigma_hat[i])
})
str(y_sim)

res1 <- createDHARMa(observedResponse = data_in$overlap,
                     simulatedResponse = y_sim,
                     fittedPredictedResponse = rowMeans(mu_hat))
plot(res1, rank = F)
plotResiduals(res1, form = data_in$par_values[, "intercept"], rank = F)
plotResiduals(res1, form = data_in$par_values[, "wind"], rank = F)
plotResiduals(res1, form = data_in$par_values[, "steps"], rank = F)
plotResiduals(res1, form = data_in$par_values[, "vfi"], rank = T)
plotResiduals(res1, form = data_in$par_values[, "tfi"], rank = F)
plotResiduals(res1, form = data_in$par_values[, "slope"], rank = F)
# For data far from the high density of points, in the extremes, the fit is bad,
# with mu underestimated. This might be due to truncation or to the lack of
# knots in those areas, making the function approach zero (intended.)

# Partial predictions.

# set the non-plotted dimensions at the maximum, based on the map iteration.
lp <- as.matrix(m1, "lp__")[, 1]
id_map <- which.max(lp)

start <- data_in$par_values[which.max(mu_hat[, id_map]), ]
fn <- function(x) {
  # x <- start

  # evaluate the bases at x
  XX <- matrix(x, nrow = 1)
  w <- sapply(1:n_knots, function(k) {
    multi_normal_lupdf(XX, knots[k, ], V, log = F)
  })

  # compute mu
  mm <- w %*% b_hat[, id_map]
  return(mm)
}

opt <- optim(start, fn = fn, control = list(fnscale = -1))
opt

new_data <- do.call("rbind", lapply(par_names, function(v) {
  make_newdata(varying = v, data = d, mle = opt$par, ci = NULL)
}))

# evaluate bases at new_data
nd_matrix <- as.matrix(new_data[, par_names])
nd_weights <- sapply(1:n_knots, function(k) {
  multi_normal_lupdf(nd_matrix, knots[k, ], V, log = F)
})

# mu
mu_pred <- nd_weights %*% b_hat
mu_pred_ci <- apply(mu_pred, 1, hdmean, ci = 0.95) %>% t %>% as.data.frame

# plot
pred1 <- cbind(new_data, mu_pred_ci)
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

# Very bad fit for VFI, which seems to have the mode at zero. It is likely
# because I did not put knots there.
apply(knots, 2, range)
par_names

# The conditional distribution at the peak is truncated normal, which is
# probably because of basis choice. But let's see the marginals with an
# MCMC.

# function to sample at log scale (did not work)
log_like <- function(x) {
  # x <- start
  # unlog the params
  x <- params_unlog(x)

  # evaluate the bases at x
  XX <- matrix(x, nrow = 1)
  w <- sapply(1:n_knots, function(k) {
    multi_normal_lupdf(XX, knots[k, ], V, log = F)
  })

  # compute mu
  mm <- w %*% b_hat[, id_map]
  return(log(mm))
}

# function to sample at raw scale, returning -Inf if any parameter
# is out of bounds
log_like2 <- function(x) {
  if(any(x[2:n_par] < 0)) return(-Inf)
  # evaluate the bases at x
  XX <- matrix(x, nrow = 1)
  w <- sapply(1:n_knots, function(k) {
    multi_normal_lupdf(XX, knots[k, ], V, log = F)
  })

  # compute mu
  mm <- w %*% b_hat[, id_map]
  return(log(mm))
}

# MCMC at log scale
sampling_iters <- 20000
adapt_iters <- 2000
mcmc1 <- MCMC(log_like,
              n = sampling_iters + adapt_iters,
              adapt = adapt_iters,
              scale = as.numeric(apply(data_in$par_values_log, 2, sd)),
              init = params_log(start), acc.rate = 0.234)

mcmc2 <- MCMC(log_like2,
              n = sampling_iters + adapt_iters,
              adapt = adapt_iters,
              scale = as.numeric(apply(data_in$par_values, 2, sd)),
              init = start, acc.rate = 0.234)

draws1 <- params_unlog_matrix(mcmc1$samples)
colnames(draws1) <- par_names
par(mfrow = c(3, 2))
for(par in par_names) {
  plot(density(draws1[-(1:adapt_iters), par]), main = par)
}
par(mfrow = c(1, 1))

draws2 <- mcmc2$samples
colnames(draws2) <- par_names
par(mfrow = c(3, 2))
for(par in par_names) {

  if(par == "intercept") {
    den <- density(draws2[-(1:adapt_iters), par])
  } else {
    den <- density(draws2[-(1:adapt_iters), par], from = 0)
  }

  den$y <- den$y * (max(d$overlap) / max(den$y))
  plot(d$overlap ~ d$par_values[, par], main = par,
       col = rgb(0, 0, 0, 0.05))
  lines(den, main = par, lwd = 2, col = "red")

}
par(mfrow = c(1, 1))

# bad fit. It is too smooth, too normal. Maybe we need a more narrow kernel to
# make the knots





# TAREAS ------------------------------------------------------------------

# CHEQUEAR EL MODELO GRAFICANDO.

# Muestrear la likelihood con y sin incertidumbre, tanto en noise como
# en la estimación de mu.

# posibles mejoras:
# * usar una previa t sobre los coefs, para regularizar menos los extremos.
# * coeficientes más similares según proximidad (GP sobre los coefs)
# * usar todos los datos pero setear el overlap de los out_of_bounds en
#   0.001.

# Usar una distrib beta? es difícil, porque habría que restringir a mu
# para que quede entre 0 y 1. O darle muuuuy buenos inits (difícil).



# GP con intercept 0 en overlap scale -------------------------------------

# Y si probamos un GP con km fijando el intercept en 0 pero en escala del
# overlap? Creo que nunca hicimos eso.
library(DiceKriging)

m0 <- km(formula = ~ 1, coef.trend = 0,
         design = data$par_values,
         response = data$overlap,
         nugget.estim = TRUE,
         #multistart = ceiling(n_cores / 2),
         control = list(maxiter = 5e4))

# ver qué onda

fitted <- predict(m0, type = "UK",
              newdata = data$par_values,
              se.compute = TRUE, light.return = TRUE)$mean

plot(data$overlap ~ fitted)


fitted_all <- predict(m0, type = "SK",
                  newdata = d$par_values,
                  se.compute = TRUE, light.return = TRUE)$mean

plot(d$overlap ~ fitted_all)
# un sobreajuste de la gran concha de la lora.
# pero habría que plotearlo antes de descartarlo.
# no, no va a cero cuando no hay datos. Está mal.

# Go for a GAM with bounded data = 0.
