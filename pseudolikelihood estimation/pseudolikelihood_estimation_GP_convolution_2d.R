# Try to fit a GP by convolution. The intercept will be
# fixed at zero, so when there is no data the GP tends to zero. Fit at the
# overlap scale.

### IDEA:
# to avoid the problem of space-filling,
# fit a convoluted GP but only considering pair-wise interactions. this way
# it would be like a spline, but with constrained shape.

# pair-wise interactions:
n_par <- 6
(n_interact <- ncol(combn(1:n_par, 2)))
side_dim <- 3
side_dim ^ 2
(n_interact * side_dim ^ 2)
n_par * 10 # 60 marginal basis

# 6 marginal sigmas
# 1 interaction sigma (for all pairs)

# se podría ajustar con lm() y ya? sin restricciones?
# raro.

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

props_in_range <- function(p, range) {
  p * (range[2] - range[1]) + range[1]
}


unit_to <- function(x) (x - min(x)) / (max(x) - min(x))
unit_from <- function(x_unit, x) x_unit * max(x - min(x)) + min(x)

dnorm_un <- function(x, mu, sigma) {
  exp(-((x - mu) / sigma) ^ 2)
}


# Function to make marginal bases (gaussian).
# data: matrix[N, P] or data.frame containing the data to evaluate, with cases
#   in rows and variables in columns.
# knots: matrix[K, ncol(data)] defining the knots for each variable.
# sigmas: vector[ncol(data)] defining the sigma for each variable.
bases_marginal <- function(data, knots, sigmas, names = par_names) {

  # ## test
  # v1 <- seq(-7.5, 7.5, length.out = 100)
  # v2 <- seq(0, 12, length.out = 100)
  # dd <- expand.grid(
  #   intercept = v1,
  #   tfi = v2
  # )
  # data <- dd
  # names <- c("intercept", "tfi")
  # knots <- knots_interact[, names]
  # sigmas <- sigmas_interact[names]
  # ##

  if(is.null(dim(data))) data <- matrix(data, nrow = 1)
  k_marginal <- nrow(knots)
  v <- ncol(data)

  bb <- do.call("cbind", lapply(1:v, function(vv) {
    xx <- do.call("cbind", lapply(1:k_marginal, function(k) {
      dnorm_un(data[, vv], knots[k, vv], sigmas[vv])
    }))
    return(xx)
  }))

  # name the bases
  if(!is.null(names)) {
    gg <- expand.grid(kk = 1:k_marginal, nn = names)
    nn <- paste(gg$nn, gg$kk, sep = "_")
    colnames(bb) <- nn
  }

  return(bb)
}

# Function to make 2D bases (gaussian) by convolving marginal bases.
# data: matrix[N, P] or data.frame containing the data to evaluate, with cases
#   in rows and variables in columns.
# knots: matrix[K, ncol(data)] defining the knots for each variable.
# sigmas: vector[ncol(data)] defining the sigma for each variable.
bases_interact <- function(data, knots, sigmas, names = par_names) {

  # ## test
  # v1 <- seq(-7.5, 7.5, length.out = 100)
  # v2 <- seq(0, 12, length.out = 100)
  # dd <- expand.grid(
  #   intercept = v1,
  #   tfi = v2
  # )
  # data <- dd
  # names <- c("intercept", "tfi")
  # knots <- knots_interact[, names]
  # sigmas <- sigmas_interact[names]
  # ##

  if(is.null(dim(data))) data <- matrix(data, nrow = 1)
  k_marginal <- nrow(knots)
  v <- ncol(data)

  # marginal bases for every variable
  bm <- bases_marginal(data, knots, sigmas, names)

  # turn into list by variable
  ii <- rep(1:v, each = k_marginal)
  blist <- lapply(1:v, function(vv) bm[, ii == vv])

  # pair-wise interactions
  v_pairs <- combn(1:v, 2) %>% t
  n_pairs <- nrow(v_pairs)

  # ids for all pair-wise basis in each interaction
  ids_convolve <- as.matrix(expand.grid(
    v1 = 1:k_marginal,
    v2 = 1:k_marginal
  ))
  k_interact <- k_marginal ^ 2

  bb <- do.call("cbind", lapply(1:n_pairs, function(v) {
    v1 <- v_pairs[v, 1]
    v2 <- v_pairs[v, 2]

    xx <- sapply(1:k_interact, function(k) {
      blist[[v1]][, ids_convolve[k, 1]] * blist[[v2]][, ids_convolve[k, 2]]
    })

    return(xx)
  }))

  # name the bases
  if(!is.null(names)) {
    pairs_names <- paste(names[v_pairs[, 1]], names[v_pairs[, 2]], sep = ".")
    nums_names <- paste(ids_convolve[, 1], ids_convolve[, 2], sep = ".")
    cc <- expand.grid(num = nums_names, var = pairs_names)
    cn <- paste(cc$var, cc$num, sep = "_")
    colnames(bb) <- cn
  }

  return(bb)
}


# Understanding gaussian bases --------------------------------------------

# use intercept and steps as reference
data_in <- d[d$overlap > 0.2, ]

# ranges <- apply(data_in$par_values, 2, range)
ranges <- rbind(
  apply(data_in$par_values, 2, hdi_lower, ci = 0.98),
  apply(data_in$par_values, 2, hdi_upper, ci = 0.98)
)
widths <- apply(ranges, 2, diff)

x1 <- data_in$par_values[, "intercept"]
x2 <- data_in$par_values[, "steps"]

plot(x1 ~ unit_to(x1))
plot(unit_from(ppoints(20), x1) ~ ppoints(20))

# marginal bases. Use sigma = 1 / n_knots
n_knots <- 10

knots1 <- unit_from(ppoints(n_knots), x1)
knots2 <- unit_from(ppoints(n_knots), x2)

s_factor <- 1
s1 <- widths["intercept"] / (n_knots * s_factor)
s2 <- widths["steps"] / (n_knots * s_factor)

# marginal bases for x1 and x2 evaluated over sequences

nseq <- 100
xseq1 <- seq(min(x1), max(x1), length.out = nseq)
xseq2 <- seq(min(x2), max(x2), length.out = nseq)

X1 <- do.call("cbind", lapply(1:n_knots, function(i) {
  dnorm_un(xseq1, knots1[i], s1)
}))

X2 <- do.call("cbind", lapply(1:n_knots, function(i) {
  dnorm_un(xseq2, knots2[i], s2)
}))

matplot(xseq1, X1, type = "l")
matplot(xseq2, X2, type = "l")

# marginal curves
b <- rnorm(n_knots) %>% abs
f1 <- X1 %*% b
plot(f1 ~ xseq1, type = "l", ylim = c(0, max(f1)))
# when s_factor is small (0.1), the kernel is too wide, and possible curves
# are all very very smooth and quadratic.


# Bivariate bases: convolution.
rownames(X1) <- xseq1
rownames(X2) <- xseq2
i <- sample(1:n_knots, 1)
j <- sample(1:n_knots, 1)
Z <- outer(X1[, i], X2[, j], "*")
image(xseq1, xseq2, Z)
# outer and image understand themselves, although outer puts the first argument
# in rows, while image() takes the rows in the x-axis.


# redefine bases
n_knots <- 4
knots1 <- props_in_range(seq(0, 1, length.out = n_knots), ranges[, "intercept"])
knots2 <- props_in_range(seq(0, 1, length.out = n_knots), ranges[, "steps"])

s_factor <- 1.1
s1 <- widths["intercept"] / (n_knots * s_factor)
s2 <- widths["steps"] / (n_knots * s_factor)

xseq1 <- seq(ranges[1, "intercept"], ranges[2, "intercept"], length.out = 100)
xseq2 <- seq(ranges[1, "steps"], ranges[2, "steps"], length.out = 100)

xgrid <- expand.grid(x1 = xseq1, x2 = xseq2)
kgrid <- expand.grid(x1 = knots1, x2 = knots2)

XX1 <- do.call("cbind", lapply(1:n_knots, function(i) {
  dnorm_un(xgrid$x1, knots1[i], s1)
}))

XX2 <- do.call("cbind", lapply(1:n_knots, function(i) {
  dnorm_un(xgrid$x2, knots2[i], s2)
}))

# convolve
XX_conv <- do.call("cbind", lapply(1:n_knots, function(i) {
  do.call("cbind", lapply(1:n_knots, function(j) {
    XX1[, i] * XX2[, j]
  }))
}))

b <- abs(rnorm(n_knots ^ 2))
# b <- rep(1, n_knots ^ 2)
f2 <- XX_conv %*% b
image(xseq1, xseq2, matrix(f2, nrow = nseq))
points(x2 ~ x1, data = kgrid, pch = 19)
# rug(x1, side = 1)
# rug(x2, side = 2)

# try with k_side = 4 for bivariate bases


# Knots definition ---------------------------------------------------

k <- 10
thres <- 0.2
vv <- "vfi"

knots_98 <- props_in_range(seq(0, 1, length.out = k),
                           range = c(hdi_lower(d$par_values[, vv], ci = 0.98),
                                     hdi_upper(d$par_values[, vv], ci = 0.98)))

knots_98_1 <- props_in_range(
  seq(0, 1, length.out = k),
  range = c(hdi_lower(d$par_values[d$overlap > thres, vv], ci = 0.98),
            hdi_upper(d$par_values[d$overlap > thres, vv], ci = 0.98))
)

plot(d$overlap ~ d$par_values[, vv], pch = 19, col = rgb(0, 0, 0, 0.1),
     xlab = vv, ylab = "overlap")
# abline(v = knots_98)
abline(v = knots_98_1, col = "red")
abline(h = thres, lwd = 2, col = "blue")

# for steps it is necessary to use threshold = 0.2.
# otherwise, too many knots are necessary.

# use a uniform sequence of points within the 98 % hdi of points > threshold.
knots_bounds <- sapply(par_names, function(par) {
  c(lower = hdi_lower(d$par_values[d$overlap > thres, par], ci = 0.98),
    upper = hdi_upper(d$par_values[d$overlap > thres, par], ci = 0.98))
})
widths <- apply(knots_bounds, 2, diff)

sigma_factor_marg <- 1
sigma_factor_int <- 1.1
k_marginal <- 10
k_interact_side <- 4
k_interact <- k_interact_side ^ 2
sigmas_marginal <- widths / (k_marginal * sigma_factor_marg)
sigmas_interact <- widths / (k_interact_side * sigma_factor_int)

knots_marginal <- sapply(par_names, function(par) {
  props_in_range(seq(0, 1, length.out = k_marginal),
                 range = knots_bounds[, par])
})

knots_interact <- sapply(par_names, function(par) {
  props_in_range(seq(0, 1, length.out = k_interact_side),
                 range = knots_bounds[, par])
})


# Create bases for fitting ------------------------------------------------

data_in <- d[d$overlap > thres, ]

# evaluate bases on data (subset)
bm <- bases_marginal(data_in$par_values, knots_marginal, sigmas_marginal,
                     par_names)
bi <- bases_interact(data_in$par_values, knots_interact, sigmas_interact,
                     par_names)

# (on all data)
bm_all <- bases_marginal(d$par_values, knots_marginal, sigmas_marginal,
                         par_names)
bi_all <- bases_interact(d$par_values, knots_interact, sigmas_interact,
                         par_names)

# evaluate bases on sobol sequence over the parameter space
ranges <- apply(d$par_values, 2, range)
pp <- sobol(3000, dim = n_par, init = TRUE, seed = 2134)
data_fill <- sapply(1:n_par, function(v) {
  qunif(pp[, v], ranges[1, v], ranges[2, v])
})
colnames(data_fill) <- par_names

bm_fill <- bases_marginal(data_fill, knots_marginal, sigmas_marginal,
                          par_names)
bi_fill <- bases_interact(data_fill, knots_interact, sigmas_interact,
                          par_names)

# Prior check -------------------------------------------------------------

p_marg <- ncol(bm)
p_int <- ncol(bi)
v_id <- rep(1:n_par, each = k_marginal)

# marginal bases coefficients

# check marginal predictions
bmarg <- rtruncnorm(p_marg, a = 0, mean = 0, sd = 0.2)
fmarg <- bm_fill %*% bmarg

par(mfrow = c(3, 2))
for(par in par_names) {
  plot(fmarg ~ data_fill[, par], ylab = "y", xlab = par,
       pch = 19, col = rgb(0, 0, 0, 0.1))
  abline(h = c(0, 1), col = "red")
  abline(v = knots_bounds[, par], col = "blue")
}
par(mfrow = c(1, 1))

# check partial predictions
bmarg <- rtruncnorm(p_marg, a = 0, mean = 0, sd = 1)

par(mfrow = c(3, 2))
for(par in 1:n_par) {
  ff <- bm_fill[, v_id == par] %*% bmarg[v_id == par]
  plot(ff ~ data_fill[, par], ylab = "y", xlab = par_names[par],
       pch = 19, col = rgb(0, 0, 0, 0.1))
  abline(h = c(0, 1), col = "red")
  abline(v = knots_bounds[, par], col = "blue")
}
par(mfrow = c(1, 1))

# sigma for marginal bases could be near 0.2. 1 is terrible.
# Normal(0, 0.3) would be a reasonable prior.


# interaction bases coefficients

# check marginal predictions
bint <- rtruncnorm(p_int, a = 0, mean = 0, sd = 0.2)
fint <- bi_fill %*% bint

par(mfrow = c(3, 2))
for(par in par_names) {
  plot(fint ~ data_fill[, par], ylab = "y", xlab = par,
       pch = 19, col = rgb(0, 0, 0, 0.1))
  abline(h = c(0, 1), col = "red")
  abline(v = knots_bounds[, par], col = "blue")
}
par(mfrow = c(1, 1))

# check pseudo-partial predictions
bint <- rtruncnorm(p_int, a = 0, mean = 0, sd = 0.2)
# bint <- rep(1, p_int)
par(mfrow = c(3, 2))
for(par in 1:n_par) {
  ids <- grep(par_names[par], colnames(bi_fill))
  ff <- bi_fill[, ids] %*% bint[ids]
  plot(ff ~ data_fill[, par], ylab = "y", xlab = par_names[par],
       pch = 19, col = rgb(0, 0, 0, 0.1))
  abline(h = c(0, 1), col = "red")
  abline(v = knots_bounds[, par], col = "blue")
}
par(mfrow = c(1, 1))


# m1: bases marg y bases int, k = 4, solo sigma_marg y sigma_int -----------

# 100 knots, previa jerárquica en los betas, sigmas = widths / 2.

smodel1 <- stan_model(
  file.path("pseudolikelihood estimation",
            "pseudolikelihood_estimation_GP_convolution_2d.stan"),
  verbose = F
)

sdata1 <- list(
  N = nrow(data_in), P_marg = ncol(bm), P_int = ncol(bi),
  y = data_in$overlap, X_marg = bm, X_int = bi,
  y_lower = thres,
  prior_sigma_b_marg_sd = 0.3,
  prior_sigma_b_int_sd = 0.2
)

m1 <- sampling(
  smodel1, data = sdata1, refresh = 20,
  chains = 4, cores = 4, iter = 2000
  # chains = 1, cores = 1, iter = 10
)
# con k_marg = 10, k_int = 4, sigma_factor_int = 1.1 y prior_sd en 0.3 y 0.2:
# 567.631 / 60 =  9.46 min, pero no converge ni a palos.


# m2: only marginal bases, linear model, separate sigmas ------------------

smodel2 <- stan_model(
  file.path("pseudolikelihood estimation",
            "pseudolikelihood_estimation_GP_marginals_sums.stan"),
  verbose = F
)

# matrix to repeat sigmas by dimension
aaa <- factor(rep(par_names, each = k_marginal), levels = par_names)
S <- model.matrix(~ aaa - 1)

sdata2 <- list(
  N = nrow(data_in), K = ncol(bm), X = bm, D = n_par, S = S,
  y = data_in$overlap,
  y_lower = thres,
  prior_sigma_b_sd = 0.3
)

m2 <- sampling(
  smodel2, data = sdata2, refresh = 20,
  chains = 4, cores = 4, iter = 2000
  # chains = 1, cores = 1, iter = 10
)
# muy rápido, sin problemas
# 114.112 / 60 = 1.9 min


sm2 <- summary(m2, pars = c("sigma_b", "sigma"))[[1]]
print(sm2)
mu_hat <- as.matrix(m2, "mu") %>% t
mu_mean <- rowMeans(mu_hat)

plot(mu_mean ~ data_in$overlap); abline(0, 1, col = 2)
# muchísimo mejor que antes!

b_hat <- as.matrix(m2, "b") %>% t
str(b_hat)

# check fit inside and outside fitted data
data_in$overlap_pred <- mu_mean
data_out <- d[d$overlap <= thres, ]
bm_out <- bases_marginal(data_out$par_values, knots_marginal,
                         sigmas_marginal, par_names)
data_out$overlap_pred <- rowMeans(bm_out %*% b_hat)

yy <- range(c(data_in$overlap, data_in$overlap_pred,
              data_out$overlap, data_out$overlap_pred))

vv <- "vfi"
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
sigma_hat <- as.matrix(m2, "sigma")
y_sim <- sapply(1:length(sigma_hat), function(i) {
  rtruncnorm(nrow(data_in), a = sdata2$y_lower, b = Inf, mean = mu_hat[, i],
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
# not so good

# Partial predictions.

# set the non-plotted dimensions at the maximum, based on the map iteration.
lp <- as.matrix(m2, "lp__")[, 1]
id_map <- which.max(lp)

start <- data_in$par_values[which.max(mu_hat[, id_map]), ]
fn <- function(x) {
  # x <- start

  # evaluate the bases at x
  XX <- matrix(x, nrow = 1)
  bb <- bases_marginal(XX, knots_marginal, sigmas_marginal, names = NULL)

  # compute mu
  mm <- bb %*% b_hat[, id_map]
  return(mm)
}

opt <- optim(start, fn = fn, control = list(fnscale = -1))
opt

new_data <- do.call("rbind", lapply(par_names, function(v) {
  make_newdata(varying = v, data = d, mle = opt$par, ci = NULL)
}))

# evaluate bases at new_data
nd_matrix <- as.matrix(new_data[, par_names])
nd_bases <- bases_marginal(nd_matrix, knots_marginal, sigmas_marginal,
                           names = par_names)

# mu
mu_pred <- nd_bases %*% b_hat
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
# terrible? Or is it the partial prediction?


# predictions for all data points (conditional)
bm_all <- bases_marginal(d$par_values, knots_marginal, sigmas_marginal)
d$overlap_pred <- rowMeans(bm_all %*% b_hat)

d_long <- cbind(overlap = d$overlap,
                overlap_pred = d$overlap_pred,
                as.data.frame(d$par_values))
d_long <- pivot_longer(d_long, all_of(which(names(d_long) %in% par_names)),
                       values_to = "varying_val", names_to = "varying_var")

ggplot() +
  geom_point(data = d_long, mapping = aes(x = varying_val, y = overlap),
             inherit.aes = F, alpha = 0.05) +
  geom_point(data = d_long, mapping = aes(x = varying_val, y = overlap_pred),
             inherit.aes = F, alpha = 0.05, color = "red") +
  facet_wrap(vars(varying_var), ncol = 2, scales = "free_x",
             strip.position = "bottom") +
  ylab("overlap") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text.x = element_text(margin = margin(t = 0), vjust = 1),
        axis.title.x = element_blank(),
        panel.spacing.y = unit(5, "mm"))

# no está tan terriblemente mal. O sea, está muy mal donde no hay datos,
# pero fuera de eso, no. El problema es que al no haber interacciones
# cualquier valor solito puede tener sentido, y no es para nada así en la realidad.

# m3: only marginal bases, multiplicative model, separate sigmas ----------

# bases matrix into array
bases_dim <- rep(1:n_par, each = k_marginal)
bm_fill_array <- array(NA, dim = c(nrow(bm_fill), ncol(bm_fill) / n_par, n_par))
for(dd in 1:n_par) {
  bm_fill_array[, , dd] <- bm_fill[, bases_dim == dd]
}

# prior check
sss <- rtruncnorm(n_par, a = 0, mean = 0, sd = 1)
b <- matrix(abs(rnorm(n_par * k_marginal, sd = rep(sss, each = k_marginal))),
                ncol = n_par)

smooths <- do.call("cbind", lapply(1:n_par, function(dd) {
  bm_array[, , dd] %*% b[, dd]
}))
fmarg <- apply(smooths, 1, prod)

par(mfrow = c(3, 2))
for(par in par_names) {
  plot(fmarg ~ data_in$par_values[, par], ylab = "y", xlab = par,
       pch = 19, col = rgb(0, 0, 0, 0.1))
  abline(h = c(0, 1), col = "red")
  abline(v = knots_bounds[, par], col = "blue")
}
par(mfrow = c(1, 1))
sss

# model

smodel3 <- stan_model(
  file.path("pseudolikelihood estimation",
            "pseudolikelihood_estimation_GP_marginals_prods.stan"),
  verbose = F
)

bm_array <- array(NA, dim = c(nrow(bm), ncol(bm) / n_par, n_par))
for(dd in 1:n_par) {
  bm_array[, , dd] <- bm[, bases_dim == dd]
}

# array for stan!
bm_array_stan <- aperm(bm_array, c(3, 1, 2))
bm_array_stan %>% str

sdata3 <- list(
  N = nrow(data_in), Kd = k_marginal, D = n_par,
  X = bm_array_stan,
  y = data_in$overlap,
  y_lower = thres,
  prior_sigma_b_sd = 1
)

m3 <- sampling(
  smodel3, data = sdata3, refresh = 20,
  chains = 4, cores = 4, iter = 2000
  # chains = 1, cores = 1, iter = 10
)
# Chain 4:                823.839 seconds (Total)
# Chain 4:
#   Warning messages:
#   1: There were 3982 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
# https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
# 2: Examine the pairs() plot to diagnose sampling problems
#
# 3: The largest R-hat is NA, indicating chains have not mixed.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#r-hat
# 4: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#bulk-ess
# 5: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# https://mc-stan.org/misc/warnings.html#tail-ess
sm3 <- summary(m3)[[1]]
hist(sm3[, "Rhat"])
# quizás se pueda???
View(sm3)

# Está mal pero no taaan mal. quizás se pueda simplificar un poco.

# Old code que puede servir ----------------------------------------------




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





# TAREAS y NOTAS -------------------------------------------------------------


# Ya hecho:

# - Ajustar modelo solo con 10 bases marginales, sumadas, variando el sigma
#   por dimensión.
#   Esto anduvo mal (m2). Habría que ver qué pasa si incluimos en el ajuste
#   todos los datos. Pero parece que el problema es que los efectos sean
#   puramente aditivos.
#

# Por hacer:

#  - Probar modelo solo con bases de interacción? con K = 4? Regularizando más el
#    sigma? (m1 simplificado)

# - Marginales pero multiplicando las bases, como si fueran realmente densities.
#   m3. intentar hacerle converger. Quizás con simplexes?


# Tareas extra.

# Muestrear la likelihood con y sin incertidumbre, tanto en noise como
# en la estimación de mu.

# posibles mejoras:

# * usar una previa t sobre los coefs, para regularizar menos los extremos.
# * coeficientes más similares según proximidad (GP sobre los coefs)
# * usar todos los datos pero setear el overlap de los out_of_bounds en
#   0.001.
# * usar una distrib beta? es difícil, porque habría que restringir a mu
#   para que quede entre 0 y 1. O darle buenos inits (difícil).

## Abandono todo esto porque está demasiado complicado. Veamos
## cuánto tarda en ajustarse un GP común. Puedo combinar el GP con una
## density ajustada a los datos para que esa density limite dónde anular
## la log posterior.

# GP again, how much does it take? ----------------------------------------

library(GauPro)

gp <- gpkm(
  X = d$par_values,
  Z = d$overlap,
  parallel = F, useC = TRUE, nug.max = 100,
  kernel = "matern52"
) # 5 min para 1800 obs. No es tanto.

# define loglik_high associated with this new GP and the new loglik_threshold,
# to be used in the next wave
fitted_ll <- loglik_model$pred(particles_data_join$par_values[rows_fit, ],
                               se.fit = F)
loglik_high <- max(fitted_ll)
loglik_threshold_new <- get_loglik_threshold(loglik_high,
                                             loglik_tolerance,
                                             loglik_lower)

# find MLE in the fitted GP
if(gp_optimize) {
  # max fitted parameters value
  id_fitted <- particles_data_join$particle_id[rows_fit]
  id_max <- id_fitted[which.max(fitted_ll)]
  id_filter <- which(particles_data_join$particle_id == id_max)
  params_max_fitted <- particles_data_join$par_values[id_filter, ]

  # bounds
  lowers <- apply(particles_data_join$par_values[rows_fit, ], 2, min)
  uppers <- apply(particles_data_join$par_values[rows_fit, ], 2, max)

  message("Optimizing pseudo-likelihood surface")
  op <- optim(params_max_fitted, gp_predict, loglik_model = loglik_model,
              se.fit = FALSE,
              control = list(fnscale = -1, maxit = 1e5),
              method = "L-BFGS-B", lower = lowers, upper = uppers)
}

