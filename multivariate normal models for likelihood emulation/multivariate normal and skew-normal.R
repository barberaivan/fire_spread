# Calibrating multivariate normal and skew-normal functions in Stan.

# On how to parameterize a correlation matrix.
# https://yingqijing.medium.com/multivariate-normal-distribution-and-cholesky-decomposition-in-stan-d9244b9aa623library(tidyverse); 


library(ggplot2); library(magrittr); theme_set(theme_bw())
library(viridis)
library(sn)
library(rstan)
library(extraDistr)  # half t density (dht)
library(DHARMa)
library(mnormt) # for pd.solve
library(trialr)
library(truncnorm)  # check truncated models

focal_folder <- "multivariate normal models for likelihood emulation"

# functions ---------------------------------------------------------------

# equal-tailed credible interval and mean
etimean <- function(x, ci = 0.95, name = "mu") {
  out <- 1 - ci
  q <- quantile(x, probs = c(out / 2, 1 - out / 2), method = 8)
  result <- c(q[1], mean(x), q[2])
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

# unnormalized log-density of a multivariate normal distribution, computed from
# the location vector xi and vcov matrix Omega. 
multi_normal_lupdf <- function(X, xi, Omega) {
  N <- nrow(X); K <- ncol(X)
  log_den <- numeric(N)
  P <- pd.solve(Omega, silent = TRUE, log.det = F) # same computation as in {sn}
  
  for(n in 1:N) {
    x_centred = X[n, ] - xi
    log_den[n] = -0.5 * t(x_centred) %*% P %*% x_centred;
  }
    
  return (log_den)
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

# Tests -------------------------------------------------------------------


# # TEST: these densities should differ from those computed with {sn} only by a 
# # constant (because ours is unnormalized)
# N <- 100; K <- 20
# x <- matrix(rnorm(N * K), N, K)
# xi <- rnorm(K)
# sigma <- runif(K, 1, 4)
# # alpha <- rep(0, K)
# alpha <- rnorm(K, sd = 2)
# # Rho <- make_corr(runif(cor_num(K), -1.5, 1.5), K)
# Rho <- rlkjcorr(1, K, eta = 1)
# Omega <- Rho * (sigma %*% t(sigma))
# P <- pd.solve(Omega)
# x[1, , drop = F] %*% P
# 
# d_hand <- multi_normal_lupdf(x, xi, Omega)
# d_sn <- dmsn(x, xi, Omega, alpha) %>% log
# 
# # check if there are Inf in my data
# sum(d_hand > -Inf) == N
# range(d_hand)
# # compare with slope = 1
# par(mfrow = c(1, 2))
# hist(d_hand, breaks = 10)
# plot(d_hand ~ d_sn)
# abline(c(coef(lm(d_hand ~ d_sn))[1], 1) )
# par(mfrow = c(1, 1))

# the sn function fails to compute log-densities that are too low.

# Watching the sn code for dmsn we could use some strategy to improve the 
# mu computation in stan. Only if it's necessary.

# ## LKJ correlation matrices
# K <- 8
# sims <- sapply(1:1000, function(i) {
#   r <- rlkjcorr(1, K, eta = 1) 
#   vecs <- as.numeric(r[lower.tri(r)])
# })
# sims2 <- sapply(1:1000, function(i) {
#   r <- rlkjcorr(1, K, eta = 0.5) 
#   vecs <- as.numeric(r[lower.tri(r)])
# })
# 
# plot(density(sims, from = -1, to = 1))
# lines(density(sims2, from = -1, to = 1), col = 2)
# # eta = 1 or 0.1 seem OK

# Simulate and fit multivariate normal curve -------------------------------

# compile stan model
smodel1 <- stan_model(file.path(focal_folder, "multivariate normal.stan"), 
                      verbose = TRUE)

N <- 3000
K <- 8

xi <- rnorm(K)
sigma_prior_sd <- 5
sigma_prior_nu <- 3
sigma <- 0.1 + rht(K, nu = sigma_prior_nu, sigma = sigma_prior_sd)
# prior for sigma
# curve(dht(x, nu = 3, sigma = sigma_prior_sd), from = 0, to = 100)
Rho <- rlkjcorr(1, K, eta = 0.8)
Omega <- Rho * (sigma %*% t(sigma))
x <- matrix(runif(N * K, -5, 5), N, K)

mu <- multi_normal_lupdf(x, xi, Omega)

# simulate tau, the residual standard deviation based on sd(mu) to get
# an R2 near 0.8
tau <- sd(mu) * 0.5
y <- rnorm(N, mu, tau)
(R2 <- var(mu) / (var(mu) + var(y - mu)))

# # visualize mu 
# ext <- 10
# dpred <- expand.grid(x1 = seq(-ext, ext, length.out = 150),
#                      x2 = seq(-ext, ext, length.out = 150))
# dpred$mu <- multi_normal_lupdf(as.matrix(dpred), xi, Omega) #%>% exp
# 
# ggplot(dpred, aes(x = x1, y = x2, z = mu)) +
#   geom_contour_filled()

sdata1 <- list(
  # data
  X = x, y = y, N = N, K = K,
  
  # prior parameters
  xi_prior_sd = 50,
  sigma_prior_sd = sigma_prior_sd,
  sigma_prior_nu = sigma_prior_nu,
  tau_prior_sd = tau * 2,
  tau_prior_nu = 3,
  sigma_lower = 0.1
)
m1 <- sampling(smodel1, data = sdata1, seed = 2123, 
               cores = 1, chains = 1, iter = 2000,
               control = list(max_treedepth = 25))
# 2154.93 / 60 = 36 min, una cadena; con N = 3000; K = 8
# Anduvo muy bien. Pasar a la skew mañana.
# Oarece que el timing cambia linealmente con el N
# 901.077 / 60 = 15 min con N = 1000, K = 8
# 8 divergent transitions, pero el ajuste es buenísimo

# Lo que más afecta el timing es la dimensión del problema, terriblemente. con 
# K = 3 no tarda nada

mu_hat <- as.matrix(m1, pars = "mu") %>% t
tau_hat <- as.matrix(m1, pars = "tau")

# dharma analysis
n_sim <- length(tau_hat)
y_sim <- sapply(1:n_sim, function(i) {
  y_ <- rnorm(N, mu_hat[, i], tau_hat[i])
}) # OK, columns are replicates
res <- createDHARMa(y_sim, y)
plot(res)

# comparar mu estimado y real con int de confianza
# luego probar qué pasa con 8 dim
mu_ci <- apply(mu_hat, 1, etimean) %>% t %>% as.data.frame()
mu_ci$mu_obs <- mu

ggplot(mu_ci, aes(x = mu_obs, y = mu_mean, ymin = mu_lower, ymax = mu_upper))+ 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  geom_ribbon(color = NA, alpha = 0.3, fill = "green") +
  geom_point(alpha = 0.1) +
  coord_fixed()
  

# multi-normal curve, truncated-normal likelihood -------------------------

# compile stan model
smodel <- stan_model(
  file.path(focal_folder, "multivariate normal - truncated.stan"), 
  verbose = TRUE
)

# non-truncated model
smodel_nt <- stan_model(
  file.path(focal_folder, "multivariate normal.stan"), 
  verbose = TRUE
)

N <- 2000 * 2 # for truncated data, make more
K <- 3

xi <- rnorm(K)
sigma_prior_sd <- 5
sigma_prior_nu <- 3
sigma <- 0.1 + rht(K, nu = sigma_prior_nu, sigma = sigma_prior_sd)
# prior for sigma
# curve(dht(x, nu = 3, sigma = sigma_prior_sd), from = 0, to = 100)
Rho <- rlkjcorr(1, K, eta = 1)
Omega <- Rho * (sigma %*% t(sigma))
x <- matrix(runif(N * K, -5, 5), N, K)

mu_et_al <- multi_normal_lupdf(x, xi, Omega)

# simulate tau, the residual standard deviation based on sd(mu) to get
# an R2 near 0.8
tau <- sd(mu) * 0.5
y <- rnorm(N, mu, tau)
(R2 <- var(mu) / (var(mu) + var(y - mu)))

# truncate y
L <- median(y)
ids_use <- which(y >= L)
N_use <- length(ids_use)
y_trunc <- y[ids_use]
mu_trunc <- mu[ids_use]

sdata <- list(
  # data
  X = x[ids_use, ], y = y[ids_use], N = N_use, K = K,
  
  # prior parameters
  xi_prior_sd = 50,
  sigma_prior_sd = sigma_prior_sd,
  sigma_prior_nu = sigma_prior_nu,
  tau_prior_sd = tau * 2,
  tau_prior_nu = 3,
  sigma_lower = 0.1,
  
  # truncation for y
  L = L
)

# fit model
m <- sampling(smodel, data = sdata, seed = 2123, 
              cores = 3, chains = 3, iter = 3000, refresh = 50,
              control = list(max_treedepth = 25, adapt_delta = 0.95))
sm <- summary(m, pars = c("xi", "Omega", "sigma", "tau"), 
              prob = c(0.025, 0.975))[[1]]
sm
traceplot(m, pars = c("sigma"), inc_warmup = TRUE)
pairs(m, pars = c("xi", "tau"))
pairs(m, pars = c("sigma", "tau"))
pairs(m, pars = c("Rho"))
pairs(m, pars = c("Omega"))
# N = 2000, K = 3, 3000 iter, adapt = 0.95: 483 / 60 = 8.05 min
#   12 div transitions, good neff (but < 1000)
#   A pesar de los problemas, anda muy bien.

# Tiene problemas porque hay un rho que es -0.87, y lo estima bien. Eso lleva a 
# que los sigmas asociados tengan una corr posterior muy alta.

# dharma analysis
mu_hat <- as.matrix(m, pars = "mu") %>% t
tau_hat <- as.matrix(m, pars = "tau")

n_sim <- length(tau_hat)
y_sim <- sapply(1:n_sim, function(i) {
  y_ <- rtruncnorm(N_use, a = L, b = Inf, mean = mu_hat[, i], sd = tau_hat[i])
}) # OK, columns are replicates
res <- createDHARMa(y_sim, y_trunc)
plot(res)

# comparar mu estimado y real con int de confianza
# luego probar qué pasa con 8 dim
mu_ci <- apply(mu_hat, 1, etimean) %>% t %>% as.data.frame()
mu_ci$mu_obs <- mu_trunc

ggplot(mu_ci, aes(x = mu_obs, y = mu_mean, ymin = mu_lower, ymax = mu_upper))+ 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  geom_ribbon(color = NA, alpha = 0.3, fill = "green") +
  geom_point(alpha = 0.1) +
  coord_fixed()


# the same for non-truncated model.
m_nt <- sampling(smodel_nt, data = sdata, seed = 2123, 
                 cores = 3, chains = 3, iter = 3000, refresh = 50,
                 control = list(max_treedepth = 25, adapt_delta = 0.95))
# N = 2000, K = 3:  266.907 / 60 = 4.45 min
# 2 diver transitions

# dharma analysis
mu_hat_nt <- as.matrix(m_nt, pars = "mu") %>% t
tau_hat_nt <- as.matrix(m_nt, pars = "tau")

n_sim <- length(tau_hat)
y_sim <- sapply(1:n_sim, function(i) {
  y_ <- rnorm(N_use, mu_hat_nt[, i], tau_hat_nt[i])
}) # OK, columns are replicates
res <- createDHARMa(y_sim, y_trunc)
plot(res)
# se ven problemas notables: underdispersion in low values y overdispersion in 
# high values.

# chech mu estimation
mu_ci_nt <- apply(mu_hat_nt, 1, etimean) %>% t %>% as.data.frame()
mu_ci_nt$mu_obs <- mu_trunc

ggplot(mu_ci_nt, aes(x = mu_obs, y = mu_mean, ymin = mu_lower, ymax = mu_upper))+ 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  geom_ribbon(color = NA, alpha = 0.3, fill = "green") +
  geom_point(alpha = 0.1) +
  coord_fixed()
# sobreestima mu para la mierda!!

# Compare fitted mu with fitted mu_nt
plot(mu_ci$mu_mean ~ mu_ci_nt$mu_mean, col = rgb(0, 0, 0, 0.1))
abline(0, 1)
# los mu estimados con el modelo truncado son mucho más bajos.
# Quizás tarde el doble en ajustar los modelos, pero sí o sí tengo que considerar
# la truncation. Si no, la likelihood va a ser muuuucho más amplia que lo que debería.


# multi-skew-normal -------------------------------------------------------

# compile stan model
smodel <- stan_model(
  file.path(focal_folder, "multivariate skew-normal - truncated.stan"), 
  verbose = F
)

# simulate data
# As extreme x values may create shifted values below -37, which
# returns -inf in Stan, we choose higher values.

N_try <- 2000
N_large <- N_try * 20
K <- 8

xi <- rnorm(K)
sigma_prior_sd <- 5
sigma_prior_nu <- 3
sigma <- 1 + rht(K, nu = sigma_prior_nu, sigma = sigma_prior_sd)
# prior for sigma
# curve(dht(x, nu = 3, sigma = sigma_prior_sd), from = 0, to = 100)
alpha <- rnorm(K, 0, 3)
Rho <- rlkjcorr(1, K, eta = 2)
Omega <- Rho * (sigma %*% t(sigma))

x <- matrix(runif(N_try * K, -2, 2), N_try, K)

# compute mu and remove those with cumulative probability = 0
mu_et_al <- multi_skew_normal_lupdf(x, xi, Omega, alpha)
ids_ok <- which(mu_et_al$shifted > -30)
N <- ifelse(length(ids_ok) > N_try, N_try, length(ids_ok))
ids_use <- sample(ids_ok, N, replace = F)
mu <- mu_et_al$mu[ids_use]
x <- x[ids_use, ]
N
# simulate tau, the residual standard deviation based on sd(mu) to get
# an R2 near 0.8
tau <- sd(mu) * 0.15
y <- rnorm(N, mu, tau)
(R2 <- var(mu) / (var(mu) + var(y - mu)))


# filter again, by truncation
L <- median(y)
trunc_use <- which(y > L)
N <- length(trunc_use)
y <- y[trunc_use]
x <- x[trunc_use, ]
mu <- mu[trunc_use]

(R2 <- var(mu) / (var(mu) + var(y - mu))) # 0.69
N

# visualize mu
hist(mu)
hist(y)

# dd <- data.frame(y = y, mu = mu, x1 = x[, 1], x2 = x[, 2])
# dpred <- expand.grid(x1 = seq(min(x[, 1]), max(x[, 1]), length.out = 150),
#                      x2 = seq(min(x[, 2]), max(x[, 2]), length.out = 150))
# dpred$mu <- b0 + multi_skew_normal_lupdf(as.matrix(dpred), xi, Omega, alpha)$mu #%>% exp
# ggplot(dpred, aes(x = x1, y = x2, z = mu)) +
#   geom_contour_filled(bins = 30) +
#   geom_point(dd, mapping = aes(color = mu))

sdata <- list(
  # data
  # 
  X = x, y = y, N = N, K = K,
  
  # prior parameters
  xi_prior_sd = 2,
  sigma_prior_sd = sigma_prior_sd,
  sigma_prior_nu = sigma_prior_nu,
  alpha_prior_sd = 3,
  tau_prior_sd = tau * 2,
  tau_prior_nu = 3,
  sigma_lower = 0.5,

  L = L
)


# initialize at reasonable values
init_f_random <- function(chain_id) {
  list(
    xi = rnorm(K, 0, 0.5),
    alpha = rnorm(K, 0, 0.5), 
    sigma = rtruncnorm(K, a = 0, mean = 1, sd = 0.5),
    tau = rtruncnorm(1, a = 0, mean = 1, sd = 0.5),
    Lcorr = t(chol(rlkjcorr(1, K, eta = 30)))
  )
}

nchains <- 1
init_list <- lapply(1:nchains, function(i) init_f_random(i))

m <- sampling(smodel, data = sdata, refresh = 20,
              cores = nchains, 
              chains = nchains, 
              iter = 10, 
              warmup = 5,
              init = init_f, #init_f_true,
              control = list(max_treedepth = 25, adapt_delta = 0.9))
# With higher K, the runif(-2, 2) default initialization fails. 
# tight initial values must be supplied.


sm <- summary(m, pars = c("xi", "Rho", "sigma", "tau"))[[1]]
sm[, c("n_eff", "Rhat")]

traceplot(m, pars = c("xi"), inc_warmup = T)
traceplot(m, pars = c("Rho", "sigma"), inc_warmup = T)
traceplot(m, pars = c("alpha"), inc_warmup = T)

pairs(m, pars = c("Rho", "sigma"))
pairs(m, pars = c("sigma", "alpha", "tau"))
pairs(m, pars = c("alpha", "xi"))

# Evaluate model

mu_hat <- as.matrix(m, pars = "mu") %>% t
tau_hat <- as.matrix(m, pars = "tau")

# dharma analysis
n_sim <- length(tau_hat)
y_sim <- sapply(1:n_sim, function(i) {
  y_ <- rtruncnorm(N, a = L, b = Inf, mean = mu_hat[, i], sd = tau_hat[i])
}) # OK, columns are replicates
res <- createDHARMa(y_sim, y)
plot(res)

mu_ci <- apply(mu_hat, 1, etimean) %>% t %>% as.data.frame()
mu_ci$mu_obs <- mu

ggplot(mu_ci, aes(x = mu_obs, y = mu_mean, ymin = mu_lower, ymax = mu_upper))+ 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  geom_ribbon(color = NA, alpha = 0.3, fill = "green") +
  geom_point(alpha = 0.1) +
  coord_fixed()

# multi-skew-normal b0 -----------------------------------------------------

# truncated skew with b0
smodel_b0 <- stan_model(
  file.path(focal_folder, "multivariate skew-normal - truncated - b0.stan"), 
  verbose = F
)

# simulate data
# As extreme x values may create shifted values below -37, which
# returns -inf in Stan, we choose higher values.

N_try <- 2000
N_large <- N_try * 20
K <- 8

xi <- rnorm(K)
sigma_prior_sd <- 5
sigma_prior_nu <- 3
sigma <- 1 + rht(K, nu = sigma_prior_nu, sigma = sigma_prior_sd)
# prior for sigma
# curve(dht(x, nu = 3, sigma = sigma_prior_sd), from = 0, to = 100)
alpha <- rnorm(K, 0, 3)
Rho <- rlkjcorr(1, K, eta = 2)
Omega <- Rho * (sigma %*% t(sigma))
b0 <- rnorm(1, 0, 3)

x <- matrix(runif(N_try * K, -2, 2), N_try, K)

# compute mu and remove those with cumulative probability = 0
mu_et_al <- multi_skew_normal_lupdf(x, xi, Omega, alpha)
ids_ok <- which(mu_et_al$shifted > -30)
N <- ifelse(length(ids_ok) > N_try, N_try, length(ids_ok))
ids_use <- sample(ids_ok, N, replace = F)
mu <- b0 + mu_et_al$mu[ids_use]
x <- x[ids_use, ]
N
# simulate tau, the residual standard deviation based on sd(mu) to get
# an R2 near 0.8
tau <- sd(mu) * 0.15
y <- rnorm(N, mu, tau)
(R2 <- var(mu) / (var(mu) + var(y - mu)))


# filter again, by truncation
L <- median(y)
trunc_use <- which(y > L)
N <- length(trunc_use)
y <- y[trunc_use]
x <- x[trunc_use, ]
mu <- mu[trunc_use]

(R2 <- var(mu) / (var(mu) + var(y - mu))) # 0.69

# visualize mu
# hist(mu)
# hist(y)

# dd <- data.frame(y = y, mu = mu, x1 = x[, 1], x2 = x[, 2])
# dpred <- expand.grid(x1 = seq(min(x[, 1]), max(x[, 1]), length.out = 150),
#                      x2 = seq(min(x[, 2]), max(x[, 2]), length.out = 150))
# dpred$mu <- b0 + multi_skew_normal_lupdf(as.matrix(dpred), xi, Omega, alpha)$mu #%>% exp
# ggplot(dpred, aes(x = x1, y = x2, z = mu)) +
#   geom_contour_filled(bins = 30) +
#   geom_point(dd, mapping = aes(color = mu))

sdata <- list(
  # data
  X = x, y = y, N = N, K = K,
  
  # prior parameters
  xi_prior_sd = 2,
  sigma_prior_sd = sigma_prior_sd,
  sigma_prior_nu = sigma_prior_nu,
  alpha_prior_sd = 3,
  tau_prior_sd = tau * 2,
  tau_prior_nu = 3,
  sigma_lower = 0.5,
  b0_prior_sd = 4,
  
  L = L
)

# initialize at reasonable values
init_f_random <- function(chain_id) {
  list(
    xi = rnorm(K, 0, 0.5),#xi,
    alpha = rnorm(K, 0, 0.5), #alpha, 
    sigma = 1 + rht(K, nu = sigma_prior_nu, sigma = sigma_prior_sd),
    # sigma = rtruncnorm(K, a = sdata$sigma_lower, mean = 3, sd = 3),
    tau = rtruncnorm(1, a = 0, mean = 1, sd = 0.5),#tau,
    Lcorr = t(chol(rlkjcorr(1, K, eta = 30))),#t(chol(Rho)),
    b0 = rnorm(1, 0, 3)
  )
}

nchains <- 4
init_list <- lapply(1:nchains, function(i) init_f_random(i))
it <- 2000; wu <- 800

m <- sampling(smodel_b0, data = sdata, refresh = 20,
              cores = nchains, 
              chains = nchains, 
              iter = it, 
              warmup = wu,
              init = init_list,
              control = list(max_treedepth = 25, adapt_delta = 0.9))
# Ahora arranca casi siempre, aunque algunas veces hay que hacer varios intentos

# 2000 iter, 8000 warm, adapt = 0.9, N = 2000, R2 = 0.8, K = 8
# 1590.71 / 60 = 26.5 min, sin problemas
 
# lo mismo, con R2 = 0.65:
# 1754.41 / 60 = 29 min

sm <- summary(m, pars = c("xi", "alpha", "sigma", "tau", "b0", "Rho"))[[1]]
sm[, c("n_eff", "Rhat")]
hist(sm[, "n_eff"] / ((it - wu) * nchains))

# Evaluate model
mu_hat <- as.matrix(m, pars = "mu") %>% t
tau_hat <- as.matrix(m, pars = "tau")

# dharma analysis
n_sim <- length(tau_hat)
y_sim <- sapply(1:n_sim, function(i) {
  y_ <- rtruncnorm(N, a = L, b = Inf, mean = mu_hat[, i], sd = tau_hat[i])
}) # OK, columns are replicates
res <- createDHARMa(y_sim, y)
plot(res)

mu_ci <- apply(mu_hat, 1, etimean) %>% t %>% as.data.frame()
mu_ci$mu_obs <- mu

ggplot(mu_ci, aes(x = mu_obs, y = mu_mean, ymin = mu_lower, ymax = mu_upper))+ 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  geom_ribbon(color = NA, alpha = 0.3, fill = "green") +
  geom_point(alpha = 0.1) +
  coord_fixed()


# multi-normal curve, truncated-normal likelihood, with b0 -----------------

# compile stan model
smodel_nb0 <- stan_model(
  file.path(focal_folder, "multivariate normal - truncated - b0.stan"), 
  verbose = F
)

N <- 2000 * 2 # for truncated data, make more
K <- 8

xi <- rnorm(K)
sigma_prior_sd <- 5
sigma_prior_nu <- 3
sigma <- 0.1 + rht(K, nu = sigma_prior_nu, sigma = sigma_prior_sd)
# prior for sigma
# curve(dht(x, nu = 3, sigma = sigma_prior_sd), from = 0, to = 100)
Rho <- rlkjcorr(1, K, eta = 1)
Omega <- Rho * (sigma %*% t(sigma))
b0 <- rnorm(1, 0, 4)

x <- matrix(runif(N * K, -5, 5), N, K)

mu <- b0 + multi_normal_lupdf(x, xi, Omega)

# simulate tau, the residual standard deviation based on sd(mu) to get
# an R2 near 0.8
tau <- sd(mu) * 0.10
y <- rnorm(N, mu, tau)
(R2 <- var(mu) / (var(mu) + var(y - mu)))

# truncate y
L <- median(y)
ids_use <- which(y >= L)
N_use <- length(ids_use)
y <- y[ids_use]
mu <- mu[ids_use]
x <- x[ids_use, ]
# new R2:
(R2 <- var(mu) / (var(mu) + var(y - mu)))

sdata <- list(
  # data
  X = x, y = y, N = N_use, K = K,
  
  # prior parameters
  xi_prior_sd = 50,
  sigma_prior_sd = sigma_prior_sd,
  sigma_prior_nu = sigma_prior_nu,
  tau_prior_sd = tau * 2,
  tau_prior_nu = 3,
  sigma_lower = 0.1,
  b0_prior_sd = 4, 
  
  # truncation for y
  L = L
)

init_f_random <- function(chain_id) {
  list(
    xi = rnorm(K, 0, 0.5),#xi,
    alpha = rnorm(K, 0, 0.5), #alpha, 
    sigma = 1 + rht(K, nu = sigma_prior_nu, sigma = sigma_prior_sd),
    # sigma = rtruncnorm(K, a = sdata$sigma_lower, mean = 3, sd = 3),
    tau = rtruncnorm(1, a = 0, mean = 1, sd = 0.5),#tau,
    Lcorr = t(chol(rlkjcorr(1, K, eta = 30))),#t(chol(Rho)),
    b0 = rnorm(1, 0, 3)
  )
}

nchains <- 4
init_list <- lapply(1:nchains, function(i) init_f_random(i))

m <- sampling(smodel_nb0, data = sdata, refresh = 20,
              cores = nchains, 
              chains = nchains, 
              iter = 2000, 
              warmup = 800,
              init = init_list,
              control = list(max_treedepth = 25, adapt_delta = 0.9))
# está terriblemente lento, pero antes no había ocurrido eso.
# o sea, puede que dependa de los datos.

# 16702.1 / 3600 = 4.63 h; 41 divergent transitions.
# R2 = 0.75; N = 2000
# A pesar de los problemas, anduvo perfecto.

# Quizás el problema fue la escala de los parámetros. Por ej, tau ~ 600.
# los xi casi que copiaron la previa, con posteriores absurdamente amplias.
# anduvo bien, pero costó. sería más fácil escalar todo.
# En fin, anduvo bien. El problema fue que al generar datos con x en -5 , 5,
# los mu se hicieron muuuuuy bajos (eso pasa al aumentar las dims.)
# entonces tau también fue grande, terriblemente.

sm <- summary(m, pars = c("xi", "sigma", "tau", "b0", "Omega"), 
              prob = c(0.025, 0.975))[[1]]
sm
traceplot(m, pars = c("sigma"), inc_warmup = TRUE)
pairs(m, pars = c("xi", "tau"))
pairs(m, pars = c("sigma", "tau"))
pairs(m, pars = c("Rho"))
pairs(m, pars = c("Omega"))

# dharma analysis
mu_hat <- as.matrix(m, pars = "mu") %>% t
tau_hat <- as.matrix(m, pars = "tau")

n_sim <- length(tau_hat)
y_sim <- sapply(1:n_sim, function(i) {
  y_ <- rtruncnorm(N_use, a = L, b = Inf, mean = mu_hat[, i], sd = tau_hat[i])
}) # OK, columns are replicates
res <- createDHARMa(y_sim, y)
plot(res)

# comparar mu estimado y real con int de confianza
# luego probar qué pasa con 8 dim
mu_ci <- apply(mu_hat, 1, etimean) %>% t %>% as.data.frame()
mu_ci$mu_obs <- mu

ggplot(mu_ci, aes(x = mu_obs, y = mu_mean, ymin = mu_lower, ymax = mu_upper))+ 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  geom_ribbon(color = NA, alpha = 0.3, fill = "green") +
  geom_point(alpha = 0.1) +
  coord_fixed()




# Tarea: ------------------------------------------------------------------

# Explorar correlación entre b0, xi, alpha...
# Tomar el tiempo de la multi-normal con b0 y la multi-skew con b0.
# ambas truncadas.

# NOTAS -------------------------------------------------------------------

# ESTIMABILIDAD

# Al simular datos truncados, el R2 puede quedar muy bajo, y en ese caso, el
# modelo puede ser casi inestimable. Supongo que agregando un b0 eso podría 
# mejorar.
 
# RESTRICCIONES Y PREVIAS

# Idealmente, las predictoras deberían estar estandarizadas en el training set,
# para que la previa sobre sigma permita funciones planas. 

# el xi para el intercept y para los parámetros de vegetación debería estar 
# restringido al rango de los datos. Sin estas restricciones corremos el riesgo 
# de que la likelihood suba estúpidamente fuera de rango.
# para los demás parámetros de propagación (predictoras), el xi
# puede quedar fuera del rango de los datos hacia el lado restringido, para 
# permitir que la likelihood tenga un pico alrededor de cero. 
# Pensar bien esta situación de simetrías. quizás sea imposible restringir
# bien la forma en el caso multivariado, ya que los xi y los alpha significan
# cosas misteriosas.

# Estas previas serían buenas:
# xi: centrada en la media de los datos y con truncamiento dependiendo de la 
#     predictora
# sigma: truncada por abajo, porque desconfiamos de likelihoods super puntudas
# alpha: regularizadora centrada en cero, ya que no creemos en formas 
#        terriblemente asimétricas que anulan la p en xi.
# tau: regularizadora. El paisaje nos deja saber qué valores son razonables.
# rho: lkj(eta = 3)? Algo regularizador, no hay motivos para pensar en alta 
#      corr si excluimos el fwi.
# b0: en escala logit, puede estar centrado en 0.5, y quizás debería ser amplio.
#     hay que ver cuánto escala el modo en función de los alpha. 
#     cuando el modo se corre mucho de xi, b0 va a estar lejos del modo, 
#     pudiendo tomar valores muy bajos.

# AJUSTE

# cuando no usamos truncación pero los datos están truncados, los dharma son muy
# altos para las predicciones bajas del modelo. 

# con datos skew, encontré que los dharma eran muy altos en las predicciones
# altas, ajustando un modelo non-skew.
