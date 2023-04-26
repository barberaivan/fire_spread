# Calibrating multivariate normal and skew-normal functions in Stan.

# On how to parameterize a correlation matrix.
# https://yingqijing.medium.com/multivariate-normal-distribution-and-cholesky-decomposition-in-stan-d9244b9aa623library(tidyverse); 


library(ggplot2); library(magrittr); theme_set(theme_bw())
library(viridis)
library(sn)
library(cmdstanr)   ## changed rstan to use profile() function
library(extraDistr)  # half t density (dht)
library(DHARMa)
library(mnormt) # for pd.solve
library(trialr)

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
  
  for(n in 1:N) {
    x_centred <- X[n, ] - xi
    # normal log-density
    log_den <- -0.5 * t(x_centred) %*% P %*% x_centred;
    # normal cumulative probability for shifted value
    log_cum <- pnorm(alpha %*% (x_centred / sigma)) %>% log
    # skew normal log-density
    log_den_skew[n] <- log(2) + log_den + log_cum
  }
  
  return (log_den_skew)
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


# Tests -------------------------------------------------------------------

# inverse function: detect dimension of the matrix from the length of the 1D 
# vector of correlation coefficients
cor_dim <- function(l) 1 / 2 * (sqrt(8 * l + 1) + 1)
# using wolframalpha: https://www.wolframalpha.com/widgets/view.jsp?id=c86d8aea1b6e9c6a9503a2cecea55b13

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

N <- 3000
K <- 3

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
range(mu)

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
  X = t(x), y = y, N = N, K = K,
  
  # prior parameters
  xi_prior_sd = 50,
  sigma_prior_sd = sigma_prior_sd,
  sigma_prior_nu = sigma_prior_nu,
  tau_prior_sd = tau * 2,
  tau_prior_nu = 3
)

# compile stan model
smodel1 <- cmdstan_model(file.path(focal_folder, "multivariate normal_profiling.stan"))

# sample
m1 <- smodel1$sample(data = sdata1, 
                     chains = 1, iter_warmup = 5, iter_sampling = 100)

m1$profiles()[[1]] # calcular mu es lo que lleva todo el tiempo.
# particularmente, el loop es lo que tarda una locura.
tt <- m1$profiles()[[1]][, "total_time"]
which.min(tt)
tt / tt[1] # 0.2185
tt[which.min(tt)] / tt[1] # 0.2139

# check function 6 is the same as f1
sss <- m1$summary(variables = c("mu", "mu_12"))
all.equal(sss$mean[1:3000], sss$mean[3001:6000]) # OK
all.equal(sss$median[1:3000], sss$median[3001:6000]) # OK
all.equal(sss$q95[1:3000], sss$q95[3001:6000]) # OK

plot(sss$mean[1:3000], sss$mean[3001:6000]);abline(0, 1)
plot(sss$median[1:3000], sss$median[3001:6000]);abline(0, 1)
plot(sss$q95[1:3000], sss$q95[3001:6000]);abline(0, 1)

# La función nro 10 demostró ser la más rápida, aunque no es una locura.
# 2000 iter pueden tardar around 30 min para K = 8


### TEST CLEAN STAN CODE
model_clean <- stan_model(file.path(focal_folder, "multivariate normal.stan"))
sdata <- list(
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
m_clean <- sampling(model_clean, data = sdata, seed = 2123, 
                    cores = 1, chains = 1, iter = 2000,
                    control = list(max_treedepth = 25))
# 2154.93 / 60 = 36 min, una cadena; con N = 3000; K = 8
# Anduvo muy bien. Pasar a la skew mañana.
# Oarece que el timing cambia linealmente con el N
# 901.077 / 60 = 15 min con N = 1000, K = 8
# 8 divergent transitions, pero el ajuste es buenísimo
# 
# Lo que más afecta el timing es la dimensión del problema, terriblemente. con 
# K = 3 no tarda nada
m1 <- m_clean

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

  
