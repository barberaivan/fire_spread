# Modelo completamente jerárquico.
#   Todos los parámetros dependen del FWI,
#   se estima una matriz de varianza-covarianza para los errores,
#   se incluyen los steps de muchos fuegos que no se usaron en propagación,
#   area tiene likelihood normal truncada.

# Tareas
#   - Implementar gibbs sampler para regresión comunacha.
#   - Implementar transform para los steps (ver jacobian).
#   - Implementar update para la vcov matrix a estimar.
#   - Implementar gibbs para normal truncada.



# packages ----------------------------------------------------------------

library(tidyverse)
library(rstan)
library(LaplacesDemon) # rinvwishart
library(invgamma)
library(MASS)  # mvrnorm
library(trialr) # rlkjcorr
library(truncnorm)
library(truncreg) # tnorm regression

# functions ---------------------------------------------------------------


# Inverse whisart ---------------------------------------------------------

d <- 20
N <- 2000
nu <- d + 1
s2 <- 10 ^ 2
vcovs <- array(NA, dim = c(d, d, N))
corrs <- array(NA, dim = c(d, d, N))

for(i in 1:N) {
  vcovs[, , i] <- rinvwishart(nu, diag(rep(s2, d)))
  corrs[, , i] <- cov2cor(vcovs[, , i])
}

sigmas <- apply(vcovs, 3, function(x) sqrt(diag(x))) |> as.vector()
rr <- apply(corrs, 3, function(x) x[lower.tri(x)]) |> as.vector()

par(mfrow = c(1, 2))
plot(density(sigmas, from = 0, to = 100))# hist(sigmas)
hist(rr)
par(mfrow = c(1, 1))


# Use prior nu = d+1, S = I over correlation matrix.
# It makes a Uniform(-1, 1) prior over correlations.


# Correlation matrix updates ----------------------------------------------

eta = 0.1
d <- 6
(R <- rlkjcorr(1, d, eta))
Rvec <- R[lower.tri(R)]
dc <- length(Rvec)

N <- 100
Y <- rmvn(N, rep(0, d), R)

Rsample <- cor(Y)
Rsample_vec <- Rsample[lower.tri(Rsample)]

A <- t(Y) %*% Y

nu_prior <- d + 1 # uniform
S_prior <- diag(rep(1, d))

nu_posterior <- nu_prior + N
S_posterior <- S_prior + t(Y) %*% Y

niter <- 5000
samples_R <- array(NA, dim = c(d, d, niter))
samples_vec <- matrix(NA, dc, niter)

for(i in 1:niter) {
  # samples_R[, , i] <- rinvwishart(nu_posterior, S_posterior) |> cov2cor()
  Ri <- rinvwishart(nu_posterior, S_posterior) |> cov2cor()
  samples_vec[, i] <- Ri[lower.tri(Ri)]
}

par(mfrow = c(3, 5))
for(i in 1:dc) {
  plot(density(samples_vec[i, ], from = -1, to = 1, n = 2^10),
       main = NA, xlab = NA)
  abline(v = Rvec[i], lty = 2, col = 2)
  abline(v = Rsample_vec[i], lty = 2, col = 4)
}
par(mfrow = c(1, 1))
# perfect, works

# Linear regression -------------------------------------------------------

# taking code from
# https://www.r-bloggers.com/2021/05/bayesian-linear-regression-with-gibbs-sampling-using-r-code/

# real values
b_real <- rnorm(2, 4)
s2_real <- rnorm(1, 0, 1) ^ 2
N <- 500
x <- rnorm(N)
X <- cbind(rep(1, N), x)
y <- rnorm(N, X %*% b_real, sqrt(s2_real))

K <- length(b_real)

# initialization for priors
b0     <- rep(0, K)    # priors for B
S0 <- diag(K) * 10  # priors for sigma2 (matrix)
t0 <- 1; d0 <- 0.1     # priors for shape, scale of IG (non informative)
#curve(dinvgamma(x, t0, d0), to = 50, n = 300)

# starting values
b <- b0; s2 <- 1

niter <- 20000

# output (parameters sampled) matrix
samples <- matrix(NA, niter, K+1)
colnames(samples) <- c("beta0", "beta1", "sigma")

# Iterated sampling
for (i in 1:niter) {
  # i = 1
  # print(paste(i, "-th iteration of total", niter, "…….."))

  #———————————————
  # step 2 : Draw B conditional on sigma2
  #———————————————
  V <- solve(solve(S0) + (1 / s2) * (t(X) %*% X))
  M <- V %*% (solve(S0) %*% b0 + (1 / s2) * t(X) %*% y)
  b <- mvrnorm(1, M, V) |> t() |> t()

  #———————————————
  # step 3 : Sample sigma2 conditional on B
  #———————————————
  resids <- y - X %*% b # residuals

  # posterior df and scale matrix
  t1 <- t0 + N
  d1 <- d0 + t(resids) %*% resids

  # Draw sigma2 from IG(T1,D1)
  s2 <- rinvgamma(1, t1, d1)

  # save samples
  samples[i, ] <- c(as.numeric(b), s2)
}

# compare with mle

mm <- lm(y ~ x)

par(mfrow = c(2, 2))
plot(density(samples[-(niter/2), 1]), main = "b0")
abline(v = b_real[1], lty = 2, col = 2)
abline(v = coef(mm)[1], lty = 2, col = 4)

plot(density(samples[-(niter/2), 2]), main = "b1")
abline(v = b_real[2], lty = 2, col = 2)
abline(v = coef(mm)[2], lty = 2, col = 4)

plot(density(samples[-(niter/2), 3], from = 0), main = "s2")
abline(v = s2_real, lty = 2, col = 2)
abline(v = sigma(mm)[1] ^ 2, lty = 2, col = 4)
par(mfrow = c(1, 1))


# Gibbs for truncated normal likelihood -----------------------------------

# From William Griffiths

# Define yut (y untruncated) as the cuantiles of the non-truncated distribution
# having the same cumulative probability as yt (y truncated) in the truncated
# distribution.
# Eq. 12 in Griffiths
# yt: trunctated variable
# a: lower bound
# b: upper bound
untruncate <- function(yt, mu, sigma, a, b) {
  v1 = pnorm((yt - mu) / sigma)
  v2 = pnorm((a - mu) / sigma)
  v3 = pnorm((b - mu) / sigma)
  p = (v1 - v2) / (v3 - v2)
  # avoid inf
  p <- ifelse(p == 1, 1 - 1e-6, p)
  p <- ifelse(p == 0, 1e-6, p)
  yut = mu + sigma * qnorm(p)
  return(yut)
}
# Them, sample with gibbs as if the variable was not truncated.
# (Is it so simple?)

smodel <- stan_model("hierarchical model/truncated_normal.stan")

# real values
b_real <- rnorm(2, 4)#c(0, 0)# rnorm(2, 4)
s2_real <- rnorm(1, 0, 1) ^ 2
N0 <- 1000
x0 <- rnorm(N0)
X0 <- cbind(rep(1, N0), x0)
# y <- rtruncnorm(N, a = 0, b = Inf, X %*% b_real, sqrt(s2_real))
y0 <- rnorm(N0, X0 %*% b_real, sqrt(s2_real))
a <- quantile(y0, prob = 0.75)

# replace dataset
X <- X0[y0>a, ]
x <- x0[y0>a]
y <- y0[y0>a]
N <- length(y)

K <- length(b_real)

# initialization for priors
s0 <- 10
b0 <- rep(0, K)    # priors for B
S0 <- diag(K) * s0^2  # priors for sigma2 (matrix)
t0 <- 0.1; d0 <- 0.1     # priors for shape, rate of IG (non informative)
# curve(dinvgamma(x, t0, d0), to = 50, n = 300)
# sss <- rinvgamma(1e5, t0, d0)
# quantile(sss, prob = 0.95)

mm <- truncreg(y ~ x, point = a, direction = "left")
b <- coef(mm)[1:2]
s2 <- coef(mm)[3]

niter <- 5000

# output (parameters sampled) matrix
samples <- matrix(NA, niter, K+1)
colnames(samples) <- c("beta0", "beta1", "sigma")

# Iterated sampling
for(i in 1:niter) {
  # i = 1
  # print(i)

  #———————————————
  # step 1 : Define yut from b and sigma
  #———————————————
  yut <- untruncate(y, X %*% b, sqrt(s2), a, Inf)
  if(!all(is.finite(yut))) {
    print(paste("non-finte yut"))
  }
  #———————————————
  # step 2 : Draw B conditional on sigma2
  #———————————————
  V <- solve(solve(S0) + (1 / s2) * (t(X) %*% X))
  M <- V %*% (solve(S0) %*% b0 + (1 / s2) * t(X) %*% yut)
  b <- mvrnorm(1, M, V) |> as.numeric()

  #———————————————
  # step 3 : Sample sigma2 conditional on B
  #———————————————
  resids <- yut - X %*% b # residuals

  # posterior df and scale matrix
  t1 <- t0 + N
  d1 <- d0 + t(resids) %*% resids

  # Draw sigma2 from IG(T1,D1)
  s2 <- rinvgamma(1, t1, d1)

  # save samples
  samples[i, ] <- c(b, s2)
}

# check stanfit
sdata <- list(y = y, x = x, N = N, s0 = s0, t0 = t0, d0 = d0, a = a)
mstan <- rstan::sampling(smodel, data = sdata)
samples2 <- as.matrix(mstan)

par(mfrow = c(2, 2))
plot(density(samples[-(niter/2), 1]), main = "b0")
lines(density(samples2[, "b0"]), col = "green")
abline(v = b_real[1], lty = 2, col = 2)
abline(v = coef(mm)[1], lty = 2, col = 4)

plot(density(samples[-(niter/2), 2]), main = "b1")
lines(density(samples2[, "b1"]), col = "green")
abline(v = b_real[2], lty = 2, col = 2)
abline(v = coef(mm)[2], lty = 2, col = 4)

plot(density(samples[-(niter/2), 3], from = 0), main = "s2")
lines(density(samples2[, "sigma2"]), col = "green")
abline(v = s2_real, lty = 2, col = 2)
abline(v = coef(mm)[3]^2, lty = 2, col = 4)
par(mfrow = c(1, 1))

# No, my sampler does not work, it differs from stan quite much. Use metropolis
# updates for area regression.