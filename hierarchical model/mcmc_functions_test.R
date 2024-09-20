source("hierarchical model/mcmc_functions")

# Stan code ---------------------------------------------------------------

stan_trunc <- stan_model("hierarchical model/truncated_normal.stan")
stan_jacob <- stan_model("hierarchical model/jacobian.stan")

# update_lm ---------------------------------------------------------------

# taking code from
# https://www.r-bloggers.com/2021/05/bayesian-linear-regression-with-gibbs-sampling-using-r-code/

# real values
b_real <- rnorm(2, 4)
s2_real <- rnorm(1, 0, 1) ^ 2
N <- 500
x <- rnorm(N)
X <- cbind(rep(1, N), x)
tXX <- t(X) %*% X
y <- rnorm(N, X %*% b_real, sqrt(s2_real))

K <- length(b_real)

# initialization for priors
b0 <- rep(0, K)       # priors for B
S0 <- diag(K) * 10    # priors for sigma2 (matrix)
S0_inv <- solve(S0)
t0 <- 1; d0 <- 0.1    # priors for shape, scale of IG (non informative)
#curve(dinvgamma(x, t0, d0), to = 50, n = 300)

# starting values
b <- b0; s2 <- 1

niter <- 20000

# output (parameters sampled) matrix
samples <- matrix(NA, niter, K+1)
colnames(samples) <- c("beta0", "beta1", "sigma")

# initialize
samples[1, 1:K] <- b
samples[1, 3] <- s2

# sample
for(i in 2:niter) {
  cc <- update_lm(y = y, X = X, tXX = tXX, s2 = samples[i-1, 3],
                  b0 = b0, S0_inv = S0_inv,
                  t0 = t0, d0 = d0)
  # save samples
  samples[i, ] <- cc
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

# update_corr -------------------------------------------------------------

eta <- 0.1
d <- 6
(R <- rlkjcorr(1, d, eta))
Rvec <- R[lower.tri(R)]
dc <- length(Rvec)

N <- 100
Y <- mvrnorm(N, rep(0, d), R)

Rsample <- cor(Y)
Rsample_vec <- Rsample[lower.tri(Rsample)]

niter <- 5000
samples_R <- array(NA, dim = c(d, d, niter))
samples_vec <- matrix(NA, dc, niter)

for(i in 1:niter) {
  Ri <- update_corr(Y)
  samples_vec[, i] <- Ri[lower.tri(Ri)]
}

# plot
par(mfrow = c(3, 5))
for(i in 1:dc) {
  plot(density(samples_vec[i, ], from = -1, to = 1, n = 2^10),
       main = NA, xlab = NA)
  abline(v = Rvec[i], lty = 1, col = 2)
  abline(v = Rsample_vec[i], lty = 2, col = 4)
}
par(mfrow = c(1, 1))
# perfect, works


# update_truncnorm --------------------------------------------------------

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

# priors
s0 <- 10
b0 <- rep(0, K)        # priors for B
S0 <- diag(K) * s0^2   # priors for sigma2 (matrix)
t0 <- 0.1; d0 <- 0.1   # priors for shape, rate of IG (non informative)
# curve(dinvgamma(x, t0, d0), to = 50, n = 300)
# sss <- rinvgamma(1e5, t0, d0)
# quantile(sss, prob = 0.95)

# fit model with stan, to get sd for proposals
sdata <- list(y = y, x = x, N = N, s0 = s0, t0 = t0, d0 = d0, a = a)
mstan <- rstan::sampling(stan_trunc, data = sdata)
samples1 <- as.matrix(mstan, pars = c("b0", "b1", "sigma2"))
inits <- colMeans(samples1)
samples2 <- as.matrix(mstan, pars = c("b0", "b1", "sigma"))
sd_jump <- apply(samples2, 2, sd) * 3

# also compare with MLE
mm <- truncreg(y ~ x, point = a, direction = "left")

# output (parameters sampled) matrix
niter <- 20000
samples <- matrix(NA, niter, K+1)
colnames(samples) <- c("beta0", "beta1", "sigma2")
samples[1, ] <- inits

for(i in 2:niter) {
  cc <- update_truncnorm(y, x,
                         coef = samples[i-1, ],
                         L = a, b0, S0, t0, d0,
                         sd_jump)
  samples[i, ] <- cc
}

par(mfrow = c(2, 2))
plot(density(samples[-(niter/2), 1]), main = "b0")
lines(density(samples1[, "b0"]), col = "green")
abline(v = b_real[1], lty = 1, col = 2)
abline(v = coef(mm)[1], lty = 2, col = 4)

plot(density(samples[-(niter/2), 2]), main = "b1")
lines(density(samples1[, "b1"]), col = "green")
abline(v = b_real[2], lty = 1, col = 2)
abline(v = coef(mm)[2], lty = 2, col = 4)

plot(density(samples[-(niter/2), 3], from = 0), main = "s2")
lines(density(samples1[, "sigma2"]), col = "green")
abline(v = s2_real, lty = 2, col = 2)
abline(v = coef(mm)[3]^2, lty = 1, col = 4)
par(mfrow = c(1, 1))


# jacobian adjustment for exp_logit transform -----------------------------

U <- 90
L <- 60
mu <- 1
sigma <- 3
N <- 1e5

steps_logit <- rnorm(N, mu, sigma)
steps <- plogis(steps_logit) * (U - L) + L # in [10, 11]
steps_log <- log(steps)

initf <- function() {
  return(list(steps_log = runif(1, L, U) |> log()))
}

sdata <- list(L = L, U = U, mu = mu, sigma = sigma)

m1 <- sampling(stan_jacob, data = sdata, init = initf,
               warmup = 10000, iter = 20000, chains = 4, cores = 4,
               control = list(adapt_delta = 0.999))

xx <- as.matrix(m1)

# Transformation of samples vs stan sampling
par(mfrow = c(1, 3))
plot(density(steps_logit), main = "steps_logit")
lines(density(xx[, "steps_logit"]), col = "blue")
lines(density(xx[, "steps_logit2"]), col = "red", lty = 2)

plot(density(steps, from = L, to = U),
     xlim = c(L, U), main = "steps")
lines(density(xx[, "steps"], from = L, to = U),
      col = "blue")

plot(density(steps_log, from = log(L), to = log(U)),
     xlim = log(c(L, U)), main = "steps_log")
lines(density(xx[, "steps_log"], from = log(L), to = log(U)),
      col = "blue")
par(mfrow = c(1, 1))

# Samples densities vs. computed densities

steps_log_seq <- seq(min(steps_log), max(steps_log), length.out = 1000)
steps_log_d <-
  dnorm(exp_logit(steps_log_seq, L, U), mu, sigma) *
  abs(exp_logit_d(steps_log_seq, L, U))

par(mfrow = c(1, 2))
plot(density(steps_logit), main = "steps_logit", lwd = 2)
curve(dnorm(x, mu, sigma), add = T, col = "red")

plot(steps_log_d ~ steps_log_seq, type = "l", main = "steps_log", lwd = 2)
lines(density(steps_log, adjust = 0.1), col = "red")
par(mfrow = c(1, 1))

# perfect.