# Notes -------------------------------------------------------------------

# Most of the code here is inherited from
# "/home/ivan/Insync/Fire spread modelling/fire_spread/multivariate normal models for likelihood emulation/multivariate normal and skew-normal.R"
# The stan model is inherited from
# "/home/ivan/Insync/Fire spread modelling/fire_spread/multivariate normal models for likelihood emulation/multivariate skew-normal - truncated - b0.stan"

# This file, ended at _02 simply edits the previous one.
# Here I use other example database, with a slightly different structure.
# That's why I made another script.

# Tried to fit the SMVN but the sampler is not starting, even when I initialize
# from parameters fitted with sn::selm.


# Packages ----------------------------------------------------------------

library(ggplot2); library(magrittr); theme_set(theme_bw())
library(viridis)
library(sn)
library(rstan)
library(extraDistr)  # half t density (dht)
library(DHARMa)
library(mnormt) # for pd.solve
library(trialr) # rlkjcorr
library(truncnorm)  # check truncated models
library(bayestestR)      # highest density intervals


# Data --------------------------------------------------------------------

dat <- readRDS(file.path("files", "pseudolikelihood_estimation", "smc_waves_2008_3.rds"))

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


# function to get unconstrained xi (position) from the real xi
# (used to provide starting values from a selm fit)
unconstrain_xi <- function(xi, xi_lower, xi_upper, prior_xi_logit_sd) {
  if(any(xi < xi_lower) | any(xi > xi_upper)) stop("xi out of range")
  p <- (xi - xi_lower) / (xi_upper - xi_lower)
  xi_raw <- qlogis(p) / prior_xi_logit_sd
  return(xi_raw)
}

# Prepare data ------------------------------------------------------------

# work at log scale, because steps has a heavily skew distribution,
# and a log-normal would be better.

dat$par_values_log <- dat$par_values
for(i in 2:ncol(dat$par_values)) {
  dat$par_values_log[, i] <- log(dat$par_values[, i])
}

plot(dat$overlap ~ dat$par_values[, "intercept"])
low <- 0.2 # threshold
y_lower <- log(low)
data <- dat[dat$overlap_log > y_lower, ]
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

nchains <- 1
wu <- 500
sa <- 2000

# initialize at reasonable values
init_f_random <- function(chain_id) {
  list(
    xi_logit_raw = rnorm(K, 0, 0.5),
    alpha_raw = rnorm(K, 0, 0.5),
    b0_raw = rnorm(1, 0, 1),
    sigma = rtruncnorm(K, a = 0.01, mean = sds, sd = sds / 2),
    tau = rtruncnorm(1, a = 0, mean = sd(y), sd = sd(y) / 2),
    Lcorr = t(chol(rlkjcorr(1, K, eta = 30)))
  )
}

# initialize at reasonable values
init_f_random <- function(chain_id) {
  list(
    xi_logit_raw = rep(0, K),#rnorm(K, 0, 0.001),
    alpha_raw = rep(0, K),#rnorm(K, 0, 0.001),
    b0_raw = rnorm(1, 0, 0.001),
    sigma = sds,#rtruncnorm(K, a = 0.01, mean = sds, sd = sds / 5),
    tau = sd(y) / 4, #rtruncnorm(1, a = 0, mean = sd(y), sd = sd(y) / 5),
    Lcorr = t(chol(diag(K)))#t(chol(rlkjcorr(1, K, eta = 30)))
  )
}
init_list <- lapply(1:nchains, function(i) init_f_random(i))

## Initialize at parameters from estimated SMVN
ids_boot <- sample(1:nrow(data), size = 5000, prob = data$overlap, replace = TRUE)
data_resample <- data[ids_boot, ]
smvn_brute <- sn::selm.fit(y = data_resample$par_values_log,
                           x = matrix(1, nrow(data_resample)), family = "SN")
dp <- smvn_brute$param$dp
smvn_dist <- sn::makeSECdistr(dp = dp, family = "SN",
                              compNames = colnames(dat$par_values))
# plot(smvn_dist)
# r2 <- sn::rmsn(n = 3e4, dp = f2$param$dp)

# turn fitted parameters to unconstrained scale to use as inits.
init_list <- list(
  xi_logit_raw = unconstrain_xi(as.numeric(dp$beta),
                                sdata$xi_lower, sdata$xi_upper,
                                sdata$prior_xi_logit_sd),
  alpha_raw = dp$alpha / sdata$prior_alpha_sd,
  b0_raw = ((max(y) - sdata$prior_b0_mean) / sdata$prior_b0_sd) * 0.5,
  sigma = sqrt(diag(dp$Omega)),
  tau = sd(y) / 6,
  Lcorr = t(chol(cov2cor(dp$Omega)))
)

m <- sampling(smodel, data = sdata, refresh = 20,
              cores = nchains,
              chains = nchains,
              iter = 10,
              warmup = 5,
              init = list(init_list),# init_f_random, #init_f_true,
              control = list(max_treedepth = 25, adapt_delta = 0.9))
# No arranca el muestreo, está re difícil.


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




# MVN, porque skew no funciona --------------------------------------------

# compile stan model
smodel_mvn <- stan_model(
  file.path("pseudolikelihood estimation",
            "pseudolikelihood_estimation_mvn_01.stan"), # truncated, mvn (symmetrical)
  verbose = F
)

# parameters will be heavily constrained for the surface to be fitted.
sdata_mvn <- list(
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

  # residual standard deviation (student_t)
  prior_tau_sd = sd(y) * 2,
  prior_tau_nu = 3,

  # prior eta parameter for correlation matrix
  prior_corr_eta = 1
)

nchains <- 4
wu <- 500
sa <- 2000

# initialize at reasonable values (random no arrancaba)
init_f_random <- function(chain_id) {
  list(
    xi_logit_raw = rep(0, K),#rnorm(K, 0, 0.001),
    alpha_raw = rep(0, K),#rnorm(K, 0, 0.001),
    b0_raw = rnorm(1, 0, 0.001),
    sigma = sds*10,#rtruncnorm(K, a = 0.01, mean = sds, sd = sds / 5),
    tau = sd(y), #rtruncnorm(1, a = 0, mean = sd(y), sd = sd(y) / 5),
    Lcorr = t(chol(rlkjcorr(1, K, eta = 30)))
  )
}
init_list <- lapply(1:nchains, function(i) init_f_random(i))

m <- sampling(smodel_mvn, data = sdata_mvn, refresh = 50,
              cores = nchains,
              chains = nchains,
              iter = sa,
              warmup = wu,
              init = init_list,# init_f_random, #init_f_true,
              control = list(max_treedepth = 25, adapt_delta = 0.9))

# Con MVN inicia... pero no funcionó. No converge, todas divergent transitions.


