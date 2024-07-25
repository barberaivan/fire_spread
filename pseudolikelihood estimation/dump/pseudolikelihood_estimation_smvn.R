
# Notes -------------------------------------------------------------------

# Most of the code here is inherited from
# "/home/ivan/Insync/Fire spread modelling/fire_spread/multivariate normal models for likelihood emulation/multivariate normal and skew-normal.R"
# The stan model is inherited from
# "/home/ivan/Insync/Fire spread modelling/fire_spread/multivariate normal models for likelihood emulation/multivariate skew-normal - truncated - b0.stan"


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

wdata <- readRDS(file.path("files", "wave_test.rds"))

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



# Prepare data ------------------------------------------------------------

dat <- wdata$w6$particles_data
low <- max(wdata$similarity_bounds[, "ll"])
y_lower <- qlogis(plogis(low) + 0.05)
dat <- dat[dat$ll > y_lower, ]
dat <- dat[dat$in_bounds == 1 & dat$ll > y_lower, ]
nrow(dat) # 212 may be too little

hist(dat$ll %>% plogis) # terrible

# reasonable values for sigma
x_ranges <- sapply(1:ncol(dat$par_values), function(c) {
  diff(range(dat$par_values[, c]))
})

sds <- apply(dat$par_values, 2, sd)

# data out
X = dat$par_values
y = dat$ll
N = nrow(dat)
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
  "pseudolikelihood_estimation_smvn_01.stan", # truncated, skew, mvn
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
    tau = sd(y), #rtruncnorm(1, a = 0, mean = sd(y), sd = sd(y) / 5),
    Lcorr = t(chol(rlkjcorr(1, K, eta = 30)))
  )
}
init_list <- lapply(1:nchains, function(i) init_f_random(i))

m <- sampling(smodel, data = sdata, refresh = 20,
              cores = nchains,
              chains = nchains,
              iter = 10,
              warmup = 5,
              init = init_list,# init_f_random, #init_f_true,
              control = list(max_treedepth = 25, adapt_delta = 0.9))



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

dat <- wdata$w6$particles_data
low <- max(wdata$similarity_bounds[, "ll"])
y_lower <- qlogis(plogis(low) + 0.05)
dat <- dat[dat$ll > y_lower, ]
dat <- dat[dat$in_bounds == 1 & dat$ll > y_lower, ]
nrow(dat) # 212 may be too little

hist(dat$ll %>% plogis) # terrible

# reasonable values for sigma
x_ranges <- sapply(1:ncol(dat$par_values), function(c) {
  diff(range(dat$par_values[, c]))
})

sds <- apply(dat$par_values, 2, sd)

# data out
X = dat$par_values
y = dat$ll
N = nrow(dat)
K = ncol(X)

# compile stan model
smodel_mvn <- stan_model(
  "pseudolikelihood_estimation_mvn_01.stan", # truncated, mvn (symmetrical)
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
pairs(m, pars = c("xi"))
# todo como el culo. todas transiciones divergentes, re mal mixeo.


# MVN sin truncar? quitando los out of bounds y dejando el resto ----------

dat <- wdata$w6$particles_data
dat <- dat[dat$in_bounds == 1, ]
nrow(dat) # 2821 if we use all

hist(dat$ll %>% plogis) # terrible

# reasonable values for sigma
x_ranges <- sapply(1:ncol(dat$par_values), function(c) {
  diff(range(dat$par_values[, c]))
})

sds <- apply(dat$par_values, 2, sd)

# data out
X = dat$par_values
y = dat$ll
N = nrow(dat)
K = ncol(X)

# compile stan model
smodel_mvn <- stan_model(
  "pseudolikelihood_estimation_mvn_02.stan", # truncated, mvn (symmetrical)
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

# initialize at reasonable values
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
# 1572.82 / 60 = 26 min
# saveRDS(m, "llmodel_trials_mvn_all_in_bounds.rds")
pairs(m, pars = c("xi"))
# solo tuvo una divergent transition.
# pero da ocote.
# estima un intercept re alto, al borde de lo que le permito.



# Simplest: just a quadratic model, no correlation ------------------------

# compile stan model
smodel_quad <- stan_model(
  "pseudolikelihood_estimation_quad_01.stan", # truncated, mvn (symmetrical)
  verbose = F
)


dat <- wdata$w6$particles_data
low <- max(wdata$similarity_bounds[, "ll"])
y_lower <- qlogis(plogis(low) + 0.05)
dat <- dat[dat$in_bounds == 1, ]# & dat$ll > y_lower, ]
nrow(dat) # 2821 if we use all

hist(dat$ll %>% plogis) # terrible

# reasonable values for sigma
x_ranges <- sapply(1:ncol(dat$par_values), function(c) {
  diff(range(dat$par_values[, c]))
})

sds <- apply(dat$par_values, 2, sd)
precs <- 1 / sds ^ 2

# data out
X = dat$par_values
y = dat$ll
N = nrow(dat)
K = ncol(X)

# parameters will be heavily constrained for the surface to be fitted.
sdata_quad <- list(
  # data
  X = X, y = y, N = N, K = K,

  # lower truncation for the response
  L = y_lower,
  truncated = 0,

  # force b0 to be near the optimum
  prior_b0_mean = max(y),
  prior_b0_sd = diff(range(y)), # sd(y) * 4

  # constrain xi (~mu) within the 90 % HDI of the observations
  xi_lower = apply(X, 2, hdi_lower),
  xi_upper = apply(X, 2, hdi_upper),
  prior_xi_logit_sd = 1.5,

  # prior for sigma (student_t)
  prior_prec_sd = precs * 4,
  prior_prec_nu = 3,

  # residual standard deviation (student_t)
  prior_tau_sd = sd(y) * 2,
  prior_tau_nu = 3,

  # prior eta parameter for correlation matrix
  prior_corr_eta = 1
)

nchains <- 4
wu <- 500
sa <- 2000

# initialize at reasonable values
init_f_random <- function(chain_id) {
  list(
    xi_logit_raw = rep(0, K),#rnorm(K, 0, 0.001),
    b0_raw = rnorm(1, 0, 0.001),
    precision = 1 / (sds*10) ^ 2,#rtruncnorm(K, a = 0.01, mean = sds, sd = sds / 5),
    tau = sd(y) #rtruncnorm(1, a = 0, mean = sd(y), sd = sd(y) / 5),
  )
}
init_list <- lapply(1:nchains, function(i) init_f_random(i))

mq <- sampling(smodel_quad, data = sdata_quad, refresh = 200,
              cores = nchains,
              chains = nchains,
              iter = sa,
              warmup = wu,
              init = init_list,# init_f_random, #init_f_true,
              control = list(max_treedepth = 25, adapt_delta = 0.9))
# muestrea rapido, pero...
# no converge ni a palos (usando el truncado con in_bounds)

# usando el no truncado con todos los datos in_bounds...
# a pesar del alto N, vuela comparado con el MVN (este no tiene correlación)

pairs(mq, pars = c("xi"))
pairs(mq, pars = c("precision"))
# estima precisiones bajísimas.
# ahora al menos anda, estima cosas sin drama, pero sigue estimando valores
# absurdos para le intercept
1/sqrt(0.04)

plot(plogis(dat$ll[dat$in_bounds == 1]) ~ dat$par_values[dat$in_bounds == 1, "intercept"])

# Fixing mu at the best iteration -----------------------------------------

smodel_quad <- stan_model(
  "pseudolikelihood_estimation_quad_02.stan", # truncated, mvn (symmetrical)
  verbose = F
)

dat <- wdata$w6$particles_data
low <- max(wdata$similarity_bounds[, "ll"])
y_lower <- qlogis(plogis(low) + 0.05)
# dat <- dat[dat$in_bounds == 1, ]# & dat$ll > y_lower, ]
dat <- dat[dat$ll > y_lower, ]
nrow(dat) # 2821 if we use all; 781 si usamos solo los que superan ese umbral de y_lower

xi_fixed <- dat$par_values[which.max(dat$ll), ] #
xi_fixed <- wdata$w6$loglik_optim$par # max fitted by GP

# fit a quadratic using all the in_bounds
# compile stan model

# reasonable values for sigma
x_ranges <- sapply(1:ncol(dat$par_values), function(c) {
  diff(range(dat$par_values[, c]))
})

sds <- apply(dat$par_values, 2, sd)
precs <- 1 / sds ^ 2

# data out
X = dat$par_values
y = dat$ll
N = nrow(dat)
K = ncol(X)

# parameters will be heavily constrained for the surface to be fitted.
sdata_quad <- list(
  # data
  X = X, y = y, N = N, K = K,

  # lower truncation for the response
  L = y_lower,
  truncated = 1, # OJO CON LA TRUNCATION!!

  xi = xi_fixed,

  # force b0 to be near the optimum
  prior_b0_mean = max(y),
  prior_b0_sd = diff(range(y)), # sd(y) * 4

  # constrain xi (~mu) within the 90 % HDI of the observations
  xi_lower = apply(X, 2, hdi_lower),
  xi_upper = apply(X, 2, hdi_upper),
  prior_xi_logit_sd = 1.5,

  # prior for sigma (student_t)
  prior_prec_sd = rep(10, K),#precs * 4,
  prior_prec_nu = 3,

  # residual standard deviation (student_t)
  prior_tau_sd = sd(y) * 2,
  prior_tau_nu = 3,

  # prior eta parameter for correlation matrix
  prior_corr_eta = 1
)

nchains <- 4
wu <- 500
sa <- 2000

# initialize at reasonable values
init_f_random <- function(chain_id) {
  list(
    b0_raw = rnorm(1, 0, 0.001),
    precision = 1 / (sds*10) ^ 2,#rtruncnorm(K, a = 0.01, mean = sds, sd = sds / 5),
    tau = sd(y) #rtruncnorm(1, a = 0, mean = sd(y), sd = sd(y) / 5),
  )
}
init_list <- lapply(1:nchains, function(i) init_f_random(i))

mqf <- sampling(smodel_quad, data = sdata_quad, refresh = 200,
               cores = nchains,
               chains = nchains,
               iter = 3000,
               warmup = 2000,
               init = init_list,# init_f_random, #init_f_true,
               control = list(max_treedepth = 25, adapt_delta = 0.98))

## USANDO TODOS LOS IN-BOUNDS.
# este es rapido. por ahi podemos animarnos a meterle las correlaciones.
# anduvooooooooo. 153 s
pairs(mqf, pars = "precision")
# saveRDS(mqf, "llmodel_trials_quad_all_in_bounds_xi_fixed.rds")

# sigma en mu de ~ 7
ss <- summary(mqf)[[1]]
View(ss)

b0_hat <- -3.6
prec1_hat <- 0.0183643353


# plot
plot(plogis(dat$ll[dat$in_bounds == 1]) ~ dat$par_values[dat$in_bounds == 1, "intercept"],
     col = rgb(0, 0, 0, 0.1), pch = 19)
curve(plogis(b0_hat - prec1_hat * (x - xi_fixed[1]) ^ 2), add = TRUE)


# dfull
df <- wdata$w6$particles_data

plot(df$ll[df$in_bounds == 0] ~ df$par_values[df$in_bounds == 0, "intercept"],
     col = rgb(1, 0, 0, 0.2), pch = 19)
points(df$ll[df$in_bounds == 1] ~ df$par_values[df$in_bounds == 1, "intercept"],
     col = rgb(0, 0, 0, 0.1), pch = 19)
curve(b0_hat - prec1_hat * (x - xi_fixed[1]) ^ 2, add = TRUE, lwd = 2)
## -------------


# ahora comento su uso todos los ll > y_lower
# va a los pedos, truncada.
# 233 div transitions.
# tarda 50 s, pero las transiciones persisten.

# ahora le remil agrandé la previa a las precisiones, a ver si encuentra valores más grandes.

pairs(mqf, pars = "precision") # estima parecido.

ss <- summary(mqf)[[1]]
View(ss)
b0_hat <- -2.1
prec1_hat <- 0.0080849577

df <- wdata$w6$particles_data

plot(df$ll[df$in_bounds == 0] ~ df$par_values[df$in_bounds == 0, "intercept"],
     col = rgb(1, 0, 0, 0.2), pch = 19)
points(df$ll[df$in_bounds == 1] ~ df$par_values[df$in_bounds == 1, "intercept"],
       col = rgb(0, 0, 0, 0.1), pch = 19)
curve(b0_hat - prec1_hat * (x - xi_fixed[1]) ^ 2, add = TRUE, lwd = 2)

# overlap scale
plot(plogis(df$ll[df$in_bounds == 0]) ~ df$par_values[df$in_bounds == 0, "intercept"],
     col = rgb(1, 0, 0, 0.2), pch = 19)
points(plogis(df$ll[df$in_bounds == 1]) ~ df$par_values[df$in_bounds == 1, "intercept"],
       col = rgb(0, 0, 0, 0.1), pch = 19)
curve(plogis(b0_hat - prec1_hat * (x - xi_fixed[1]) ^ 2), add = TRUE, lwd = 2)

b0s <- as.matrix(mqf, pars = "b0") %>% as.numeric
p1s <- as.matrix(mqf, pars = "precision[1]") %>% as.numeric

for(i in 1:length(b0s)) {
  curve(plogis(b0s[i] - p1s[i] * (x - xi_fixed[1]) ^ 2), add = TRUE,
        col = rgb(0, 0, 1, 0.03))
}



# NOTES -------------------------------------------------------------------

# maybe it's better to parameterize sigma through precision, so priors near zero
# make a flat curve, and not a terribly steep one. Moreover, as P ~ 0, sigma
# would need to go to infinity, which could happen and wouldn't make sense.

# anyway, the quadratic fits terribly.

# If we fit a good GP at the last wave, its maximum is a good candidate. Then,
# we can fix the mode of the MVN there, and fit the curve using all the in_bounds
# data.

# está re lleno de out_bounds en valores de ll muy altos. eso está raro.
# revisar las waves. Quizás el problema sea la definicion de in_ por varianza.


# ver si este problema ocurre con otros fuegos. O armar un modelo de eff fijos y ya?
# haciendo eso quizás esté a tiempo de mostrar algo para la RAE.

# parece que este es un fuego poco informativo, por ahi valga la pena probar
# con otros.

# parece que no voy a llegar. mejor mañana hacer la simulación con todos los fuegos
# a la vez, estimando modelo de eff fijo.