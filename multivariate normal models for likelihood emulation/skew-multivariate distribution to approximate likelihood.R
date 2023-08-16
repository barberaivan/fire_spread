# For context, read the "notes on estimation" file.
# Here I explore the skew-multivariate-(normal/t) distributions to approximate
# the likelihood.


# Packages ----------------------------------------------------------------

library(sn)
library(Matrix)     # to find nearest positive-definite matrix
library(matrixcalc) # test positive-definiteness (is.positive.definite)
library(tidyverse); theme_set(theme_bw())
library(viridis)
library(mgcv)       # compare with a gam
library(scam)       # shape-constrained gam (convex)
library(GA)         # global optimization
# see sn::pkg_sn-intro, many example codes.


# Exploring a 3D t --------------------------------------------------------

dp3 <- list(xi = c(3, -2, 0),                     # location
            Omega = diag(1:3) + outer(1:3,1:3)/2, # scale (~vcov)
            alpha = c(5, -2, 6),                  # slant 
            nu = 6)                               # df, must be constant
st3 <- makeSECdistr(dp = dp3, family = "ST", name = "Multiv.ST",
                    compNames = c("X", "Y", "Z"))
show(st3) # of just type ’st3’
plot(st3)
# compute log-density at the mode
mode <- modeSECdistr(object = st3) 
dmst(x = mode, dp = dp3, log = TRUE)

dmst(x = mode, log = TRUE, 
     dp = list(
       xi = c(3, -2, 0),                     # location
       Omega = diag(1:3) + outer(1:3,1:3)/2, # scale (~vcov)
       alpha = c(5, -2, 6),                  # slant 
       nu = 3)
     )

# get marginal distribution
marg3 <- marginalSECdistr(object = st3, comp = 3)
# plot(marg3)



# Plotting a univariate t ------------------------------------------------

dp_uni_1 <- c(xi = 0,     # location
            omega = 1,  # scale 
            alpha = 2, # slant 
            nu = Inf)     # df

dp_uni_2 <- c(xi = 0,     # location
              omega = 1,  # scale 
              alpha = 0, # slant 
              nu = Inf)     # df


t_uni_1 <- makeSECdistr(dp = dp_uni_1, family = "ST")
t_uni_2 <- makeSECdistr(dp = dp_uni_2, family = "ST")

# curve(dst(x, dp = dp_uni_1), from = xl[1], to = xl[2])

# par(mfrow = c(1, 2))
# 
# plot(t_uni_1)
# abline(v = modeSECdistr(object = t_uni_1), col = "blue")
# abline(v = mean(t_uni_1), col = "red")
# 
# plot(t_uni_2)
# abline(v = modeSECdistr(object = t_uni_2), col = "blue")
# abline(v = mean(t_uni_2), col = "red")
# 
# par(mfrow = c(1, 1))

# location is somewhat like the symmetry axis;
# scale is just the scale (it doesn't affect the shape)
# con nu ~15 o 20 ya es muy parecida a una normal.
# alpha (slant) ~ 1 es bastante simétrico. En 5 ya es muy asimétrico
# 2 es relativamente bajo, poca asimetría. quizás 1 es un buen punto de partida.


# functions ---------------------------------------------------------------

# link function for correlation parameters. They scale from (-Inf, Inf) to 
# (-1, 1) and vice versa. They are like logit and inv_logit but scaling the 
# "probability" to [-1, 1]
inv_logit_cor <- function(x) plogis(x) * 2 - 1
logit_cor <- function(x) qlogis((x + 1) / 2)

# function to compute the number of correlation parameters as a function of the
# multivariate distribution dimension
cor_num <- function(d) (d ^ 2 - d) / 2

# inverse function: detect dimension of the matrix from the length of the 1D 
# vector of correlation coefficients
cor_dim <- function(l) 1 / 2 * (sqrt(8 * l + 1) + 1)
# using wolframalpha: https://www.wolframalpha.com/widgets/view.jsp?id=c86d8aea1b6e9c6a9503a2cecea55b13

# function to compute the total number parameters to define a multivariate 
# skew-t distribution as a function of is dimension
coef_num <- function(d) d * 3 + 1 + cor_num(d)

# positive-definiteness check made by makeSECdistr {sn}. It differs from 
# is.positive.definite(), from {matrixcalc}
is_positive_definite <- function(m) {
  min(eigen(m, symmetric = TRUE, only.values = TRUE)$values) > 0
}

is_positive_semi_definite <- function(m) {
  min(eigen(m, symmetric = TRUE, only.values = TRUE)$values) >= 0
}

# make_positive_definite makes a vcov positive definite. It's supposed that if
# the correlation matrix is pd, then the vcov contructed with any set of 
# scales would be. But it isn't so in practice.
# https://stats.stackexchange.com/questions/23317/create-positive-definite-3x3-covariance-matrix-given-specified-correlation-value

# If the positive-definiteness fails after max_iter trials, the result is NA.
# This is useful for optim(), because otherwise the iteration number for finding
# a positive-definite matrix took too long. It's safer to just ignore problematic
# parameter vector
make_positive_definite <- function(m, corr, max_iter = 10e3) {
  if(is_positive_definite(m)) {
    return(m)
  } else {
    temp <- nearPD(m, corr = corr, keepDiag = TRUE, maxit = max_iter)
    mat <- as.matrix(temp$mat)
    if(!is_positive_definite(mat)) {
      return(NA)
    } else return(mat)
  }
}

# make pd by adding a small value to the diagonal. (COMPLETE LATER)
make_positive_definite <- function(m, corr, max_iter = 10e3) {
  if(is_positive_definite(m)) {
    return(m)
  } else {
    temp <- nearPD(m, corr = corr, keepDiag = TRUE, maxit = max_iter)
    mat <- as.matrix(temp$mat)
    if(!is_positive_definite(mat)) {
      return(NA)
    } else return(mat)
  }
}

# Function to make positive-definite vcov matrix. The positive-definiteness
# correction is applied first at the corr matrix and then, if necessary, at the 
# vcov matrix. Theoretically, the second stage shouldn't be needed, but, in
# practice, it was. Forcing the positive-definiteness only at the vcov level
# led to frequent failures. 
make_vcov <- function(rho_vector, sigma) {
  
  d <- cor_dim(length(rho_vector))
  
  Rho <- matrix(NA, d, d)
  diag(Rho) <- 1
  Rho[lower.tri(Rho)] <- rho_vector
  Rho <- as.matrix(forceSymmetric(Rho, "L"))
  Rho <- make_positive_definite(Rho, corr = TRUE)
  if(anyNA(Rho)) return(Rho)
  
  else {
    Omega <- as.matrix((sigma %*% t(sigma)) * Rho)
    Omega <- make_positive_definite(Omega, corr = FALSE)
    return(Omega) # may be NA if is not positive definite after many iterations
  }
}

# make vcov not forcing positive-definiteness
make_vcov_raw <- function(rho_vector, sigma) {
  
  d <- cor_dim(length(rho_vector))
  
  Rho <- matrix(NA, d, d)
  diag(Rho) <- 1
  Rho[lower.tri(Rho)] <- rho_vector
  Rho <- as.matrix(forceSymmetric(Rho, "L"))
  Omega <- as.matrix((sigma %*% t(sigma)) * Rho)
  return(Omega)
}


# Function to test whether the make_vcov function always returns a pd matrix.
# to be used, the make_vcov function should be twicked so it returns
# T or F according to pdness.
pd_test <- function() {
  d <- sample(3:20, 1)
  rho <- runif(cor_num(d), -1, 1)
  sigma <- rexp(d)
  return(make_vcov(rho, sigma))
}
# sum(replicate(1e3, pd_test(), simplify = T))

# coef_to_dp converts the raw model coefficients (all at the unconstrained scale) 
# to the dp parameterization of the multivariate skew-t distribution.
# dp stands for "direct parameterization", used in the sn package.
# coef is the unconstrained parameters vector that will be passed to optim(),
# excluding the intercept ("alpha") and the residual standard deviation 
# ("tau").
# d is the dimension of the multivariate density function to fit.
coef_to_dp <- function(coef, d) {
  
  # coef <- rnorm(coef_num(9)); d <- 9 # test
  
  # number of correlation parameters
  n_cor <- cor_num(d)
  
  # get positions (ordered as xi, rho, sigma, alpha, nu)
  xi_ids <- 1:d
  rho_ids <- (xi_ids[d] + 1) : (xi_ids[d] + n_cor)
  sigma_ids <- (rho_ids[length(rho_ids)] + 1) : (rho_ids[length(rho_ids)] + d)
  alpha_ids <- (sigma_ids[d] + 1) : (sigma_ids[d] + d)
  nu_id <- length(coef)
  
  # extract correlation parameters vector
  rho <- inv_logit_cor(coef[rho_ids])  # logit_cor link
  
  # extract sigma vector
  sigma <- exp(coef[sigma_ids]) # log link for sigma (sds)
  
  # make vcov matrix forcing positive-definiteness after many but not infinite
  # iterations of near PD
  Omega <- make_vcov(rho, sigma) 
  if(anyNA(Omega)) warning("Omega is not positive definite.")
  
  # alpha is not constrained
  
  # df
  nu <- exp(coef[nu_id])  # log link for nu
  
  dp <- list(
    xi = coef[xi_ids],
    Omega = Omega, 
    alpha = coef[alpha_ids],
    nu = nu
  )
  
  return(dp)
}

# dp_to_coef makes the inverse transformation. It is useful to set initial values
# for optim at the scale of interest. It may use the parameters (xi, rho, sigma,
# alpha, nu) to create the coef vector (unconstrained) or to create a dp object.
# parameters should be provided in their untransformed scale.
dp_to_coef <- function(dp) {
  
  # # TEST
  # dp <- dpp
  
  # Get rho
  Rho <- cov2cor(dp$Omega)
  rho <- Rho[lower.tri(Rho)]
  
  # Get sigma
  sigma <- sqrt(diag(dp$Omega))
  
  coef <- c(
    dp$xi,
    logit_cor(rho),
    log(sigma),
    dp$alpha,
    log(dp$nu)
  )
  
  return(coef)
}

# params_to_dp creates a dp list from the parameter values in their original 
# scale (constrained when appropriate), but using sigma and rho values to create
# the Omega matrix
params_to_dp <- function(xi, rho, sigma, alpha, nu)  {
  
  # # TEST
  # d = 3
  # xi = rnorm(d)
  # rho = runif(cor_num(d), -1, 1)
  # sigma = rexp(d)
  # alpha = rnorm(d)
  # nu = rexp(1)

  Omega <- make_vcov(rho, sigma)
  if(anyNA(Omega)) warning("Omega is not positive definite.")
  
  # create list
  dp <- list(
    xi = xi, 
    Omega = Omega,
    alpha = alpha, 
    nu = nu
  )
  
  return(dp)
}

# Function to compute the mean of the mst-function. b0 is the maximum value 
# of the function, as the raw mst curve (at log scale) has maximum at zero.
mst_mu <- function(dp, b0, X) {
  # scale maximum log-density to zero 
  mode <- modeSECdistr(dp, family = "ST")
  ld_max <- dmst(x = mode, dp = dp, log = TRUE)
  # dmst(x = mode, dp = dp, log = TRUE) - ld_max # OK, should be 0
  
  # Compute mean
  mu <- b0 + dmst(x = X, dp = dp, log = TRUE) - ld_max
  
  return(mu)
} 


# Function to fit multivariate skew-t function (mst) to a simulated-likelihood
# dataset. This function is to be passed to optim using a wrapper.
# coef_extended is the unconstrained parameter vector to evaluate. The last 2 
# are the intercept and the (log) residual standard deviation, corresponding to 
# a Normal data model.
# X is a matrix with the fire-spread-parameters values, and y is the observed
# similarity between simulated and observed fires.
mst_loglik <- function(coef_extended, X, y) {

  # separate parameters
  last <- length(coef_extended)
  dp <- coef_to_dp(coef_extended[1:(last - 2)], d = ncol(X)) # remove b0 and tau
  if(anyNA(dp$Omega)) return(NA)
  b0 <- coef_extended[last - 1]
  tau <- exp(coef_extended[last])
  
  # compute mu
  mu <- mst_mu(dp = dp, b0 = b0, X = X)
  
  # Compute -log-likelihood for the simulated likelihood. Negate to find
  # minimum (default optim() behaviour)
  ll <- sum(dnorm(y, mu, tau, log = TRUE)) * (-1)
  
  return(ll)
}

# Function to fit mst model.
# X is the predictors matrix, y is the response variable,
# init is the initial value for the parameters, and 
# max_iter is the maximum number of iterations for optim(), and trace
# the detail level for printing the progress (see ?optim).
# value: list with many objects...
mst_fit <- function(init, X, y, max_iter = 1e5, trace = 1) {
  
  # wrapper for optim()
  mst_loglik_wrapper <- function(params) {
    mst_loglik(params, X = X, y = y)
  }
  
  # fit model and record runtime
  mbm = microbenchmark::microbenchmark(
    test = {fit <- optim(init, mst_loglik_wrapper, 
                         control = list(maxit = max_iter, trace = trace))},
    times = 1
  )

  # tidy parameters
  coef_mle <- fit$par
  last <- length(coef_mle)
  
  dp <- coef_to_dp(coef_mle[1:(last - 2)], d = ncol(X))
  b0 <- coef_mle[last - 1]
  tau <- exp(coef_mle[last])
  
  # compute mode and ld_max 
  mode <- modeSECdistr(dp, "ST")
  ld_max <- dmst(mode, dp = dp, log = TRUE)
  
  # tidy result
  out <- list(
    optim_result = fit,
    coef_mle = fit$par,
    params_mle = list(dp = dp, b0 = b0, tau = tau, 
                      ld_max = ld_max, mode = mode),
    y = y,
    X = X,
    runtime = mbm
  )
  
  if(fit$convergence > 0) print("WARNING: Optimization did not converge.")
  
  return(out)
}


# Same function but with another method, usgin gradient.
# The wrapper doesn't return NA or Inf, but a high_value (to be setted)
mst_fit_grad <- function(init, X, y, max_iter = 1e5, trace = 1, method = "BFGS",
                         high_value = 1e5) {
  
  # wrapper for optim()
  mst_loglik_wrapper <- function(params) {
    res <- mst_loglik(params, X = X, y = y)
    if(is.na(res)) {
      return(high_value)
    } else {
      return(res)
    }
  }
  
  # fit model and record runtime
  mbm = microbenchmark::microbenchmark(
    test = {fit <- optim(init, mst_loglik_wrapper, method = method,
                         control = list(maxit = max_iter, trace = trace))},
    times = 1
  )
  
  # tidy parameters
  coef_mle <- fit$par
  last <- length(coef_mle)
  
  dp <- coef_to_dp(coef_mle[1:(last - 2)], d = ncol(X))
  b0 <- coef_mle[last - 1]
  tau <- exp(coef_mle[last])
  
  # compute mode and ld_max 
  mode <- modeSECdistr(dp, "ST")
  ld_max <- dmst(mode, dp = dp, log = TRUE)
  
  # tidy result
  out <- list(
    optim_result = fit,
    coef_mle = fit$par,
    params_mle = list(dp = dp, b0 = b0, tau = tau, 
                      ld_max = ld_max, mode = mode),
    y = y,
    X = X,
    runtime = mbm
  )
  
  if(fit$convergence > 0) print("WARNING: Optimization did not converge.")
  
  return(out)
}


# simulate function. It may work with the parameters (dp, b0 and tau)
# or with a fitted model. It returns a matrix, where columns are replications.
# If a model is provided, the parameters are overridden.
mst_simulate <- function(dp, b0, tau, model = NULL, 
                         newdata = X, nsim = 100) {
  
  if(!is.null(model)) {
    b0 <- model$params_mle$b0
    dp <- model$params_mle$dp
    tau <- model$params_mle$tau
  }
  
  mu <- mst_mu(dp, b0, newdata)
  y_sim <- matrix(
    rnorm(nrow(newdata) * nsim, mu, tau), 
    nrow(newdata), nsim
  )
  
  return(y_sim)
}

# dharma function. it returns the dharma residuals and a DHARMa-like plot.
mst_dharma <- function(model) {
  
  mu <- mst_mu(b0 = model$params_mle$b0, 
               dp = model$params_mle$dp, 
               X = model$X)
  res <- pnorm(model$y, mean = mu, sd = model$params_mle$tau)
  
  res_ord <- res[order(res)]
  mu_ecdf <- ecdf(mu)
  mu_rank <- mu_ecdf(mu)
  
  par(mfrow = c(1, 2))
  plot(res_ord ~ ppoints(length(res)), 
       xlab = "Expected quantiles", ylab = "Observed quantiles",
       xlim = c(0, 1), ylim = c(0, 1),
       col = rgb(0, 0, 0, 0.1), pch = 19)
  abline(0, 1)
  plot(res ~ mu_rank, 
       xlab = "Fitted values (rank-transformed)",
       ylab = "DHARMa residuals", 
       xlim = c(0, 1), ylim = c(0, 1),
       col = rgb(0, 0, 0, 0.1), pch = 19)
  par(mfrow = c(1, 1))
  
  return(res)
}

# function to make parameter bounds (in the constrained scale) given the problem
# dimension.
bounds_unique <- data.frame(
  #       xi,  corr, sigma, alpha, nu, b0, tau
  lower = c(-50, -1, 0, -20, 1.5, -100, 0),
  upper = c(50, 1, 100, 20, 100, 100, 100)
)
 
make_bounds <- function(d, bounds = bounds_unique) {
  reps <- c(d, cor_num(d), d, d, 1, 1, 1)
  bd <- data.frame(
    lower = rep(bounds$lower, times = reps),
    upper = rep(bounds$upper, times = reps),
    param = rep(c("xi", "cor", "sigma", "alpha", "nu", "b0", "tau"),
                times = reps)
  )
  return(bd)
}

# Transform parameters from the constrained scale to the unconstrained one.
# these params include b0 and tau
unconstrain <- function(params, d) {
  
  # param order: xi, rho, sigma, alpha, nu, b0, tau
  n_cor <- cor_num(d)
  
  # get positions (ordered as xi, rho, sigma, alpha, nu)
  last <- length(params)
  
  xi_ids <- 1:d
  rho_ids <- (xi_ids[d] + 1) : (xi_ids[d] + n_cor)
  sigma_ids <- (rho_ids[length(rho_ids)] + 1) : (rho_ids[length(rho_ids)] + d)
  alpha_ids <- (sigma_ids[d] + 1) : (sigma_ids[d] + d)
  nu_id <- alpha_ids[d] + 1
  b0_id <- last - 1
  tau_id <- last

  # unconstrained parameters
  coeff <- c(
    params[xi_ids],
    params[rho_ids] %>% logit_cor,
    params[sigma_ids] %>% log,
    params[alpha_ids], 
    params[nu_id] %>% log,
    params[b0_id],
    params[tau_id] %>% log
  )
  
  return(coeff)
}

# wrapper for ga(). In this case, params are in the constrained scale, so they are
# unconstrained before passing them to mst_loglik. In addition, the log_lik is 
# returned, not the -log_lik
mst_loglik_ga <- function(params_const, d) {
  coeff <- unconstrain(params_const, d = d)
  # turn into unconstrained scale
  return(mst_loglik(coeff, X = X, y = y) * (-1))
}

# Simulate data and fit model ---------------------------------------------

d <- 8
N <- 2000
X <- matrix(rnorm(N * d), N, d)

xi <- rnorm(d)
rho <- runif(cor_num(d), -1, 1)
sigma <- rexp(d)
alpha <- rnorm(d)
nu <- runif(1, 1.5, 20)

dp <- params_to_dp(xi, rho, sigma, alpha, nu)
b0 <- rnorm(1)
tau <- rexp(1)

coeff <- dp_to_coef(dp)
coef_extended <- c(coeff, b0, log(tau))

# Compute mean
mu <- mst_mu(dp, b0, X)

# simulate data
y <- rnorm(N, mu, tau)


mst_loglik(coef_extended, X, y)


# Fit model
start <- coef_extended + rnorm(coef_num(d) + 2, sd = 0.1)
fit1 <- mst_fit(start, X, y)
mst_dharma(fit1)
plot(fit1$coef_mle ~ coef_extended); abline(0, 1)
fit1$runtime
fit1$optim_result$convergence

# Fit using gradient method (BFGS)
fit2 <- mst_fit_grad(coef_extended + rnorm(coef_num(d) + 2, sd = 0.1), 
                     X, y)
mst_dharma(fit1)
plot(fit1$coef_mle ~ coef_extended); abline(0, 1)
fit1$runtime
fit1$optim_result$convergence


# Fit model with GA. it requires using bounds on the parameter space.
# Use coeff in the constrained scale.
params_extended <- c(xi, rho, sigma, alpha, nu, b0, tau)
mst_loglik_ga(params_extended, d) # 884.5099

bb <- data.frame(
  #       xi,  corr, sigma, alpha, nu, b0, tau
  lower = c(-10, -1, 0, -10, 1.5, -20, 0),
  upper = c(10, 1, 20, 10, 20, 20, 20)
)

ga1 <- ga(type = "real-valued", fitness = mst_loglik_ga, d = d, 
          lower = make_bounds(d, bb)$lower, upper = make_bounds(d, bb)$upper,
          popSize = 5000, maxiter = 10, parallel = 14,
          pmutation = 0.5,
          elitism = round(5000 * 0.1),
          suggestions = matrix(params_extended, nrow = 1))

plot(ga1)

plot(as.numeric(ga1@solution) ~ params_extended); abline(0, 1)
plot(unconstrain(as.numeric(ga1@solution), d) ~ unconstrain(params_extended, d)); abline(0, 1)
ga1@summary

# GAM fitting a MST --------------------------------------------------------

d <- 2
xnames <- paste("x", 1:d, sep = "")
N <- 2000
X <- matrix(rnorm(N * d, sd = 5), N, d)

xi <- rnorm(d)
rho <- -0.9#runif(cor_num(d), -1, 1)
sigma <- rexp(d)
alpha <- rep(0, d)#rnorm(d)
nu <- 1000#runif(1, 1.5, 20)

dp <- params_to_dp(xi, rho, sigma, alpha, nu)
distr <- makeSECdistr(dp, family = "ST")
# plot(distr)

b0 <- rnorm(1)
tau <- 0.1#rexp(1)

coeff <- dp_to_coef(dp)
# coef_to_dp(coeff, d = ncol(X))

coef_extended <- c(coeff, b0, log(tau))

# Compute mean
mu <- mst_mu(dp, b0, X)

# simulate data
y <- rnorm(N, mu, tau)


## Tidy data for GAM

dat <- as.data.frame(cbind(y, X))
names(dat) <- c("y", xnames)

## Fit mst model
# fit1 <- mst_fit(coef_extended + rnorm(coef_num(d) + 2, sd = 0.00001), 
#                 X, y)
### # TIRA UN ERROR QUE NO ENTIENDOOOOOO
# mst_dharma(fit1)
# plot(fit1$coef_mle ~ coef_extended);abline(0, 1)
# fit1$runtime


## Fit GAM
k_marg <- 10
k_int <- 4
# time_start <- Sys.time()
# gam1 <- gam(
#   y ~ s(x1, bs = "cr", k = k_marg) + 
#       s(x2, bs = "cr", k = k_marg) + 
#       ti(x1, x2, bs = "cr", k = k_int),
#   method = "REML", data = dat
# )
# time_end <- Sys.time()
# (time_gam <- time_end - time_start)
# plot(gam1)

# # try a shape-constrained gam (scam)
gam1 <- scam(
  y ~ s(x1, bs = "cv", k = k_marg) +
      s(x2, bs = "cv", k = k_marg) +
      ti(x1, x2, bs = "ts", k = k_int, m = 1),
  data = dat
)

# Compare GAM with real curve in the range

nd <- expand.grid(x1 = seq(min(dat$x1) * 3, max(dat$x1) * 3, length.out = 1000),
                  x2 = 0)

gam_pred <- predict(gam1, nd, type = "response")
# gam2_pred <- predict(gam2, nd, type = "response")
mst_pred <- mst_mu(dp, b0, as.matrix(nd))

# plot(gam_pred ~ mst_pred)

plot(mst_pred ~ nd$x1, type = "l")
lines(gam_pred ~ nd$x1, col = 2)
# lines(gam2_pred ~ nd$x1, col = 4)
abline(v = c(min(dat$x1), max(dat$x1)))



# evaluate the 2d shape (is gam good for corr structure?)
nd2 <- expand.grid(
  x1 = seq(min(dat$x1) * 3, max(dat$x1) * 3, length.out = 100),
  x2 = seq(min(dat$x2) * 3, max(dat$x2) * 3, length.out = 100)
)

gam_pred_bi <- predict(gam1, nd2, type = "response")
# gam2_pred <- predict(gam2, nd, type = "response")
mst_pred_bi <- mst_mu(dp, b0, as.matrix(nd2))

# data 
data_bi <- cbind(rbind(nd2, nd2), y = c(gam_pred_bi, mst_pred_bi)) %>% as.data.frame
data_bi$model <- rep(c("gam", "mst"), each = nrow(nd2))

ggplot(data_bi, aes(x = x1, y = x2, z = y, fill = y)) + 
  geom_tile() +
  scale_fill_viridis() +
  facet_wrap(vars(model))



# jugando con rmvn --------------------------------------------------------

# (para intentar con una normal multivariada)
d <- 3
ff <- function() is_positive_definite(make_vcov_raw(runif(cor_num(d), -0.99, 0.99), rexp(d, 0.1) + 0.001))
many_mats <- replicate(1e5, ff())
sum(unlist(many_mats))

# Incluso forzando valores razonables de rho y sigma, la vcov sigue no siendo
# siendo positive definite.

d <- 3
ffs <- function() is_positive_semi_definite(make_vcov_raw(runif(cor_num(d), -0.9, 0.9), rexp(d, 0.1) + 1))
many_mats <- replicate(1e5, ffs())
sum(unlist(many_mats))
# Y si logramos una positive-semi-definite para usar con mgcv? (dmvn)

rmvn(1, c(0, 0), m)

getAnywhere(dmvn)

# Optim gradient con Inf ----------------------------------------------------

try_inf <- function(x) {
  pp <- rbinom(1, prob = 0.9, size = 1)
  if(pp == 1) {
    return(x ^ 2)
  } else {
    return(1000)
  }
}


optim(0.3, try_inf, method = "BFGS")


# Parece que la frecuencia de valores muy altos (malos) y su magnitud
# afectan el ajuste.


# trials, plots ------------------------------------------------------------


# a <- coef_to_dp(rnorm(coef_num(9)), d = 9)
# a$Omega %>% is.positive.definite()
# a$Omega %>% is_positive_definite()
# 
# dpp <- params_to_dp(xi = rep(1, 3), 
#              rho = rep(0.5, cor_num(3)),
#              sigma = rep(3, 3),
#              alpha = rep(-2, 3),
#              nu = 19)
# dpp$Omega
# 
# dp_to_coef(dpp)
# dp_to_coef(dpp) %>% length == coef_num(3) # OK
# 
# 
# # Try with a mst distrib
# 
# dd <- 4
# dp_sim <- params_to_dp(
#   xi = rep(0, dd),#rnorm(dd, sd = 5),
#   rho = runif(cor_num(dd), -1, 1),
#   sigma = rexp(dd),
#   alpha = rnorm(dd, sd = 10),
#   nu = rexp(1, 0.1) + 0.9 # make reasonable tails
# )
# 
# d_trial <- makeSECdistr(dp = dp_sim, 
#                         family = "ST", name = "Multiv.ST",
#                         compNames = c("X1", "X2", "X3", "X4"))
# plot(d_trial)
# dp_sim$nu
# 
# # cauchy worse 1D
# d_trial <- makeSECdistr(dp = c(xi = 0, omega = 2, alpha = -1, nu = 0.7),
#                         family = "ST", name = "Cauchy worse")
# plot(d_trial)

# TAREA -------------------------------------------------------------------

# hacer una función que estime el modelo y devuelva muchas cosas con el modelo
# ajustado, como predict(), simulate(), residuals(),   
# R2, inverse_link(), plot(), emulated_loglik_max() -es decir, el modo de la
# log-density function-, posterior_flat_samples(), posterior_flat_summary().

# (hice bastante pero faltan un par)

# Problema: ajustar un mst con optim no es fácil porque hace una búsqueda muy 
# local, y se ve que la loglik tiene muchos optimos locales. El mle depende
# mucho del starting point. 

# El GAM ajusta razonable en el rango de los datos, a veces oversmoothing... la 
# linealidad que asume al extrapolar es bastante distinta a la curva de la mst.
# Me parece mejor usar una optimización global.

# Probar optimización global con GA:
# https://cran.r-project.org/web/packages/GA/vignettes/GA.html#parallel-computing

# quizás podría restringir los parámetros a valores razonables? 
# (serviría usar la previa para los xi!).
# evaluar la performance del GA contra la de optim(). Debería ser mejor y menos
# variable. NO, ANDA RE MAL

# GA puede correr en paralelo :) y sirve con funciones no diferenciables.
# Resultó muy malo para encontrar el óptimo, está re re lejos.


# Evaluar si una GAM con ti puede ajustar (en forma) una normal multivariada.
# NO, no puede. Incluso con shape constraint, no podría modelar bien la 
# correlación entre parámetros. Cuando te vas de rango, las bases de interacción
# (no restringidas) se van a la miercole, incluso con el signo incorrecto.


# SÍ O SÍ HAY QUE USAR UNA FUNCIÓN DECRECIENTE HACIA LOS EXTREMOS 
# (TIPO NORMAL MULTIVARIADA)

# Simplificar: probar con una normal multivariada

# El problema de la positive definiteness viene de la PDF de la normal 
# multivariada: La constante normalizadora requiere det(Sigma), y la parte
# rica de la función necesita Sigma ^ (-1), es decir, la matriz de precisión.
# Podría escribir la pdf a mano directamente usando matrices de precisión. 
# Mientras la diagonal y las corr condicionales sean sensatas, no debería 
# tener problema.

# Otra idea: Seguir con GA, pero híbrido, con "L-BFGS-B". Esto haría una búsqueda
# más eficiente. El truco para usar gradientes es que la función devuelva 
# valores |grandes| en vez de +-Inf. La escala de ese valor y la frecuencia 
# afectan la calidad de la estimación. Entonces, para definirlo, podemos correr
# una generación de GA con varias partículas (10000) y ver cuál es el peor valor. 
# Y hacemos que devuelva ese valor en vez de Inf. Si esto no anda, o es inestable,
# intentar ajustar una mvn 
# - IMPORTANTE - parece haber bugs, porque cuando intento probar una optim con 
# BFGS me dice esto:
# ---Error in eigen(m, symmetric = TRUE, only.values = TRUE) : 
# ---infinite or missing values in 'x'
# También me está tirando errores raros el Nelder-Mead.

d <- 8
cor_num(d) + d * 2 + 2 # = 46
# params para una normal multivariada

d <- 9
cor_num(d) + d * 2 + 2 # = 56

# hardcoding la density podría optimizarla en Stan... 
# Y así podría agregarle unas splines marginales para darle más flexibilidad,
# que estén centradas en la func determinística. O sea, que sean splines que 
# tienden a cero. 
# Una b spline de orden 0 andaría bien, tipo
# s(x, k = 50, bs = "bs", m = c(3, 0))
# ver   https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/


# ENTONCES, seguir así:

# Ajustar a mano la density MVN con precisión.
# Ver si surgen problemas.
# 
# Luego hacerlo en Stan y sumarle splines marginales centradas en cero.
# (radial basis functions)

# El problema de las marginales es que no servirían mucho para las correlaciones.
# Entonces podría intentar hacer tensor products a mano. Pero ya usaría un montón
# de parámetros...

# Probar primero una MVN, luego agregarle una spline, y luego verás



# building bases centred at 0  --------------------------------------------

x <- rnorm(100)
K <- 10
bs <- smoothCon(s(x, k = K, bs = "cr", m = 1), data = data.frame(x = x),
                diagonal.penalty = F, absorb.cons = TRUE)[[1]]
nd <- data.frame(x = seq(-10, 10, length.out = 300))
bs_mat <- PredictMat(bs, data = nd)
matplot(bs_mat, type = "l")
b <- rnorm(K-1)
cur <- bs_mat %*% b
plot(cur ~ nd$x, type = "l")


x <- rnorm(100)
K <- 10
bs <- smoothCon(s(x, k = K, bs = "ts", m = c(3, 0)), data = data.frame(x = x),
                diagonal.penalty = T, absorb.cons = TRUE)[[1]]
nd <- data.frame(x = seq(-10, 10, length.out = 300))
bs_mat <- PredictMat(bs, data = nd)
matplot(nd$x, bs_mat, type = "l")
b <- rnorm(K)
cur <- bs_mat %*% b
plot(cur ~ nd$x, type = "l")


x <- runif(200, -2, 2)
knots <- list(x = c(-5, -2, 2, 5))
K <- 10
bs <- smoothCon(s(x, k = K, bs = "bs", m = c(3, 0)), data = data.frame(x = x),
                diagonal.penalty = T, absorb.cons = TRUE,
                knots = knots)[[1]]
nd <- data.frame(x = seq(-10, 10, length.out = 300))
bs_mat <- PredictMat(bs, data = nd)
matplot(nd$x, bs_mat, type = "l")
abline(v = c(-5, 5), col = 4)
abline(v = c(-2, 2), col = 2)
b <- rnorm(K-1)
cur <- bs_mat %*% b
plot(cur ~ nd$x, type = "l")
abline(v = c(-5, 5), col = 4)
abline(v = c(-2, 2), col = 2)

# Son una caca. usar radial basis functions, que se anulan fuera de dominio.
