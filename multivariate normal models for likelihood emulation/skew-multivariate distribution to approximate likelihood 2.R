# The emulated likelihood must be decreasing in the range out of data, and this
# is not achievable by any kind of GAM (including shape-constrained ones) if
# the correlation between parameters is meant to be considered. Hence, a 
# multivariate normal shape is needed. The problem with that shape is that 
# estimating such a function is not trivial using the density functions 
# from the mgcv or sn packages because they require a positive-(semi)-definite
# covariance matrix. If the objective function to be optimized finds an invalid
# Sigma, it must return NA, Inf or some bad value. This makes the gradiend-based
# optimization techniques invalid, and non-gradient algorithms take too much and 
# are very local. 

# A way to overcome the need for a positive-definite matrix is by parameterizing
# directly the precision matrix (P = solve(Sigma)), which is used in the
# computation of the multivariate normal density. In this way, if the density
# is written by hand, there should be no problems. 

# Another option would be to trust Stan, and use it as optimizer. Here, a
# correlation or covariance matrix could be defined, and Stan should take care
# of optimizing over their valid parameters space. However, it requires to write
# the skew-multivariate-normal density by hand (although in Sigma scale, not 
# precision.)
# GO FOT IT, Stan is great for optimizing over positive-definite matrices.
 
# Consider defining the skew-distribution with centred parameterization.
# (doesn't seem easy)

# hasta ahora, lo más difícil de la skew-normal, parece ser calcular el modo. 
# Se puede hacer by hand, it's explained in the book, but it doesn't seem so 
# trivial. Another option is not using Stan, so we can use {sn} to compute the
# mode, but that requires to optimize over pd-matrices by hand.
# Another option is not to set the maximum of the raw function at zero. In this 
# way, the intercept will be correlated with the slant, but that would probably
# not be a problem.
# We can copy the R code to Stan :) That helps

# NOTE: quad_form_diag(Omega, tau) == diag(tau) %*% Omega %*% diag(tau),
# but the first is more efficient (tau are the sigmas, and Omega is the corr matrix)


# OJO, empezar simple. primero normal multivariada, en R y en Stan. Ver si eso
# funciona muy distinto en R o en Stan, evaluar si es estimable. Si funciona,
# decidir si vale la pena Stanear.
# Si todo eso anda muy bien, considerar una skew-normal.


# Packages ----------------------------------------------------------------

library(sn)
library(Matrix)     # to find nearest positive-definite matrix
library(matrixcalc) # test positive-definiteness (is.positive.definite)
library(tidyverse); theme_set(theme_bw())
library(viridis)
library(GA)         # global optimization


# functions ---------------------------------------------------------------

# Link functions used by Stan for correlation parameters:

# tanh(c(-Inf, Inf))
# atanh(c(-1, 1))
# curve(tanh(x), from = -5, to = 5)
# curve(plogis(x) * 2 - 1, add = TRUE, col = 2)

# function to compute the number of correlation parameters as a function of the
# multivariate distribution dimension
cor_num <- function(d) ncol(combn(d, 2))

# inverse function: detect dimension of the matrix from the length of the 1D 
# vector of correlation coefficients
cor_dim <- function(l) 1 / 2 * (sqrt(8 * l + 1) + 1)
# using wolframalpha: https://www.wolframalpha.com/widgets/view.jsp?id=c86d8aea1b6e9c6a9503a2cecea55b13

# function to compute the total number parameters to define a multivariate 
# skew-t distribution as a function of is dimension
# (It doen't include b0 and tau, used to fit a mst curve.)
coef_num <- function(d) d * 3 + 1 + cor_num(d)

# positive-definiteness check made by makeSECdistr {sn}. It differs from 
# is.positive.definite(), from {matrixcalc}
is_positive_definite <- function(m) {
  min(eigen(m, symmetric = TRUE, only.values = TRUE)$values) > 0
}

# make correlation matrix from unconstrainde partial correlation parameters, 
# matrix dimension and bound on partial correlation coefficients. Although 
# tanh transforms from (-Inf, Inf) to (-1, 1), when such values appear the 
# resulting matrix might be non-positive-definite.
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

# the same but already bounding the raw cor
make_corr_bounded <- function(y, d) {
  # # TEST
  # d <- 4
  # y <- rnorm(cor_num(d), sd = 1e5)
  # cor_bound <- 0.95
  # # END TEST
  
  # upper triangular matrix of Canonical Partial Correlations
  z <- matrix(0, d, d)
  z[upper.tri(z)] <- tanh(y) 
  
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


# Function to make vcov based on unconstrained partial correlation coefficients
# and log-sd vector. the sigma-bound is used to prevent sds near 0.
make_vcov_corr <- function(cor_vector, log_sigma, sigma_bound = 0.02) {
  
  # # TEST
  # d <- 4
  # cor_vector <- rnorm(cor_num(d), sd = 1e5)
  # log_sigma <- rnorm(d, sd = 1e-3)
  # sigma_bound <- 0.02
  # # END TEST
  
  d <- length(log_sigma)
  Rho <- make_corr(cor_vector, d)
  sigma <- exp(log_sigma) + sigma_bound
  Sigma <- (sigma %*% t(sigma)) * Rho
  
  # # Check
  # is_positive_definite(Sigma); isSymmetric(Sigma)
  # # End check
  
  return(Sigma)
}

# The same but with bounded params
make_vcov_corr_bounded <- function(cor_vector, sigma) {
  
  # # TEST
  # d <- 4
  # cor_vector <- rnorm(cor_num(d), sd = 1e5)
  # log_sigma <- rnorm(d, sd = 1e-3)
  # sigma_bound <- 0.02
  # # END TEST
  
  d <- length(sigma)
  Rho <- make_corr_bounded(cor_vector, d)
  Sigma <- (sigma %*% t(sigma)) * Rho
  
  # # Check
  # is_positive_definite(Sigma); isSymmetric(Sigma)
  # # End check
  
  return(Sigma)
}


# Test make_vcov_corr
vcov_test <- function() {
  d <- 4
  cor_vector <- rnorm(cor_num(d), sd = 1e5)
  log_sigma <- rnorm(d, sd = 1e4)
  make_vcov_corr(cor_vector, log_sigma) %>% is_positive_definite()
}

# sum(replicate(1e4, vcov_test()) %>% unlist) / 1e4 * 100 
# OK, works in bad conditions.

# coef_to_dp converts the raw model coefficients (all at the unconstrained scale) 
# to the dp parameterization of the multivariate skew-t distribution.
# dp stands for "direct parameterization", used in the sn package.
# coef is the unconstrained parameters vector that will be passed to optim(),
# excluding the intercept ("b0") and the residual standard deviation 
# ("tau").
# d is the dimension of the multivariate density function to fit.
# a lower bound on the degrees of freedom ("nu") is set to avoid rare shapes.
coef_to_dp <- function(coef, d, nu_bound = 1.5) {
  
  ## TEST
  # d <- 5
  # ll <- cor_num(d) + d * 3 + 1
  # coef <- rnorm(ll, sd = 2)
  ## END TEST
  names(coef) <- rep(c("xi", "rho", "sigma", "alpha", "nu"), 
                     c(d, cor_num(d), d, d, 1))
  
  # make Omega vcov
  Omega <- make_vcov_corr(coef[names(coef) == "rho"], 
                          coef[names(coef) == "sigma"])
  
  # if(!is_positive_definite(Omega)) warning("Omega is not positive definite.")
  
  # alpha is not constrained
  
  # df
  nu <- exp(coef[names(coef) == "nu"]) + 2  # log link for nu
  
  dp <- list(
    xi = coef[names(coef) == "xi"],
    Omega = Omega, 
    alpha = coef[names(coef) == "alpha"],
    nu = nu
  )
  
  return(dp)
}
# Check:
# d <- 9
# coef_to_dp(rnorm(coef_num(d), sd = 100), d)

# The same, but from bounded coeficients
coef_bounded_to_dp <- function(coef, d) {
  
  # make Omega vcov
  Omega <- make_vcov_corr_bounded(coef[names(coef) == "rho"], 
                                  coef[names(coef) == "sigma"])
  
  # if(!is_positive_definite(Omega)) warning("Omega is not positive definite.")

  dp <- list(
    xi = coef[names(coef) == "xi"],
    Omega = Omega, 
    alpha = coef[names(coef) == "alpha"],
    nu = coef[names(coef) == "nu"]
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
# are the intercept (b0) and the (log) residual standard deviation (tau), 
# corresponding to a Normal data model.
# X is a matrix with the fire-spread-parameters values, and y is the observed
# similarity between simulated and observed fires.
# It returns the log-likelihood, not negated, to be used as 
# optim(..., control = list(fnscale = -1)), and for compatibility with ga(), which
# maximizes.
mst_loglik <- function(coef_extended, X, y) {
  
  # separate parameters
  last <- length(coef_extended)
  dp <- coef_to_dp(coef_extended[1:(last - 2)], d = ncol(X)) # remove b0 and tau
  b0 <- coef_extended[last - 1]
  tau <- exp(coef_extended[last])
  
  # compute mu
  mu <- mst_mu(dp = dp, b0 = b0, X = X)
  
  # Compute -log-likelihood for the simulated likelihood. 
  ll <- sum(dnorm(y, mu, tau, log = TRUE))
  
  return(ll)
}

# Function to fit mst model.
# X is the predictors matrix, y is the response variable,
# init is the initial value for the parameters, and 
# max_iter is the maximum number of iterations for optim(), and trace
# the detail level for printing the progress (see ?optim).
# method is the method to optimize. defaults to "BFGS"
# value: list with many objects...
mst_fit <- function(init, X, y, max_iter = 1e5, trace = 1, method = "L-BFGS-B",
                    lower = NULL, upper = NULL) {
  
  # wrapper for optim()
  mst_loglik_wrapper <- function(params) {
    mst_loglik(params, X = X, y = y)
  }
  
  # fit model and record runtime
  mbm = microbenchmark::microbenchmark(
    test = {fit <- optim(init, mst_loglik_wrapper, method = method,
                         lower = lower, upper = upper,
                         control = list(maxit = max_iter, trace = trace,
                                        fnscale = -1))},
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
  
  # Compute R2
  mu <- mst_mu(dp, b0, X)
  R2 <- var(mu) / (var(mu) + var(y - mu))
  
  # tidy result
  out <- list(
    optim_result = fit,
    coef_mle = fit$par,
    params_mle = list(dp = dp, b0 = b0, tau = tau, 
                      ld_max = ld_max, mode = mode),
    y = y,
    X = X,
    fitted = mu,
    R2 = R2,
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

# Function to simulate reasonable values to simulate data and to simulate
# starting points for estimation.
simulate_coef <- function(d) {
  # unconstrained parameters
  xi <- rnorm(d)
  rho_raw <- runif(cor_num(d), -1, 1) # tanh reaches limits near [-2, 2]
  log_sigma <- log(runif(d, 0.5, 5))
  alpha <- rnorm(d, sd = 1)
  log_nu <- log(runif(1, 2, 100))
  b0 <- rnorm(1, 5)
  log_tau <- log(runif(1, 0.3, 4))
  
  # merge
  coeff_extended <- c(xi, rho_raw, log_sigma, alpha, log_nu, b0, log_tau)
  names(coeff_extended) <- rep(
    c("xi", "rho_raw", "log_sigma", "alpha", "log_nu", "b0", "log_tau"), 
    c(d, cor_num(d), d, d, 1, 1, 1)
  )
  
  return(coeff_extended)
}


# the same, but bounding parameters in the original scale.
simulate_coef_bounded <- function(d) {
  # unconstrained parameters
  xi <- runif(d, -20, 20)
  rho_partial <- runif(cor_num(d), -1.5, 1.5) # unconstrained corr coeff
  sigma <- runif(d, 0.5, 5)
  alpha <- runif(d, -8, 8)
  nu <- log(runif(1, 2, 100))
  b0 <- runif(1, -10, 10)
  tau <- runif(1, 0.3, 1)
  
  # merge
  coeff_extended <- c(xi, rho_partial, sigma, alpha, nu, b0, tau)
  names(coeff_extended) <- rep(
    c("xi", "rho", "sigma", "alpha", "nu", "b0", "tau"), 
    c(d, cor_num(d), d, d, 1, 1, 1)
  )
  
  return(coeff_extended)
}

# Simulate data and fit model ---------------------------------------------

L <- rep(c(-20, -1.5, 0.5, -8, 2, -10, 0.3))
U <- rep(c(20, 1.5, 5, 8, 100, 10, 1))


d <- 8
N <- 2000
X <- matrix(rnorm(N * d), N, d)

# unconstrained parameters
coeff_extended <- simulate_coef_bounded(d)
coeff <- coeff_extended[1:(length(coeff_extended) - 2)]
dp <- coef_bounded_to_dp(coeff, d = d)
fam <- makeSECdistr(dp, family = "ST")
# plot(fam)
modeSECdistr(dp, "ST")
# cov2cor(dp$Omega)[lower.tri(cov2cor(dp$Omega))] %>% range

b0 <- coeff_extended[names(coeff_extended) == "b0"]
tau <- coeff_extended[names(coeff_extended) == "tau"]

# Compute mean
mu <- mst_mu(dp, b0, X)

# simulate data
y <- rnorm(N, mu, tau)

(R2_obs <- var(mu) / (var(mu) + var(y - mu)))
mst_loglik(coeff_extended, X, y)

# Resample to fit snlm
size <- N * 10
ids_boot <- sample(1:length(y), size = size, 
                   prob = dnorm(y, mean = mu, sd = tau),
                   replace = TRUE)

X_boot <- X[ids_boot, ]
m_aprox <- selm(cbind(V1, V2, V3, V4, V5, V6, V7, V8) ~ 1, family = "ST", 
                data = as.data.frame(X_boot))
# Tarda bastante bastante.
summary(m_aprox, "dp")

coeff_aprox <- coef(m_aprox, vector = FALSE, param.type = "DP")
coeff_aprox$beta ; dp$xi

coeff_aprox$alpha ; dp$alpha

coeff_aprox$Omega ; dp$Omega

# Fit model

# move real parameters by some reasonable sd
start <- simulate_coef_bounded(d)
fit1 <- mst_fit(start, X, y, lower = L, upper = U)
# fit2 <- mst_fit(coeff_extended, lower = L, upper = U,# + rnorm(length(coef_extended), sd = 0.05), 
#                 X, y)

# NOW IT WORKS, BOUNDIN THE PARAMETERS, BUT IT GETS STUCK IN BAD LOCAL MINIMA.
# USE GOOD STARTING POINTS AND TRY MANY (MAYBE A SOBOL SEQUENCE OVER A 
# DECENT RANGE).

# 6400 POINTS, AND CHOOSE THE BEST 10 TO OPTIMIZE, PLUS A DECENT POINT.
# MAYBE RESAMPLING THE LIKELIHOOD DATA ONE CAN GET SOME INSIGHT ON THE 
# MARGINAL MEANS AND SDS AND USE THEM AS STARTING POINTS.

# BEFORE TRY THIS: RESAMPLE DATA ACCORDING TO LIKELIHOOD, COMPUTE 
# CORRELATIONS, SDS, MEANS TO HAVE GOOD STARTING POINTS. 
# INCLUSO, AJUSTAR UN SNLM CON SN COMO UNA PRIMERA APROXIMACIÓN.
# LO MALO DE ESO ES QUE LUEGO LAS DISTRI ESTIMADAS PUEDEN ESTAR TRUNCADAS, 
# Y PARA ESO HABRIA QUE ESTIMARLO EN ESCALA LOG EN SN.


# It has to be run many times because many starting points are bad and fail.
# Maybe we need to think better starting points.
# Consider using try.catch.

# It finds local minima and gets stuck at it. Try many starting points and
# better parameter values.

# The L-BFGS-B algorithm is too slow, maybe it's worthwhile to run it in Stan, 
# because iterations run much much slower than any hard model in Stan.


mst_dharma(fit1)
plot(fit1$coef_mle ~ coeff_extended); abline(0, 1)
fit1$runtime
fit1$optim_result$convergence
mu_fit <- fit1$fitted
plot(mu_fit ~ mu, pch = 19, col = rgb(0, 0, 0, 0.1),
     xlab = "Truth", ylab = "Estimation"); abline(0, 1)
# tiró warnings sobre Omega is not positive definite, pero convirgió sin drama.


