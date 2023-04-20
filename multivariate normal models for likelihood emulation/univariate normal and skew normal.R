# Fit univariate symmetric and skew normal function with Stan.
library(tidyverse); 
library(ggplot2); library(magrittr); theme_set(theme_bw())
library(viridis)
library(sn)
library(rstan)
library(extraDistr)  # half t density (dht)
library(DHARMa)
library(brms)

focal_folder <- "multivariate normal models for likelihood emulation"

# functions ---------------------------------------------------------------

# unnormalized log-density of a normal distribution
normal_lupdf <- function(x, xi, sigma) {
  0.5 * (-(1 / sigma ^ 2) * (x - xi) ^ 2)
}

# unnormalized log-density of a skew-normal distribution
skew_normal_lupdf <- function(x, xi, sigma, alpha) {
  normal_log_density <- 0.5 * (-(1 / sigma ^ 2) * (x - xi) ^ 2)
  normal_log_cdf_shifted <- log(pnorm(alpha * ((x - xi) / sigma)))
  skew_log_den <- log(2) + normal_log_density + normal_log_cdf_shifted
  return(skew_log_den)
}

# unnormalized log-density of a skew-normal distribution
skew_normal_updf <- function(x, xi, sigma, alpha) {
  normal_density <- exp(0.5 * -(1 / sigma ^ 2) * (x - xi) ^ 2)
  normal_cdf_shifted <- pnorm(alpha * ((x - xi) / sigma))
  skew_den <- 2 * normal_density * normal_cdf_shifted
  return(skew_den)
}

# Simulate and fit normal curve -------------------------------------------

# compile stan model
smodel1 <- stan_model(file.path(focal_folder, "univariate normal.stan"), 
                      verbose = TRUE)

# xi is going to be the location parameter, which is usually mu, but we don't
# use mu not to confuse with the mean of the model.
xi <- rnorm(1)
sigma_prior_sd <- 20
sigma_prior_nu <- 3
sigma <- 0.1#rht(1, nu = sigma_prior_nu, sigma = sigma_prior_sd)
# prior for sigma
curve(dht(x, nu = 3, sigma = sigma_prior_sd), from = 0, to = 100)

N <- 2000
x <- runif(N, -5, 5)

mu <- normal_lupdf(x, xi, sigma)  

# simulate tau, the residual standard deviation based on sd(mu) to get
# an R2 near 0.8
tau <- sd(mu) * 0.5
y <- rnorm(N, mu, tau)
(R2 <- var(mu) / (var(mu) + var(y - mu)))

plot(y ~ x, col = rgb(0, 0, 0, 0.2), pch = 19)
curve(normal_lupdf(x, xi = xi, sigma = sigma), add = TRUE, lwd = 1.8)

sdata1 <- list(
  # data
  x = x, y = y, N = N,
  
  # prior parameters
  xi_prior_sd = 50,
  sigma_prior_sd = sigma_prior_sd,
  sigma_prior_nu = sigma_prior_nu,
  tau_prior_sd = 500,
  tau_prior_nu = 3
)

m1 <- sampling(smodel1, data = sdata1, seed = 2123, 
               cores = 2, chains = 2, iter = 2000)

# sm1 <- summary(m1)[[1]]
# print(sm1)
# pairs(m1, pars = c("xi", "sigma", "tau"))
samples1 <- extract(m1)

par(mfrow = c(2, 2))
plot(samples1$xi ~ samples1$sigma, pch = 19, col = rgb(0, 0, 0.8, 0.05),
     xlim = range(c(samples1$sigma, sigma)),
     ylim = range(c(samples1$xi, xi)))
points(xi ~ sigma, pch = 19, cex = 1.5, col = "red")

plot(samples1$xi ~ samples1$tau, pch = 19, col = rgb(0, 0, 0.8, 0.05),
     xlim = range(c(samples1$tau, tau)),
     ylim = range(c(samples1$xi, xi)))
points(xi ~ tau, pch = 19, cex = 1.5, col = "red")

plot(samples1$sigma ~ samples1$tau, pch = 19, col = rgb(0, 0, 0.8, 0.05),
     xlim = range(c(samples1$tau, tau)),
     ylim = range(c(samples1$sigma, sigma)))
points(sigma ~ tau, pch = 19, cex = 1.5, col = "red")
par(mfrow = c(1, 1))

# dharma analysis
n_sim <- length(samples1$xi) - 2
y_sim <- sapply(1:n_sim, function(i) {
  mu <- normal_lupdf(x, samples1$xi[i], samples1$sigma[i])
  y_ <- rnorm(N, mu, samples1$tau[i])
}) # OK, columns are replicates

res <- createDHARMa(y_sim, y)
plot(res)

plot(rowMeans(y_sim) ~ mu, pch = 19, col = rgb(0, 0, 0.8, 0.05),
     ylab = "Estimated mu", xlab = "Real mu")
abline(0, 1)


# Observations: 
# in many cases xi and sigma are correlated.

# Simulate and fit skew-normal curve ---------------------------------------

# compile stan model
smodel2 <- stan_model(file.path(focal_folder, "univariate skew-normal.stan"), 
                      verbose = TRUE)
smodel3 <- stan_model(file.path(focal_folder, "univariate skew-normal_centred.stan"), 
                      verbose = TRUE)


# problematic parameters:
# xi = -0.4135036
# sigma = 0.611926
# alpha = 3
# x [-5, 5]

# prior for sigma
sigma_prior_sd <- 3
sigma_prior_nu <- 3
curve(dht(x, nu = 3, sigma = sigma_prior_sd), from = 0, to = 100)

# prior for alpha
alpha_prior_sd <- 5

# xi is going to be the location parameter, which is usually mu, but we don't
# use mu not to confuse with the mean of the model.
xi <- rnorm(1)
sigma <- rht(1, nu = sigma_prior_nu, sigma = sigma_prior_sd)
alpha <- -7#rnorm(1, 0, alpha_prior_sd)
N <- 2000
x <- runif(N, -5, 5)

mu <- skew_normal_lupdf(x, xi, sigma, alpha)  

# simulate tau, the residual standard deviation based on sd(mu) to get
# an R2 near 0.8
tau <- sd(mu) * 0.5
y <- rnorm(N, mu, tau)
(R2 <- var(mu) / (var(mu) + var(y - mu)))

plot(y ~ x, col = rgb(0, 0, 0, 0.2), pch = 19)
curve(skew_normal_lupdf(x, xi = xi, sigma = sigma, alpha), add = TRUE, lwd = 1.8,
      col = "blue")
sc <- sd(y)
sdata2 <- list(
  # data
  x = x, y = y / sc, N = N,
  
  # prior parameters
  xi_prior_sd = 10,
  sigma_prior_sd = 3,#sigma_prior_sd,
  sigma_prior_nu = 200,#sigma_prior_nu,
  tau_prior_sd = 20,
  tau_prior_nu = 200,
  alpha_prior_sd = 20
)

m2 <- sampling(smodel2, data = sdata2, seed = 2126, 
               cores = 2, chains = 2, iter = 2000, warmup = 1000)
pairs(m2, pars = c("xi", "sigma", "alpha", "tau"))

# symmetric model to compare
m2sym <- sampling(smodel1, data = sdata2, seed = 2126, 
                  cores = 2, chains = 2, iter = 2000, warmup = 1000)


# very high posterior correlation between sigma and alpha in some replicates, 
# and it is not related to the zero-ness of alpha.
# But this makes no problem. Maybe just a bit harder to sample.
# neffs are lower than before, for sigma and alpha.

# I found a parameter vector which made the estimation quite difficult. Although
# the curve is clearly assymetric, it can't find sesible parameter values.
# Explore curves from complicated parameters
samples2 <- extract(m2, pars = c("xi", "sigma", "alpha"))
samples2sym <- extract(m2sym, pars = c("xi", "sigma"))

plot(y ~ x, col = rgb(0, 0, 0, 0.2), pch = 19, main = "alpha positive")
curve(skew_normal_lupdf(x, xi = xi, sigma = sigma, alpha), add = TRUE, lwd = 1.8,
      col = "blue")
for(i in 1:500) {
  id <- sample(1:length(samples2$alpha), 1)
  curve(sc * skew_normal_lupdf(x, samples2$xi[id], samples2$sigma[id], samples2$alpha[id]), 
        add = TRUE, lwd = 0.5, col = rgb(1, 0, 0, 0.05))
}
# curves without asymmetry to compare:
for(i in 1:500) {
  id <- sample(1:length(samples2sym$sigma), 1)
  curve(sc * skew_normal_lupdf(x, samples2sym$xi[id], samples2sym$sigma[id], 0), 
        add = TRUE, lwd = 0.5, col = rgb(1, 1, 0, 0.05))
}
# Here it's not so clear whether the symmetric model or the underestimated-symmetric
# is better, but DHARMa analysis suggests that the asymmetric is better.


# Compare parameters
head(sm2 <- summary(m2)[[1]])
summ <- as.data.frame(sm2[1:4, c("mean", "2.5%", "97.5%")])
colnames(summ) <- c("mean", "lower", "upper")
summ$observed <- c(xi, sigma, tau, alpha)
ggplot(summ, mapping = aes(x = observed, y = mean, ymin = lower, ymax = upper)) +
  geom_point() + 
  geom_errorbar() +
  geom_abline(intercept = 0, slope = 1)

# dharma analysis
samples2 <- extract(m2, pars = c("xi", "sigma", "tau", "alpha"))
n_sim <- length(samples2$xi)
y_sim <- sapply(1:n_sim, function(i) {
  mu <- skew_normal_lupdf(x, samples2$xi[i], samples2$sigma[i], samples2$alpha[i])
  y_ <- rnorm(N, mu, samples2$tau[i]) * sc
}) # OK, columns are replicates
res <- createDHARMa(y_sim, y)
plot(res)

# asymmetric model
n_sim <- length(samples2sym$xi)
y_sim <- sapply(1:n_sim, function(i) {
  mu <- skew_normal_lupdf(x, samples2sym$xi[i], samples2sym$sigma[i], 0)
  y_ <- rnorm(N, mu, samples2$tau[i]) * sc
}) # OK, columns are replicates
res <- createDHARMa(y_sim, y)
plot(res)




plot(rowMeans(y_sim) ~ mu, pch = 19, col = rgb(0, 0, 0.8, 0.05),
     ylab = "Estimated mu", xlab = "Real mu")
abline(0, 1)



# Notes on the problematic parameters:
# problematic parameter values:
# xi = -0.4135036
# sigma = 0.611926
# alpha = 3
# x [-5, 5]

# It was impossible for the chains to converge due to the scale of the data,
# sd(y): ~78.
# once I replaced y* = y / sd(y), the model converged. However, assymetry is 
# underestimated, and a quite gaussian-like function is estimated. 
# Would it be resolved with more N? It's already high! (2000)
# No, changing the N to 8000 didn't resolve it. Maybe increasing the R2?
# A curious thing is that the prior on alpha is not regularizing, so it
# has no reason to underestimate it.
# No, increasing the R2 to 0.997 did not help to estimate it

# I gave up. For some shapes, the model fails to estimate a truly asymmetric 
# curve.

# It seems that always |alpha| is underestimated and sigma is overestimated.
# with mild asymmetry, it doesn't cause problems, because fitted and real mu
# are the same. The problem is more evident with higher asymmetries.

# perhaps we should use the centred parameterization of the skew-normal, 
# but that is quite much work.


# Extras

# mu computed in R and Stan
mu_stan <- as.matrix(m2sc, "mu") %>% t
all.equal(mu_stan[, 1], 
          skew_normal_log_mu(x, samples2$xi[1], samples2$sigma[1],
                             samples2$alpha[1]))
j = 3
plot(mu_stan[, j] ~ 
     skew_normal_log_mu(x, samples2$xi[j], samples2$sigma[j], samples2$alpha[j]))
abline(0, 1)
# not exactly the same, but almost. It shouldn't be a problem


# try with a centered parameterization on sigma and alpha:
m3 <- sampling(smodel3, data = sdata2, seed = 2123, 
               cores = 2, chains = 2, iter = 2000)
pairs(m3, pars = c("xi", "sigma", "alpha", "tau"))
pairs(m3, pars = c("xi", "sigma_raw", "alpha_raw", "tau"))
# the posterior correlation remains.

# try almost flat priors on alpha and sigma: do they make a problem?
sdata4 <- list(
  # data
  x = x, y = y, N = N,
  
  # prior parameters
  xi_prior_sd = 50,
  sigma_prior_sd = 1e5,
  sigma_prior_nu = sigma_prior_nu,
  tau_prior_sd = 5,
  tau_prior_nu = 3,
  alpha_prior_sd = 1e5
)
m4 <- sampling(smodel2, data = sdata4, seed = 2123, 
               cores = 2, chains = 2, iter = 2000)
pairs(m4, pars = c("xi", "sigma", "alpha", "tau"))
# no, the flat priors do not seem to be a problem for now.





# skew-normal dp and cp ---------------------------------------------------


# perhaps the asymmetry underestimation problem could be improved by using the 
# cp parameterization. however, it would take time to implement it in Stan. 
# moreover, parameterizing on the cp scale requires a few constraints...
# Betancourt explains how to impose them in the transformed parameters block
# and how it works (rejects if violated)
# https://groups.google.com/g/stan-users/c/fAqBYSf3u6o

# By now, a moderate-asymmetry model could probably fit the data well. Unless the 
# emulated likelihood is truly asymmetric, keep it simple. 
# (it's not that I'll be fitting a symmetric normal.)

# Some functions from sn could be translated to Stan.

getAnywhere("delta.etc")
# it's translatable

getAnywhere("cp2dp")[1]
getAnywhere("cp2dpMv")
getAnywhere("msn.cp2dp")
# it seems understandable for now,
# valid cps are
# abs(gamma1) <  ~0.99
# Omega_cor has to be inversible 

## and this restriction seems to be harder to meet.
# Obar.inv.delta <- as.vector(Obar.inv %*% delta)
# delta.sq <- sum(delta * Obar.inv.delta)
# if (delta.sq >= 1) {
#   if (silent) 
#     return(NULL)
#   else stop("non-admissible CP")
# }

getAnywhere("sn.cp2dp")



# Plotting a few asymmetric curves ----------------------------------------

xmin <- -10
xmax <- 10
curve(skew_normal_lupdf(x, 0, 1, 7), from = xmin, to = xmax,
      ylim = c(-10, 2), n = 500)
abline(h = log(0.1))
curve(skew_normal_lupdf(x, 1, 0.5, 0), add = TRUE, n = 500, col = "red")
curve(skew_normal_lupdf(x, 0.0, 0.8, 4), add = TRUE, n = 500, col = "green")



# Tomorrow: Multivariate versions -----------------------------------------


