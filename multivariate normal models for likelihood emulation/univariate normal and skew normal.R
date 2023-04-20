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

# try a nice shape
xi <- 0; sigma <- 1; alpha <- 7

N <- 10000
x <- runif(N, -5, 5) # hasta 5 antes. hasta 2 estima re bien porque no genera mues tan bajos
mu <- skew_normal_lupdf(x, xi, sigma, alpha)  

# use only data with mu > log(0.1)
keep <- which(mu > log(0.1))
mu <- mu[keep]
x <- x[keep]
N <- length(mu)

# simulate tau, the residual standard deviation based on sd(mu) to get
# an R2 near 0.8
tau <- sd(mu) * 0.5
y <- rnorm(N, mu, tau)
(R2 <- var(mu) / (var(mu) + var(y - mu)))

plot(y ~ x, col = rgb(0, 0, 0, 0.2), pch = 19)
curve(skew_normal_lupdf(x, xi = xi, sigma = sigma, alpha), add = TRUE, lwd = 1.8,
      col = "blue")

sdata2 <- list(
  # data
  x = x, y = y, N = N,
  
  # prior parameters
  xi_prior_sd = 10,
  sigma_prior_sd = 20,#sigma_prior_sd,
  sigma_prior_nu = 200,#sigma_prior_nu,
  tau_prior_sd = 20,
  tau_prior_nu = 200,
  alpha_prior_sd = 20,
  
  # true parameter values
  xi_true = xi,
  sigma_true = sigma,
  alpha_true = alpha,
  tau_true = tau
)
m2 <- sampling(smodel2, data = sdata2, seed = 2126, 
               cores = 2, chains = 2, iter = 2000, warmup = 1000)
pairs(m2, pars = c("xi", "sigma", "alpha", "tau"))

# check mu is computed well.
mu_etc_stan_raw <- as.matrix(m2, pars = "mu_true_etc")[1, ] # they are all the same
mu_etc_stan <- matrix(mu_etc_stan_raw, ncol = 10)
colnames(mu_etc_stan) <- c("log_den", "log_cum", "mu",
                           "x_centred", "x_centred2", "prec",
                           "shifted", "erf", "phi",
                           "mu_phi")
# plot(mu_stan ~ mu); abline(0, 1)
test_mat <- cbind(mu_etc_stan, mu_r = mu)
View(test_mat)
plot(test_mat[, "mu_phi"] ~ test_mat[, "mu_r"]); abline(0, 1)
all.equal(test_mat[, "mu_phi"], test_mat[, "mu_r"])
# Now mu is well computed. The problem was on the imprecision of stan::erf(),
# and it was resolved with stan::Phi()

# Plot fitted mu
samples2 <- extract(m2, pars = c("xi", "sigma", "alpha"))

plot(y ~ x, col = rgb(0, 0, 0, 0.2), pch = 19)
curve(skew_normal_lupdf(x, xi = xi, sigma = sigma, alpha), add = TRUE, lwd = 1.8,
      col = "blue")
for(i in 1:500) {
  id <- sample(1:length(samples2$alpha), 1)
  curve(skew_normal_lupdf(x, samples2$xi[id], samples2$sigma[id], samples2$alpha[id]), 
        add = TRUE, lwd = 0.5, col = rgb(1, 0, 0, 0.05))
}

# Compare parameters
head(sm2 <- summary(m2, pars = c("xi", "sigma", "alpha", "tau"))[[1]])
summ <- as.data.frame(sm2[, c("mean", "2.5%", "97.5%")])
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
  y_ <- rnorm(N, mu, samples2$tau[i])
}) # OK, columns are replicates
res <- createDHARMa(y_sim, y)
plot(res)


# fitted mu and real mu
plot(rowMeans(y_sim) ~ mu, pch = 19, col = rgb(0, 0, 0.8, 0.05),
     ylab = "Estimated mu", xlab = "Real mu")
abline(0, 1)

# Extras

# mu computed in R and Stan
mu_stan <- as.matrix(m2, "mu") %>% t
all.equal(mu_stan[, 1], 
          skew_normal_lupdf(x, samples2$xi[1], samples2$sigma[1],
                             samples2$alpha[1]))
j = sample(1:2000, 1)
plot(mu_stan[, j] ~ 
     skew_normal_lupdf(x, samples2$xi[j], samples2$sigma[j], samples2$alpha[j]))
abline(0, 1)
# not exactly the same, but almost. It shouldn't be a problem. Maybe it's caused
# by contrasting sensitivity in decimals between R and Stan.


# Notes on univariate skew-normal -----------------------------------------

# All works well there have been a few cases with a high asymmetry and 
# constraining mu > log(0.1) where the posterior was bimodal. One mode fitted
# perfectly, and the other fitted a symmetric curve that wasn't really so bad. 
# A case like that could occur, but a closer look at the posterior could help 
# to solve the problem.