# the vegetation parameters can be parameterized in two main ways:
# with an intercept term, where one class is used as the reference, and
# without intercept, where there is no reference class.
# The advantage of the reference is that a random effect can be centred on it,
# but it may have problems if in certain fires the reference vegetation class is
# absent. In addition, the intercept parameterization might increase the 
# correlation between parameters.
# Here I simulate data and fit the model in two ways, to test whether the 
# intercept_free parameterization is preferable.

library(rstan)
library(tidyverse)

# stan models

# intercept-free
model_free <- "
data {
  int N;
  int K;
  vector[N] y;
  matrix[N, K] X_free;
}

parameters {
  vector[K] b_free;
  real<lower = 0> sigma;
}

transformed parameters {
  vector[K] b_int = b_free;
  b_int[2:K] = b_int[2:K] - b_free[1];
}

model {
  y ~ normal(X_free * b_free, sigma);
  b_free ~ normal(0, 30);
  sigma ~ normal(0, 5);
}
"
smodel_free <- stan_model(model_code = model_free)

# intercept model
model_int <- "
data {
  int N;
  int K;
  vector[N] y;
  matrix[N, K] X_int;
}

parameters {
  vector[K] b_int;
  real<lower = 0> sigma;
}

transformed parameters {
  vector[K] b_free = b_int;
  b_free[2:K] = b_int[2:K] + b_int[1];
}

model {
  y ~ normal(X_int * b_int, sigma);
  b_int ~ normal(0, 30);
  sigma ~ normal(0, 5);
}
"
smodel_int <- stan_model(model_code = model_int)


# data

K <- 4
b_free <- rnorm(K, 0, 10)
b_int <- b_free
b_int[2:K] <- b_free[2:K] - b_free[1]

# enes <- rep(150, K)
enes <- c(1, rep(200, K-1))

N <- sum(enes)
x_cat <- factor(rep(letters[1:K], enes))
X_free <- model.matrix(~ x_cat - 1)
X_int <- model.matrix(~ x_cat)
mu <- X_free %*% b_free
y <- rnorm(N, mu, 1)
plot(y ~ x_cat)

sdata <- list(y = y, X_free = X_free, X_int = X_int, N = N, K = K)
# remove row from the reference level
sdata <- list(y = y[-1], X_free = X_free[-1, ], X_int = X_int[-1, ], 
              N = N-1, K = K)

# fit model

m_free <- sampling(smodel_free, data = sdata, chains = 3, cores = 1)
m_int <- sampling(smodel_int, data = sdata, chains = 3, cores = 1)

pairs(m_free, pars = "b_free")
pairs(m_free, pars = "b_int")

pairs(m_int, pars = "b_free")
pairs(m_int, pars = "b_int")

# the intercept parameterization shows high posterior correlation.
# Even when there is no data for a given level, the free parameterization works
# well, copying the prior. But when the no-data level is taken as reference, 
# sampling is troublesome, and the absolute posterior correlation approaches 1. 