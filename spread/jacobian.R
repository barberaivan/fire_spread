# sample steps_log, specifying a prior at steps_logit = logit(exp(steps_log))
library(tidyverse)
library(rstan)

U <- 300
L <- 5
mu <- 1
sigma <- 1
N <- 1e5
steps_logit <- rnorm(N, mu, sigma)

steps <- plogis(steps_logit) * (U - L) + L # in [10, 11]
steps_log <- log(steps)



initf <- function() {
  return(list(steps_log = runif(1, L, U) |> log()))
}

sdata <- list(L = L, U = U, mu = mu, sigma = sigma)
smodel <- stan_model("hierarchical model/jacobian.stan")

m1 <- sampling(smodel, data = sdata, init = initf,
               warmup = 10000, iter = 20000,
               control = list(adapt_delta = 0.999))

xx <- as.matrix(m1)


par(mfrow = c(1, 3))
plot(density(steps_logit), main = "steps_logit")
lines(density(xx[, "steps_logit"]), col = "blue")
lines(density(xx[, "steps_logit2"]), col = "red")

plot(density(steps, from = L, to = U),
     xlim = c(L, U), main = "steps")
lines(density(xx[, "steps"], from = L, to = U),
      col = "blue")

plot(density(steps_log, from = log(L), to = log(U)),
     xlim = log(c(L, U)), main = "steps_log")
lines(density(xx[, "steps_log"], from = log(L), to = log(U)),
      col = "blue")
par(mfrow = c(1, 1))
# Bien! anduvo


