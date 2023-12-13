library(adaptMCMC)

# log-densities to sample from
f_certain <- function(x) {
  like <- dnorm(x)
  if(like < dnorm(0) * 0.2) return(-Inf) else return(log(like))
}

# vary sigma with mu
sigma_fun <- function(x) 0.05 * plogis(-4 * x)

f_uncertain <- function(x) {
  if(x < -5 | x > 5) return(-Inf)
  like <- rnorm(length(x), dnorm(x), sigma_fun(x))
  if(like < dnorm(0) * 0.2) return(-Inf) else return(log(like))
  return()
}

# MCMC
sampling_iters <- 1e5
adapt_iters <- 1e4
initial <- rnorm(1)

draws_certain <- MCMC(f_certain,
                      n = sampling_iters + adapt_iters,
                      adapt = adapt_iters,
                      scale = 0.3, init = initial, acc.rate = 0.234)
draws_uncertain <- MCMC(f_uncertain,
                        n = sampling_iters + adapt_iters,
                        adapt = adapt_iters,
                        scale = 0.3, init = initial, acc.rate = 0.234)

# Compare
xseq <- seq(-5, 5, length.out = 150)
upper <- max(qnorm(0.975, dnorm(0), sigma_fun(dnorm(0))))
curve(dnorm(x), to = 5, from = -5, lwd = 2, ylim = c(0, upper),
      ylab = "Density", xlab = "X")
abline(h = dnorm(0) * 0.2, col = "red", lty = 3)
curve(qnorm(0.975, dnorm(x), sd = sigma_fun(dnorm(x))), add = T, lty = 2)
curve(qnorm(0.025, dnorm(x), sd = sigma_fun(dnorm(x))), add = T, lty = 2)
lines(density(draws_certain$samples[-(1:adapt_iters), ]), col = "blue", lty = 2)
lines(density(draws_uncertain$samples[-(1:adapt_iters), ]), col = "red", lty = 2)

plot(density(draws_uncertain$samples[-(1:adapt_iters), ]), col = "red", lty = 2)

text(2, upper * 0.8, "True density")
text(-3, upper * 0.5, "Certain sampling", col = "blue")
text(3, upper * 0.5, "Uncertain sampling", col = "red")

# Ahora no funciona pero por problemas con la inequidad de Jensen.
# Qué definimos como la verdadera? La media o la mediana en la escala exp?
# Creo que si defino como verdad la media en exp, hay que corregir cuál
# es esa certain density, y ahí todo coincidiría.

# Simpler explanation ------------------------------------------------------

# Assuming that we have uncertainty about a density function that we want to
# sample from we may, for simplicity add normal noise to the log-density.
# If we are using MCMC, in every step we should sample a normal shift to add to
# the log-density. However, an unnormalized density is invariant to addition at
# the log scale. Hence, adding a constant random noise to the density is
# irrelevant.
# Certainly, it does not apply for a heterogeneous noise at the log scale.



# Como sea, da igual muestrear de una o de otra. Y con el gp, calcular el SD
# es mucho más costoso, y re cuestra muestrear las likelihoods.