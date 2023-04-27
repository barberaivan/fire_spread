library(randtoolbox)
?quasiRNG

# make a 2D sobol sequence in [0, 1]
N <- 500
K <- 2
s1 <- sobol(N, dim = K)
plot(s1[, 2] ~ s1[, 1], col = rgb(0, 0, 0, 0.3), pch = 19)

# add more points
s2 <- sobol(N, dim = K, init = F)
points(s2[, 2] ~ s1[, 1], col = rgb(1, 0, 0, 0.3), pch = 19)
# if init = F, the sobol sequence is augmented, so you keep filling the space.

# The sobol sequence in [0, 1] is useful to cover the prior distribution space.
# Assume you have 2 parameters, with prior distributions
# b1 ~ normal(0, 1)
# b2 ~ gamma(1, 1) # only positive values
# Then, using the inverse cumulative disitribution functions you can map from
# [0, 1] to quantiles:

# transform the first sequence
r1 <- s1
r1[, 1] <- qnorm(s1[, 1], 0 , 1) # q stands for quantile, the inverse-cumulative
                                 # distribution function
r1[, 2] <- qgamma(s1[, 2], 1, 1)
plot(r1[, 2] ~ r1[, 1], xlab = "b1", ylab = "b2",
     col = rgb(0, 0, 0, 0.1), pch = 19)

# transform the second sequence
r2 <- s2
r2[, 1] <- qnorm(s2[, 1], 0 , 1) # q stands for quantile, the inverse-cumulative
# distribution function
r2[, 2] <- qgamma(s2[, 2], 1, 1)
points(r2[, 2] ~ r2[, 1],
       col = rgb(1, 0, 0, 0.1), pch = 19)

# This plot describes the prior (bivariate) distribution for b1 and b2. Note that
# the points density is proportional to the joint probability.

# In Bayesian analyses, the prior is a soft restriction on the parameter space,
# so we usually are not interested in the model behaviour outside the prior.
# This is why the particles are thrown in the prior space.