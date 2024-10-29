# Log-uniform density

# Here,
# https://stats.stackexchange.com/questions/456136/how-to-estimate-the-pdf-of-the-logarithm-of-a-uniformly-distributed-random-varia# It's said that
# X ~ Unif(a, b) for 0 < a < b
# then
# Y = log(X) ~ log-uniform(a, b).

# But wikipedia calls Y ~ log-uniform() when
# X = log(Y) ~ uniform:
# https://en.wikipedia.org/wiki/Reciprocal_distribution

# With the steps parameter prior the first case is implied, but
# for consistency with the log-normal and logit-normal nomenclature,
# I call it exp-uniform:
# If X ~ Uniform(a, b) # steps estimated in the first stage
# and Y = log(X)
# Y ~ Exp-uniform(a, b) # log-steps, this is the prior for the second stage.

dexpunif <- function(x, l = 1, u = 10) { # for x \in [log(l), log(u)]
  d <- ifelse(x >= log(l) & x <= log(u),
              exp(x) / (u-l), 0)
  return(d)
}
curve(dexpunif(x, 2, 6), from = log(2-1), to = log(6+1))

yy <- runif(1e6, 2, 6)
lines(density(log(yy), from = log(2), to = log(6)), col = 2)
