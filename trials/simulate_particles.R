library(randtoolbox)   # sobol sequences

# Prior simulator function. It simulates samples from the prior in a
# space-filling manner. Samples are computed from cumulative probabilities
# obtained from Sobol sequences. Then, the quantiles associated to these
# probabilities in the prior distribution are computed.

# p is a matrix of n_particles * d sobol samples of cumulative probability.
# d is the dimension of the parameter vector.

prior_q <- function(p, mu_int = 0, sd_int = 20, sd_veg = 5,
                    r_01 = 0.05, r_z = 0.15) {

  q <- matrix(NA, nrow(p), ncol(p))
  colnames(q) <- c("intercept", "subalpine", "wet", "dry", "fwi",
                   "aspect", "wind", "elevation", "slope")
  colnames(p) <- colnames(q)

  # fill matrix with parameters samples

  # intercept
  q[, 1] <- qnorm(p[, 1], mean = mu_int, sd = sd_int)

  # vegetation parameters
  for (i in 2:4) {
    q[, i] <- qnorm(p[, i], mean = 0, sd = sd_veg)
  }

  # parameters for [-1, 1] variables
  names_01 <- which(colnames(q) %in% c("aspect", "wind", "slope"))
  for(i in names_01) {
    q[, i] <- qexp(p[, i], rate = r_01)
  }

  # parameters for standardized variables
  names_z <- which(colnames(q) %in% c("fwi", "elevation"))
  for(i in names_z) {
    q[, i] <- qexp(p[, i], rate = r_z)
  }

  return(q)
}

# Example -----------------------------------------------------------------

p <- sobol(200, d = 9)
pairs(p, col = rgb(0, 0, 0, 0.1), pch = 19) # note the space-filling property
                                            # in [0, 1]

params <- prior_q(p)
pairs(params, col = rgb(0, 0, 0, 0.1), pch = 19) # space-filling in the prior
                                                 # distribution space