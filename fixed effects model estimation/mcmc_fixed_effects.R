library(microbenchmark)
library(GauPro)
library(adaptMCMC)
library(tidyverse)
library(posterior)

# load gp
wave_last <- readRDS(file.path("files", "pseudolikelihood_estimation",
                               "wave_8_all_fires_gp-fitted-with-more-data.rds"))
gp <- wave_last$gps[[1]]

# wave with more data
w8 <- readRDS(file.path("files", "pseudolikelihood_estimation",
                        "wave_8_all_fires_gp-fitted.rds"))
data_large <- w8$particles_data

# posterior dimension
d <- ncol(data_large$par_values)
par_names <- colnames(data_large$par_values)

# define threshold
like_threshold <- 0.14
ll_threshold <- log(like_threshold)

# get parameter ranges to avoid the GP predicting at regions with no data.
# for betas, the lower limit is not necessary.
data_good <- data_large[data_large$overlap >= like_threshold, ]
nrow(data_good)

lower_int <- quantile(data_good$par_values[, "intercept"], probs = 0.01)
upper_int <- quantile(data_good$par_values[, "intercept"], probs = 0.99)
uppers_log <- apply(data_good$par_values[, -1], 2, quantile, probs = 0.95) %>% log

# x at the original scale
log_prior <- function(x,
                      mu_int = 0, sd_int = 20,
                      r_slope = 0.04,
                      r_fi = 0.15,
                      r_wind = 0.3,
                      r_fwi = 0.15) {
  sum(
    dnorm(x[1], mu_int, sd_int, log = TRUE), # intercept
    dexp(x[2], r_fi, log = TRUE), # vfi
    dexp(x[3], r_fi, log = TRUE), # tfi
    dexp(x[4], r_slope, log = TRUE), # slope
    dexp(x[5], r_wind, log = TRUE), # wind
    dexp(x[6], r_fwi, log = TRUE) # fwi
  )
}

reject_list <- -1e16
# x is at the unconstrained scale
log_posterior <- function(x, mu_int = 0, sd_int = 20,
                          r_slope = 0.04,
                          r_fi = 0.15,
                          r_wind = 0.3,
                          r_fwi = 0.15) {
  # check range of parameter vector
  in_range <- all(x[2:d] <= uppers_log) & # betas (no lower limit)
              x[1] >= lower_int & x[1] <= upper_int # intercept

  if(!in_range) return(reject_list) else {
    # transform to original scale the log-parameters
    x[2:d] <- exp(x[2:d])

    # evaluate likelihood with the GP
    overlap_pred <- gp$pred(x, se.fit = T)
    like <- rnorm(1, overlap_pred$mean, overlap_pred$se)

    # above threshold?
    if(like < like_threshold) return(reject_list) else {

      lprior <- log_prior(x, mu_int, sd_int, r_slope, r_fi, r_wind, r_fwi)
      ll <- log(like)
      return(ll + lprior)
    }
  }
}

# compute covariance from the best particles
set.seed(34234)
ids_samp <- sample(1:nrow(data_good), size = 1000, replace = T,
                   prob = data_good$overlap)

sample_zero <- data_good$par_values[ids_samp, ]
# turn into log the positive parameters
for(i in 2:d) sample_zero[, i] <- log(sample_zero[, i])
s_init <- cov(sample_zero)
# cov2cor(s_init)

# function to get good starting values
data_best <- data_good[data_good$overlap >= 0.15, ] # more restrictive for inits
make_init <- function() {
  id_samp <- sample(1:nrow(data_best), size = 1, prob = data_best$overlap)
  init <- data_best$par_values[id_samp, ]
  init[2:d] <- log(init[2:d])
  return(init)
}

# Run MCMC ----------------------------------------------------------------

sampling_iters <- 50000
adapt_iters <- 50000

rp <- MCMC.parallel(log_posterior,
                    n = sampling_iters + adapt_iters,
                    adapt = adapt_iters,
                    n.cpu = 8,
                    n.chain = 8,
                    scale = s_init, init = make_init(), acc.rate = 0.234,
                    packages = c("GauPro"))
# around 25 min for 100000 samples

# object.size(rp) / 1e6
# saveRDS(rp, file.path("fixed effects model estimation", "draws_01_8chains.rds"))
# saveRDS(rp, file.path("fixed effects model estimation", "draws_02_8chains.rds"))
saveRDS(rp, file.path("fixed effects model estimation", "draws_03_8chains.rds"))

# rp[[1]]$samples
nc <- length(rp)

# extract draws removing adaptation
draws_list <- lapply(1:length(rp), function(i) {
  # i = 1
  mcmc(rp[[i]]$samples[-(1:adapt_iters), ])
})

draws_mcmc <- mcmc.list(draws_list) # from coda
dl <- as_draws_list(draws_mcmc)     # from posterior
sm2 <- posterior::summarise_draws(dl)
saveRDS(sm2, file.path("fixed effects model estimation", "summary_03.rds"))
print(sm2) # terrible, not converging with niter = 100000, adapt = 2000
# not converging either with 1e6 iter, throwing 5e4 for adaptation
# still not converging after getting like_thres = 0.14



# Exploring posterior correlation -----------------------------------------

# the intercept marginal is bimodal. something is happening here.

# thin the chains

# extract draws removing adaptation, and removing draws with
# intercept > 10
s <- rp
thin <- 950000 / 1000
draws_list <- lapply(1:length(s), function(i) {
  # remove adaptation
  d0 <- s[[i]]$samples[-(1:adapt_iters), ]
  # d0 <- d0[d0[, "intercept"] <= 10, ]
  if(nrow(d0) > 1000) {
    thin <- floor(nrow(d0) / 1000)
    # get a sample every 50 iter
    ii <- seq(1, nrow(d0), by = thin)
    return(mcmc(d0[ii, ]))
  } else return(mcmc(d0))
})

largos <- lapply(draws_list, function(x) nrow(x))

draws_mcmc <- mcmc.list(draws_list) # from coda
dl <- as_draws_list(draws_mcmc)     # from posterior
sm2 <- posterior::summarise_draws(dl)


# make matrix
smat <- do.call("rbind", draws_list)
d <- ncol(smat)
# transform parameters to original scale
smat[, 2:d] <- exp(smat[, 2:d])

dim(smat)
GGally::ggpairs(as.data.frame(smat))

for(i in 1:d) {
  plot(density(smat[, i]), main = colnames(smat)[i])
}

plot(smat[, "intercept"] ~ smat[, "wind"], pch = 19, col = rgb(0, 0, 0, 0.1))
plot(smat[, "intercept"] ~ smat[, "vfi"], pch = 19, col = rgb(0, 0, 0, 0.1))
plot(smat[, "intercept"] ~ smat[, "tfi"], pch = 19, col = rgb(0, 0, 0, 0.1))
plot(smat[, "intercept"] ~ smat[, "slope"], pch = 19, col = rgb(0, 0, 0, 0.1))
plot(smat[, "intercept"] ~ smat[, "fwi"], pch = 19, col = rgb(0, 0, 0, 0.1))

# compare data used for setting the loglik bounds with data used to fit the gp


# Vainilla ABC: rejection sampling ----------------------------------------

# define threshold
like_threshold <- 0.14
ll_threshold <- log(like_threshold)

# get parameter ranges to avoid the GP predicting at regions with no data
lower_int <- quantile(data_good$par_values[, "intercept"], probs = 0.01)
upper_int <- quantile(data_good$par_values[, "intercept"], probs = 0.99)
uppers <- apply(data_good$par_values[, -1], 2, quantile, probs = 0.95) # not log

# prior distribution to simulate parameters or to compute them from a sobol
# sequence (type = "quantile", which computes the icdf.)
prior_dist <- function(mu_int = 0, sd_int = 20,
                       r_slope = 0.04,
                       r_fi = 0.15,
                       r_wind = 0.3,
                       r_fwi = 0.15,
                       type = "rng", # or "quantile"
                       n = 1,
                       p = NULL) {

  if(type == "rng") {
    b <- cbind(
      "intercept" = rnorm(n, mu_int, sd_int),
      "vfi" = rexp(n, r_fi),
      "tfi" = rexp(n, r_fi),
      "slope" = rexp(n, r_slope),
      "wind" = rexp(n, r_wind),
      "fwi" = rexp(n, r_fwi)
    )

    if(n == 1) {
      parnames <- colnames(b)
      b <- as.numeric(b)
      names(b) <- parnames
    }
    return(b)
  }

  if(type == "quantile") {

    if(is.null(p)) stop("Provide probability to compute quantiles.")
    if(length(p) %% 5 != 0) stop("p must be a multiple of the number of parameters (5).")

    if(length(dim(p)) < 2) p <- matrix(p, ncol = 6)
    if(length(dim(p)) == 2) {
      if(ncol(p) != 6) p <- matrix(as.numeric(p), ncol = 6)
    }

    q <- cbind(
      "intercept" = qnorm(p[, 1], mu_int, sd_int),
      "vfi" = qexp(p[, 2], r_fi),
      "tfi" = qexp(p[, 3], r_fi),
      "slope" = qexp(p[, 4], r_slope),
      "wind" = qexp(p[, 5], r_wind),
      "fwi" = qexp(p[, 6], r_fwi)
    )

    return(q)
  }
}

# matrix to fill with particles passing the threshold:
particles_box <- matrix(NA, 1e5, d)

got <- 0
cycle <- 1
while(got < 1e4) {
  print(paste("cycle", cycle))
  print(paste("got", got))

  # simulate from the prior
  try <- prior_dist(n = 1000, type = "rng")

  # check range of parameters
  betas_ok <- apply(try[, 2:d], 1, function(x) all(x < uppers))
  int_ok <- try[, 1] >= lower_int & try[, 1] <= upper_int
  keep <- betas_ok & int_ok
  try2 <- try[keep, ]
  n_pass <- sum(keep)

  if(n_pass > 0) {
    overlap_pred <- gp$pred(try2, se.fit = T)
    like <- rnorm(nrow(overlap_pred), overlap_pred$mean, overlap_pred$se)
    # like <- rep(0.04, n_pass)
    good <- like >= like_threshold
    kept <- sum(good)

    if(kept > 0) {
      particles_box[(got+1) : (got+kept), ] <- try2[good, ]
      got <- got + kept
    }
  }
  cycle <- cycle + 1
}

bb <- particles_box[complete.cases(particles_box), ]
dim(bb)
colnames(bb) <- par_names

GGally::ggpairs(as.data.frame(bb))

saveRDS(bb,  file.path("fixed effects model estimation", "vainilla_abc_samples.rds"))

# Notes -------------------------------------------------------------------

# when using like_threshold = 0.1 the GP predicted high loglik for extreme
# intercept values. This did not make sense, because the GP had a huge uncertainty
# there, where the intercepts where probably affected by the bound imposed by
# landscape size.

# I run 1e6 iterations for 8 chains, and could not find convergenge due to a
# small mode for high intercept values (run 02).

# Try to overcome this issue in run 03, where like_thres = 0.15. By bounding
# the parameters to the extreme values found at overlap >= 0.15, the upper limit
# on the intercept y far more restrictive.

plot(data_large$overlap ~ data_large$par_values[, "intercept"],
     col = rgb(0, 0, 0, 0.2), pch = 19)
abline(v = quantile(data_large$par_values[, "intercept"], prob = c(0.02, 0.98)),
       lty = 2, col = "red")
abline(v = quantile(data_large$par_values[data_large$overlap >= 0.15, "intercept"],
                    prob = c(0, 1)),
       h = 0.15,
       lty = 2, col = "blue")
abline(v = quantile(data_large$par_values[data_large$overlap >= 0.14, "intercept"],
                    prob = c(0.01, 0.99)),
       h = 0.14,
       lty = 2, col = "green")

sum(data_large$overlap >= 0.14) # 482 data points


plot(data_large$overlap ~ data_large$par_values[, "slope"],
     col = rgb(0, 0, 0, 0.2), pch = 19)
abline(v = quantile(data_large$par_values[, "slope"], prob = c(0.02, 0.98)),
       lty = 2, col = "red")
abline(v = quantile(data_large$par_values[data_large$overlap >= 0.15, "slope"],
                    prob = c(0, 1)),
       h = 0.15,
       lty = 2, col = "blue")
abline(v = quantile(data_large$par_values[data_large$overlap >= 0.14, "slope"],
                    prob = c(0.01, 0.99)),
       h = 0.14,
       lty = 2, col = "green")



# Al final hice vainilla abc porque muestrear no funcionaba.