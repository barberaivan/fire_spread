# Perform ABC Sequential Monte Carlo, using an adaptive procedure and a smooth
# acceptance kernel.

abc_kernel_try <- function(
    overlap, overlap_high, h = 0.1,
    kernel = c(
      "gaussian", "epanechnikov", "biweight", "triweight",
      "cosine", "cauchy", "exp", "smoothstep3", "smoothstep5"
      ),
    log = FALSE
  ) {
  kernel <- match.arg(kernel)
  ov_unit <- pmin(overlap / overlap_high, 1)
  d  <- 1 - ov_unit          # 0..1
  u  <- d / h                 # scaled distance
  
  w <- switch(
    kernel,
    # Infinite support (no hard cutoff); smooth, monotone
    gaussian     = exp(-0.5 * u^2),
    exp          = exp(-u),
    cauchy       = 1 / (1 + u^2),
    
    # Compact support (flat-top); zero weight beyond u>1
    epanechnikov = { z <- pmin(u, 1); (1 - z^2) * (u <= 1) },
    biweight     = { z <- pmin(u, 1); (1 - z^2)^2 * (u <= 1) },
    triweight    = { z <- pmin(u, 1); (1 - z^2)^3 * (u <= 1) },
    cosine       = { z <- pmin(u, 1); cos(pi * z / 2) * (u <= 1) },
    
    # Flat-top "smoothstep" polynomials (C¹ or C² at the join), compact support
    smoothstep3  = { z <- pmin(u, 1); 1 - (3*z^2 - 2*z^3) },                 # C¹
    smoothstep5  = { z <- pmin(u, 1); 1 - (6*z^5 - 15*z^4 + 10*z^3) }       # C²
  )
  
  if (log) return(log(pmax(w, .Machine$double.xmin))) else return(w)
}

# curve(abc_kernel_try(x, 0.6, 0.05, "gaussian"), n = 300)
# curve(abc_kernel_try(x, 0.6, 0.2, "smoothstep5"), n = 300)

# Simply use the gaussian kernel, with flat top.
# overlap: overlap between simulated and observed fires, in [0, 1]
# epsilon: 1 - overlap_high, above which the similarity is considered maximal.
#   Adapt this value in SMC.
# sigma: scale of the kernel.
abc_kernel <- function(overlap, epsilon, sigma = 0.01, log = TRUE) {
  overlap_high <- 1 - epsilon
  ov_unit <- pmin(overlap / overlap_high, 1)
  x <- (1 - ov_unit) / sigma
  ylog <- - 0.5 * x ^ 2
  if (log) ylog else exp(ylog)
}
# curve(abc_kernel(x, 0.1, log = F), n = 300)

# Compute effective sample size (ESS) from weights
ess <- function(weights) {
  w <- weights / sum(weights)   # normalize to sum 1
  return(1 / sum(w^2))
}

# Update the gamma parameter in the DE-MCMC algorithm.
update_gamma <- function(acc_rate, gamma_old, target = 0.234, eta = 0.05) {
  # multiplicative update to move acceptance towards target
  gamma_new <- gamma_old * exp(eta * (acc_rate - target))
  
  # optional: clip to reasonable bounds
  gamma_new <- pmax(gamma_new, 0.01)
  gamma_new <- pmin(gamma_new, 5)
  
  return(gamma_new)
}

# Create proporsals for the Differential Evolution MCMC (ter Braak 2006).
# particles: matrix (N * D) with particles to update.
de_proposal <- function(particles, gamma = NULL, sd_jitter = 0.01) {
  # particles: N x dim matrix (rows = particles, cols = parameters)
  N <- nrow(particles)
  dim_theta <- ncol(particles)
  
  if (is.null(gamma)) {
    gamma <- 2.38 / sqrt(2 * dim_theta)  # default optimal for Gaussian targets
  }
  
  # Prepare proposal matrix
  theta_prop <- matrix(NA, nrow = N, ncol = dim_theta)
  
  for (i in 1:N) {
    # indices for difference, excluding current particle
    candidates <- setdiff(1:N, i)
    a <- sample(candidates, 1)
    b <- sample(setdiff(candidates, a), 1)
    
    # DE-MCMC move with small jitter
    theta_prop[i, ] <- particles[i, ] + gamma * (particles[a, ] - particles[b, ]) +
      rnorm(dim_theta, mean = 0, sd = sd_jitter)
  }
  
  return(theta_prop)
}

# Perform mcmc move, using proposals from DE-MCMC
mcmc_move <- function(
    particles_old, particles_new, overlaps_old,
    epsilon_new, log_prior, simulator
  ) {
  # particles_old: N x dim matrix
  # particles_new: N x dim matrix of proposed particles
  # overlaps_old: vector of current overlaps
  # epsilon_new: new ABC threshold
  # log_prior: function(theta_matrix) returning vector of log priors
  # simulator: function that takes matrix of particles -> vector of overlaps
  # weights: optional, current SMC weights (used for ESS)
  
  N <- nrow(particles_old)
  
  # --- Step 1: Simulate proposed overlaps (vectorized or parallel inside simulator) ---
  overlap_prop <- simulator(particles_new)
  
  # --- Step 2: Compute MH log acceptance ratios ---
  log_acc <- 
    (abc_kernel(overlap_prop, epsilon_new, log = TRUE) +
      log_prior(particles_new)) -
    (abc_kernel(overlaps_old, epsilon_new, log = TRUE) +
       log_prior(particles_old))
  
  # --- Step 3: Accept/reject all at once ---
  u <- log(runif(N))
  accept <- u < log_acc
  
  # --- Step 4: Update particles and overlaps ---
  particles_accept <- particles_old
  particles_accept[accept, ] <- particles_new[accept, ]
  
  overlaps_accept <- overlaps_old
  overlaps_accept[accept] <- overlap_prop[accept]
  
  # --- Step 5: Acceptance probability ---
  acc_prob <- mean(accept)
  
  # --- Step 6: Optional: update weights incrementally ---
  if (!is.null(weights)) {
    log_inc <- abc_kernel(overlaps_accept, epsilon_new, sigma, log = TRUE) -
      abc_kernel(overlaps_old, epsilon_new, sigma, log = TRUE)
    log_w <- log(weights) + log_inc
    log_w <- log_w - max(log_w)
    weights <- exp(log_w)
    weights <- weights / sum(weights)
  }
  
  # --- Step 7: Return list ---
  return(list(
    particles = particles_accept,
    overlaps = overlaps_accept,
    acc_prob = acc_prob,
    weights = weights
  ))
}


get_weights <- function(epsilon_new, epsilon_old, overlaps) {
  # log-weights. As always resample, the computation is simpler, just a quotient
  log_w <- 
    abc_kernel(overlaps, epsilon_new, log = TRUE) -
    abc_kernel(overlaps, epsilon_old, log = TRUE)
  w <- exp(log_w); 
  return(w / sum(w))
}

# Adaptive SMC-ABC sampler. Similar as Del Moral 2011, but with smooth 
# abc kernel, M = 1, and resampling particles in all iterations.
smc_adapt <- function(
    prior_sim,         # function: prior_sim(N) -> list of theta
    simulator,         # function: simulator(theta) -> overlap in (0,1)
    N = 1000,          # number of particles
    n_iter = 100,       # max iterations
    alpha = 0.75,       # ESS reduction factor
    log_prior,         # log density, taking a matrix of particles
    verbose = TRUE
) {
  # --- Initialization ---
  particles <- prior_sim(N) # returns a N * D matrix
  overlaps <- simulator(particles) # runs in parallel over rows
  overlap_max <- max(overlaps)
  epsilon <- 1 # first epsilon, the largest disimilarity possible
  weights <- rep(1/N, N)
  ess <- N
  epsilon_seq <- epsilon
  accept_seq <- 0
  
  posterior_samples <- list()
  
  for (iter in 1:n_iter) {
    if (verbose) cat(sprintf("Iteration %d: epsilon = %.4f, ov_max = %.4f, ESS = %.1f\n",
                             iter, epsilon, overlap_max, ess))
    
    # --- Step 1: Adapt epsilon to hit target ESS and new weights ---
    target_ess <- alpha * N
    f <- function(eps) {
      w <- get_weights(eps, epsilon, overlaps)
      (1 / sum(w^2)) - target_ess
    }
    
    # Find new epsilon by bisection (between 0 and current epsilon)
    lower <- 0
    upper <- epsilon
    while (upper - lower > 1e-6) {
      mid <- (lower + upper) / 2
      if (f(mid) > 0) lower <- mid else upper <- mid
    }
    epsilon_new <- lower
    # check for stability
    if (epsilon_new < 0.01 * epsilon) epsilon_new <- 0.01 * epsilon
    epsilon_seq <- c(epsilon_seq, epsilon_new)    
    
    # New weights 
    weights <- get_weights(epsilon_new, epsilon, overlaps)
    
    # --- Step2: resample according to weights ---
    idx <- sample.int(N, N, replace = TRUE, prob = weights)
    particles <- particles[idx, ]
    overlaps <- overlaps[idx]
    weights <- rep(1/N, N)

    # --- Step 3: MCMC move (DE-MC) ---
    # This step does not require to update weights, as the MCMC leaves the target
    # invariant. It only returns new particles and overlaps.
    
    ## End MCMC
    epsilon <- epsilon_new
    
    posterior_samples[[iter]] <- list(
      particles = particles,
      overlaps = overlaps,
      overlap_high = 1 - epsilon,
      epsilon = epsilon,
      accept = accept_rate
    )
    
    # al final de este loop necesitamos
    # - particle matrix (after move)
    # - overlaps nuevos
  }
  
  out <- list(
    particles = particles,
    weights = weights,
    epsilon_seq = epsilon_seq,
    posterior_samples = posterior_samples
  )
  return(out)
}


# --- Example Usage ---
# Define a toy prior and simulator
prior_sim <- function(n) {
  lapply(1:n, function(x) rnorm(2, mean = 0, sd = 1))  # 2D Gaussian prior
}

simulator <- function(theta) {
  true_theta <- c(1, -1)  # True parameters
  # Overlap = exp(-||theta - true_theta||^2 / 2)
  exp(-sum((theta - true_theta)^2) / 2)
}

# Run SMC
set.seed(123)
result <- adaptive_smc(
  prior_sim = prior_sim,
  simulator = simulator,
  N = 500,
  n_iter = 20,
  sigma = 0.05,
  alpha = 0.9,
  verbose = TRUE
)

# Plot results
plot(result$epsilon_seq, type = 'b', xlab = "Iteration", ylab = "overlap_high",
     main = "Adaptive Tolerance Schedule")