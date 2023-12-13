// the simplest GP by convolution, with pre-computed weights, independent,
// positive coefficients and truncated normal likelihood.

data {
  int<lower=0> N;
  int<lower=0> K; // number of knots
  vector[N] y;
  matrix[N, K] X;
  real y_lower; // truncation limit
  real prior_sigma_b_sd;
}

parameters {
  vector<lower=0>[K] b_raw;
  real<lower=0> sigma;
  real<lower=0> sigma_b;
}

transformed parameters {
  vector<lower=0>[K] b = b_raw * sigma_b;
  vector[N] mu = X * b;
}

model {
  sigma ~ std_normal();
  sigma_b ~ normal(0, prior_sigma_b_sd);
  b_raw ~ std_normal();
  for(n in 1:N)
    y[n] ~ normal(mu[n], sigma)T[y_lower, ];
}