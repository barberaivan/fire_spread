// Linear model using only marginal gaussian bases, each dimension having its
// own smoothing parameter.

data {
  int<lower=0> N;
  int<lower=0> Kd; // number of knots by dimension
  int<lower=0> D; // number of dimensions
  vector[N] y;
  matrix[N, Kd] X[D]; // array of matrices
  real y_lower;   // truncation limit
  real prior_sigma_b_sd;
}

parameters {
  matrix<lower=0>[Kd, D] b_raw;
  real<lower=0> sigma;
  vector<lower=0>[D] sigma_b;
}

transformed parameters {
  matrix<lower=0>[Kd, D] b;
  matrix[N, D] smooths;
  vector[N] mu;

  for(d in 1:D) {
    b[, d] = b_raw[, d] * sigma_b[d];
    smooths[, d] = X[d] * b[, d];    // X[d] subsets the bases matrix for the dth dimension
  }

  // mu = apply(smooths, 1, prod)
  mu = smooths[, 1];
  for(d in 2:D)
    mu = mu .* smooths[, d];
}

model {
  sigma ~ std_normal();
  sigma_b ~ normal(0, prior_sigma_b_sd);
  to_vector(b_raw) ~ std_normal();
  for(n in 1:N)
    y[n] ~ normal(mu[n], sigma)T[y_lower, ];
}