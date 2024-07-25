// Linear model using only marginal gaussian bases, each dimension having its
// own smoothing parameter.

data {
  int<lower=0> N;
  int<lower=0> K; // number of knots (total)
  int<lower=0> D; // number of dimensions
  vector[N] y;
  matrix[N, K] X;
  matrix[K, D] S; // design matrix to multiply each bunch of coefs by their
                  // corresponding variance
  real y_lower;   // truncation limit
  real prior_sigma_b_sd;
}

parameters {
  vector<lower=0>[K] b_raw;
  real<lower=0> sigma;
  vector<lower=0>[D] sigma_b;
}

transformed parameters {
  vector<lower=0>[K] sigma_b_rep = S * sigma_b; // repeats each sigma K/D times
  vector<lower=0>[K] b = b_raw .* sigma_b_rep;
  vector[N] mu = X * b;
}

model {
  sigma ~ std_normal();
  sigma_b ~ normal(0, prior_sigma_b_sd);
  b_raw ~ std_normal();
  for(n in 1:N)
    y[n] ~ normal(mu[n], sigma)T[y_lower, ];
}