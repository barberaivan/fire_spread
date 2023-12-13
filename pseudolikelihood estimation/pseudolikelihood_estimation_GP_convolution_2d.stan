// shinking priors for marginal and bivariate basis separately, but not more
// than that.

data {
  int<lower=0> N;
  int<lower=0> P_marg; // number of knots
  int<lower=0> P_int; // number of knots
  vector[N] y;
  matrix[N, P_marg] X_marg;
  matrix[N, P_int] X_int;
  real y_lower; // truncation limit
  real prior_sigma_b_marg_sd;
  real prior_sigma_b_int_sd;
}

parameters {
  vector<lower=0>[P_marg] b_marg_raw;
  vector<lower=0>[P_int] b_int_raw;
  real<lower=0> sigma;
  real<lower=0> sigma_b_marg;
  real<lower=0> sigma_b_int;
}

transformed parameters {
  vector<lower=0>[P_marg] b_marg = b_marg_raw * sigma_b_marg;
  vector<lower=0>[P_int] b_int = b_int_raw * sigma_b_int;
  vector[N] mu = X_marg * b_marg + X_int * b_int;
}

model {
  sigma ~ std_normal();
  sigma_b_marg ~ normal(0, prior_sigma_b_marg_sd);
  sigma_b_int ~ normal(0, prior_sigma_b_int_sd);
  b_marg_raw ~ std_normal();
  b_int_raw ~ std_normal();

  for(n in 1:N)
    y[n] ~ normal(mu[n], sigma)T[y_lower, ];
}