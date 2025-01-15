data {
  int<lower=0> N;
  int<lower=0> N_times;
  vector[N] y;
  matrix[N, N_times] fwi_mat;

  real prior_mean_a; // model for mu
  real prior_sd_a;
  real prior_sd_b;

  real prior_mean_c; // model for sigma
  real prior_sd_c;
  real prior_sd_d;

  real prior_sd_ls;
}

parameters {
  real a_raw;
  real<lower=0> b_raw;

  real c_raw;
  real<lower=0> d_raw;

  real<lower=0> ls_raw; // lengthscale
}

transformed parameters {
  real a = a_raw * prior_sd_a + prior_mean_a;
  real b = b_raw * prior_sd_b;

  real c = c_raw * prior_sd_c + prior_mean_c;
  real d = d_raw * prior_sd_d;

  real ls = ls_raw * prior_sd_ls; // lengthscale

  vector[N] mu;
  vector[N] s;

  // weights for previous fwi values
  vector[N_times] weights_improper;
  vector[N_times] weights;
  for(t in 1:N_times)
    weights_improper[t] = exp(-((t-1) / ls) ^ 2);

  // normalize weights
  weights = weights_improper / sum(weights_improper);

  // mu
  mu = a + b * (fwi_mat * weights);
  s = exp(c + d * fwi_mat * weights);
}

model {
  // priors
  a_raw ~ std_normal();
  b_raw ~ std_normal();

  c_raw ~ std_normal();
  d_raw ~ std_normal();

  ls_raw ~ std_normal();

  // likelihood
  y ~ normal(mu, s);
}