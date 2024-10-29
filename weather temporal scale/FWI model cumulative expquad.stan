data {
  int<lower=0> N;
  int<lower=0> N_times;
  vector[N] y;
  matrix[N, N_times] fwi_mat;

  real fwi_center;
  real fwi_scale;
  real fwi_scale_log;

  real prior_sd_int_mu;
  real prior_sd_int_sigma;

  real prior_mean_int_mu;
  real prior_mean_int_sigma;

  real prior_sd_b_mu;
  real prior_sd_b_sigma;

  real prior_sd_ls;
  real time_scale_exp;     // temporal unit for correlation
  real time_scale_expquad; // temporal unit for correlation
}

transformed data {
  matrix[N, N_times] X;
  X = fwi_mat - fwi_center;
}

parameters {
  real alpha_mu_raw;
  real alpha_sigma_raw;

  real<lower=0> b_mu_raw;
  real<lower=0> b_sigma_raw;

  real<lower=0> ls_raw; // lengthscale
}

transformed parameters {
  real alpha_mu = alpha_mu_raw * prior_sd_int_mu + prior_mean_int_mu;
  real alpha_sigma = alpha_sigma_raw * prior_sd_int_sigma + prior_mean_int_sigma;

  real b_mu = b_mu_raw * prior_sd_b_mu / fwi_scale;
  real b_sigma = b_sigma_raw * prior_sd_b_sigma / fwi_scale_log;
  real ls = ls_raw * prior_sd_ls; // lengthscale

  vector[N] mu;
  vector[N] sigma;

  // intercept for non-centered predictors
  real intercept_mu = alpha_mu - b_mu * fwi_center;
  real intercept_sigma = alpha_sigma;

  // weights for previous fwi values
  vector[N_times] weights_improper;
  vector[N_times] weights;
  for(t in 1:N_times)
    weights_improper[t] = exp(-((t-1) / ls) ^ 2);

  // normalize weights
  weights = weights_improper / sum(weights_improper);

  // compute log-mean and log-sigma
  mu = alpha_mu + b_mu * (X * weights);
  sigma = exp(alpha_sigma + b_sigma * log(fwi_mat * weights));
}

model {
  // priors
  alpha_mu_raw ~ std_normal();
  alpha_sigma_raw ~ std_normal();

  b_mu_raw ~ std_normal();
  b_sigma_raw ~ std_normal();

  ls_raw ~ std_normal();

  // likelihood
  y ~ normal(mu, sigma);
}