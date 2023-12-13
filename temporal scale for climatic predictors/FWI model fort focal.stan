data {
  int<lower=0> N;
  vector[N] y;
  vector[N] fwi_fort;

  real prior_sd_int_mu;
  real prior_sd_int_sigma;

  real prior_mean_int_mu;
  real prior_mean_int_sigma;

  real prior_sd_b_mu;
  real prior_sd_b_sigma;
}

transformed data {
  vector[N] x_fort = fwi_fort - mean(fwi_fort);
  real fwi_center = mean(fwi_fort);
  real fwi_scale = sd(fwi_fort);
  real fwi_scale_log = sd(log(fwi_fort));
}

parameters {
  real alpha_mu_raw;
  real alpha_sigma_raw;

  real<lower=0> b_mu_raw;
  real<lower=0> b_sigma_raw;
}

transformed parameters {
  real alpha_mu = alpha_mu_raw * prior_sd_int_mu + prior_mean_int_mu;
  real alpha_sigma = alpha_sigma_raw * prior_sd_int_sigma + prior_mean_int_sigma;

  real b_mu = b_mu_raw * prior_sd_b_mu / fwi_scale;
  real b_sigma = b_sigma_raw * prior_sd_b_sigma / fwi_scale_log;

  vector[N] mu;
  vector[N] sigma;

  // intercept for non-centered predictors
  real intercept_mu = alpha_mu - b_mu * fwi_center;
  real intercept_sigma = alpha_sigma;

  // compute log-mean and log-sigma
  mu = alpha_mu + b_mu * x_fort;
  sigma = exp(alpha_sigma + b_sigma * log(fwi_fort));
}

model {
  // priors
  alpha_mu_raw ~ std_normal();
  alpha_sigma_raw ~ std_normal();

  b_mu_raw ~ std_normal();
  b_sigma_raw ~ std_normal();

  // likelihood
  y ~ normal(mu, sigma);
}