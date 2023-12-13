data {
  int<lower=0> N;
  int<lower=0> N_days;
  vector[N] y;
  matrix[N, N_days] fwi_previous;
  real fwi_center;
  real fwi_scale;

  real prior_sd_int;
  real prior_mean_int;
  real prior_sd_b;
  real prior_sd_ls;
  real prior_sd_sigma;
}

transformed data {
  matrix[N, N_days] x_previous;
  x_previous = fwi_previous - fwi_center;
}

parameters {
  real alpha_raw;
  real<lower=0> b_raw;
  real<lower=0> ls_raw; // lengthscale
  real<lower=0> sigma_raw;
}

transformed parameters {
  real alpha = alpha_raw * prior_sd_int + prior_mean_int;
  real b = b_raw * prior_sd_b / fwi_scale;
  real ls = ls_raw * prior_sd_ls; // lengthscale
  real sigma = sigma_raw * prior_sd_sigma;

  vector[N] mu; // to check if we did it right

  // intercept for non-centered predictors
  real intercept = alpha - b * fwi_center;

  // weights for previous fwi values
  vector[N_days] weights_improper;
  vector[N_days] weights;
  for(d in 1:N_days)
    weights_improper[d] = exp(-((d-1) / ls) ^ 2);

  // normalize weights
  weights = weights_improper / sum(weights_improper);

  // compute log-mean
  mu = alpha + b * (x_previous * weights);
}

model {
  // priors
  alpha_raw ~ std_normal();
  b_raw ~ std_normal();
  ls_raw ~ std_normal();
  sigma_raw ~ std_normal();

  // likelihood
  y ~ normal(mu, sigma);
}