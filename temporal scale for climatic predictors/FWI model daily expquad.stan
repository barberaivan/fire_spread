data {
  int<lower=0> N;
  int<lower=0> N_days;
  vector[N] y;
  matrix[N, N_days] fwi_previous;

  // reference correlation values to compute centered cumulative term
  real rho_reference;
  real ls_reference;

  real prior_sd_int;
  real prior_mean_int;
  real prior_sd_b;
  real prior_sd_ls;
  real prior_sd_sigma;
}

transformed data {

  /*
  days_seq <- 0:(120-1)
fwi_previous_center_exp <- mean(fwi_previous[rows_fit, ] %*%
                                  normalize(0.3 ^ (days_seq / 15)))
fwi_previous_center_expquad <- mean(fwi_previous[rows_fit, ] %*%
                                      normalize(exp(-(days_seq / 20) ^ 2)))
  */
  real fwi_center;
  real fwi_scale;
  matrix[N, N_days] X;

  vector[N_days] weights_improper_ref;
  vector[N_days] weights_ref;
  vector[N] fwi_ref;

  // define center and scale for previous period, using reference correlation
  // parameter
  for(d in 1:N_days)
    weights_improper_ref[d] = exp(-((d-1) / ls_reference) ^ 2);
  // normalize weights
  weights_ref = weights_improper_ref / sum(weights_improper_ref);
  fwi_ref = fwi_previous * weights_ref;

  fwi_center = mean(fwi_ref);
  X = fwi_previous - fwi_center;
}

parameters {
  real alpha_raw;
  real<lower=0> b_focal_raw;
  real<lower=0> b_previous_raw;
  real<lower=0> ls_raw; // lengthscale
  real<lower=0> sigma_raw;
}

transformed parameters {
  real alpha_mu = alpha_raw * prior_sd_int + prior_mean_int;
  real b_mu = b_focal_raw * prior_sd_b / fwi_scale;
  real b_previous = b_previous_raw * prior_sd_b / fwi_scale;
  real ls = ls_raw * prior_sd_ls; // lengthscale
  real sigma = sigma_raw * prior_sd_sigma;

  vector[N] mu; // to check if we did it right
  vector[N] mu2; // to check if we did it right

  // intercept for non-centered predictors
  real intercept = alpha - (b_focal * fwi_focal_center +
                            b_previous * fwi_previous_center);

  // weights for previous fwi values
  vector[N_days] weights_improper;
  vector[N_days] weights;
  for(d in 1:N_days)
    weights_improper[d] = exp(-((d-1) / ls) ^ 2);

  // normalize weights
  weights = weights_improper / sum(weights_improper);

  // compute log-mean
  mu = alpha + b_focal * x_focal +
               b_previous * (x_previous * weights);

  // did we center well??
  mu2 = intercept + b_focal * fwi_focal +
                    b_previous * (fwi_previous * weights);
}

model {
  // priors
  alpha_raw ~ std_normal();
  b_focal_raw ~ std_normal();
  b_previous_raw ~ std_normal();
  ls_raw ~ std_normal();
  sigma_raw ~ std_normal();

  // likelihood
  y ~ normal(mu, sigma);
}