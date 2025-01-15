data {
  int n;        // number of ignitions
  int nlag;     // lag length for FWI
  int n_notna;
  int n_na;

  vector[n] y;
  int ids_notna[n_notna];
  int ids_na[n_na];
  real cutoff;

  matrix[n, nlag] fwi_mat;
  vector[n] vfi;
  vector[n] tfi;
  vector[n] drz;
  vector[n] dhz;

  real prior_a_mean;
  real prior_a_sd;
  real prior_b_sd;
  real prior_sigma_sd;
  real prior_g_sd;
  real prior_ls_sd;
}

parameters {
  real a_raw;

  real<lower=0> b_fwi_raw;
  real<lower=0> b_vfi_raw;
  real<lower=0> b_tfi_raw;
  real b_drz_raw;
  real b_dhz_raw;

  real<lower=0> sigma_raw;
  real g_raw;
  real<lower=0> ls_raw;
}

transformed parameters {
  real a = a_raw * prior_a_sd;

  real b_fwi = b_fwi_raw * prior_b_sd;
  real b_vfi = b_vfi_raw * prior_b_sd;
  real b_tfi = b_tfi_raw * prior_b_sd;
  real b_drz = b_drz_raw * prior_b_sd;
  real b_dhz = b_dhz_raw * prior_b_sd;

  real sigma = sigma_raw * prior_sigma_sd;
  real g = g_raw * prior_g_sd;
  real ls = ls_raw * prior_ls_sd;                      // FWI temporal corr

  vector[nlag] lag_weights_un; // unnormalized
  vector[nlag] lag_weights;    // time-weights for FWI

  vector[n] mu;

  for(t in 1:nlag)
    lag_weights_un[t] = exp(-((t-1) / ls) ^ 2); // t-1 to start at zero
  // normalize
  lag_weights = lag_weights_un / sum(lag_weights_un);

  mu = a +
       b_fwi * (fwi_mat * lag_weights) +
       b_vfi * vfi + b_tfi * tfi +
       b_drz * drz + b_dhz * dhz;
}

model {
  // priors
  a_raw ~ std_normal();

  b_fwi_raw ~ std_normal();
  b_vfi_raw ~ std_normal();
  b_tfi_raw ~ std_normal();
  b_drz_raw ~ std_normal();
  b_dhz_raw ~ std_normal();

  sigma_raw ~ std_normal();
  g_raw ~ std_normal();
  ls_raw ~ std_normal();

  // likelihood
  y[ids_notna] ~ skew_normal(mu[ids_notna], sigma, g);
  target += skew_normal_lcdf(cutoff | mu[ids_na], sigma, g);
  // the second is the probability of being below log(0.09)
}