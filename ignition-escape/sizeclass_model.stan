data {
  int n;        // number of ignitions
  int nlag;     // lag length for FWI
  int K;        // size classes

  int y[n];

  matrix[n, nlag] fwi_mat_z;
  vector[n] vfi;
  vector[n] tfi;
  vector[n] drz;
  vector[n] dhz;

  real prior_a_sd;
  real prior_b_sd;
  real prior_ls_sd;
}

transformed data {
  vector[n] ones = rep_vector(1.0, n);
}

parameters {
  ordered[K-1] a_raw;

  real<lower=0> b_fwi_raw;
  real<lower=0> b_vfi_raw;
  real<lower=0> b_tfi_raw;
  real b_drz_raw;
  real b_dhz_raw;

  real<lower=0> ls_raw;
}

transformed parameters {
  ordered[K-1] a = a_raw * prior_a_sd;

  real b_fwi = b_fwi_raw * prior_b_sd;
  real b_vfi = b_vfi_raw * prior_b_sd;
  real b_tfi = b_tfi_raw * prior_b_sd;
  real b_drz = b_drz_raw * prior_b_sd;
  real b_dhz = b_dhz_raw * prior_b_sd;

  real ls = ls_raw * prior_ls_sd;                      // FWI temporal corr

  vector[nlag] lag_weights_un; // unnormalized
  vector[nlag] lag_weights;    // time-weights for FWI

  vector[n] eta;
  matrix[n, K] cmf;
  matrix[n, K] pmf;

  for(t in 1:nlag)
    lag_weights_un[t] = exp(-((t-1) / ls) ^ 2); // t-1 to start at zero
  // normalize
  lag_weights = lag_weights_un / sum(lag_weights_un);

  eta = b_fwi * (fwi_mat_z * lag_weights) +
        b_vfi * vfi + b_tfi * tfi +
        b_drz * drz + b_dhz * dhz;

  // compute cumulative and probability mass functions
  for(k in 1:(K-1)) {
    cmf[, k] = inv_logit(a[k] - eta);
  }
  cmf[, K] = ones;

  pmf[, 1] = cmf[, 1];
  for(k in 2:K) {
    pmf[, k] = cmf[, k] - cmf[, k-1];
  }
}

model {
  // priors
  a_raw ~ std_normal();

  b_fwi_raw ~ std_normal();
  b_vfi_raw ~ std_normal();
  b_tfi_raw ~ std_normal();
  b_drz_raw ~ std_normal();
  b_dhz_raw ~ std_normal();

  ls_raw ~ std_normal();

  // likelihood
  for(i in 1:n) {
    y[i] ~ categorical(pmf[i, ]');
  }
}