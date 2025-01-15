data {
  // ignition rate model
  int nfort;        // fortnights in the whole study period
  int nlag;         // lag length for FWI
  int nig;          // ignition types (human, lightning)

  int ig_counts[nfort, nig];   // 1: human; 2: lightning; 3: unknown
  
  matrix[nfort, nlag] fwi_mat;
  real study_area_size;         // km2, to have lambda = rate / (km2 * fort)
                                 
  real prior_a_sd;
  real prior_a_mean;
  real prior_b_sd;
  real prior_phi_sd;
  real prior_ls_sd;
  real prior_U_sd;
  real prior_U_mean;
}

parameters {
  vector[nig] a_raw;
  vector<lower=0>[nig] b_raw;
  vector<lower=0>[nig] U_raw;
  vector<lower=0>[nig] phi_raw;
  real<lower=0> ls_raw;
}

transformed parameters {
  // ignition rate model
  vector[nig] a = a_raw * prior_a_sd + prior_a_mean;   // intercepts
  vector[nig] b = b_raw * prior_b_sd;                  // FWI slopes
  vector[nig] U = U_raw * prior_U_sd + prior_U_mean;   // igrate upper limits
  vector[nig] phi = phi_raw * prior_phi_sd;            // dispersion
  real ls = ls_raw * prior_ls_sd;                      // FWI temporal corr

  vector[nfort] fwi; // cumulative, weighted
  matrix[nfort, nig] lambda;

  vector[nlag] lag_weights_un; // unnormalized
  vector[nlag] lag_weights;    // time-weights for FWI

  // cumulative FWI
  for(t in 1:nlag)
    lag_weights_un[t] = exp(-((t-1) / ls) ^ 2); // t-1 to start at zero
  // normalize
  lag_weights = lag_weights_un / sum(lag_weights_un);
  fwi = fwi_mat * lag_weights;
  
  // lambda for ignition rates
  for(i in 1:nig) {
    lambda[, i] = inv_logit(a[i] + b[i] * fwi) * U[i] * study_area_size;
  }
}

model {
  // Ignition rate model ------------------------------------------------

  // priors
  a_raw ~ std_normal();
  b_raw ~ std_normal();
  phi_raw ~ std_normal();
  ls_raw ~ std_normal();
  U_raw ~ std_normal();

  // likelihood
  ig_counts[, 1] ~ neg_binomial_2(lambda[, 1], phi[1]); // human
  ig_counts[, 2] ~ neg_binomial_2(lambda[, 2], phi[2]); // lightning
}