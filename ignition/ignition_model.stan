data {
  // ignition rate model
  int nfort;        // fortnights in the whole study period
  int nfort_both;   // fortnights for both ignition causes
  int npix;         // FWI pixels in the study area
  int nlag;         // lag length for FWI
  int nig;          // ignition types (human, lightning)


  int ig_counts[nfort, nig+1];   // 1: human; 2: lightning; 3: unknown
  int condition[nfort];          // 1: known, either huamn or lightning;
                                 // 2: some unknown; 3: only lightning

  matrix[nfort, nlag] fwi[npix];
  vector[npix] pix_weights;      // relative size of each pixel in the study area

  real prior_a_sd;
  real prior_a_mean;
  real prior_b_sd;
  real prior_phi_sd;
  real prior_ls_sd;

  // ignition location model
  int npoint; // ignition points
  int nland;  // sampled pixels in the whole study area
  int ny;     // years, to use the correct NDVI to compute the VFI
  int nfi;    // flammability indices (TFI, VFI)
  int ndist;  // distance variables (human settlements and roads)

  matrix[npoint, nfi] X_ig_fi;     // design matrices
  // vector[npoint] ig_fwi;           // long-term mean of summer-fwi (spatial)
  matrix[npoint, ndist] X_ig_dist;
  matrix[nland, nfi] X_pop_fi[ny]; // the VFI = f(NDVI) varies by year
  // vector[nland] pop_fwi;
  matrix[nland, ndist] X_pop_dist;
  // ig for ignited points, pop for population of pixels sample.

  int cause[npoint];    // cause of each ignition (human, light, unk)
  int fort_id[npoint];  // fortnight when it occurred
  int year_id[npoint];  // year when it occurred

  real prior_c_sd;
}

parameters {
  vector[nig] a_raw;
  vector<lower=0>[nig] b_raw;
  vector<lower=0>[nig] phi_raw;
  real<lower=0> ls_raw;

  matrix<lower=0>[nfi, nig] c_fi_raw;  // [{tfi, vfi}, {human, lightning}]
  // vector<lower=0>[nig] c_fwi_raw;      // avg fwi effect
  vector<upper=0>[ndist] c_dist_raw;   // distance effect for human ignitions
}

transformed parameters {
  // ignition rate model
  vector[nig] a = a_raw * prior_a_sd + prior_a_mean; // intercepts
  vector[nig] b = b_raw * prior_b_sd;                  // FWI slopes
  vector[nig] phi = phi_raw * prior_phi_sd;            // dispersion
  real ls = ls_raw * prior_ls_sd;                      // FWI temporal corr

  matrix[nfort, npix] fwi_mat;
  matrix[nfort, nig] lambda;

  vector[nlag] lag_weights_un; // unnormalized
  vector[nlag] lag_weights;    // time-weights for FWI

  // ignition location model
  matrix[nfi, nig] c_fi = c_fi_raw * prior_c_sd;  // VFI-TFI effects
  // vector[nig] c_fwi = c_fwi_raw * prior_c_sd;     // avg fwi effect
  vector[ndist] c_dist = c_dist_raw * prior_c_sd; // distances effects
  matrix[ny, nig] theta_land; // sums of unnormalized selection probabilities
                              // in the landscape by year and ignition type
  matrix[npoint, nig] theta_point; // unnormalized selection probabilities
                                   // by pixel and ignition type

  matrix[npoint, nig] cause_prob_all;
  // compute variables for both models

  for(t in 1:nlag)
    lag_weights_un[t] = exp(-((t-1) / ls) ^ 2); // t-1 to start at zero
  // normalize
  lag_weights = lag_weights_un / sum(lag_weights_un);

  for(p in 1:npix) {
    fwi_mat[, p] = fwi[p] * lag_weights;
  }

  for(i in 1:nig) lambda[, i] = exp(a[i] + b[i] * fwi_mat) * pix_weights;
  // lambda is the weighted mean of each pixel's lambda.

  // compute theta land
  {
    vector[nland] eta_dist = X_pop_dist * c_dist;
    vector[nland] eta_fi;

    for(y in 1:ny) {
      eta_fi = X_pop_fi[y] * c_fi[, 1];  // human
      theta_land[y, 1] = sum(exp(eta_fi + eta_dist));

      eta_fi = X_pop_fi[y] * c_fi[, 2];  // lightning
      theta_land[y, 2] = sum(exp(eta_fi));
    }
  }

  // The X for ignition points (X_ig_fi) already includes year-varying VFI.
  theta_point[, 1] = exp(X_ig_fi * c_fi[, 1] +
                         // ig_fwi * c_fwi[1] +
                         X_ig_dist * c_dist);
  theta_point[, 2] = exp(X_ig_fi * c_fi[, 2]);// +
                         // ig_fwi * c_fwi[2]);
  // human cause considers distance to settlements and roads; lightning
  // does not
}

model {
  // Ignition rate model ------------------------------------------------

  // priors
  a_raw ~ std_normal();
  b_raw ~ std_normal();
  phi_raw ~ std_normal();
  ls_raw ~ std_normal();

  // likelihood
  for(i in 1:nfort) {
    // known, human or lightning
    if(condition[i] == 1) {
     ig_counts[i, 1] ~ neg_binomial_2(lambda[i, 1], phi[1]); // human
     ig_counts[i, 2] ~ neg_binomial_2(lambda[i, 2], phi[2]); // lightning
    }

    // some unknown
    if(condition[i] == 2) {
     int nunk = ig_counts[i, 3];
     vector[nunk+1] ll;
     for(k in 0:nunk) {
       int n_human = ig_counts[i, 1] + k;
       int n_light = ig_counts[i, 2] + (nunk - k);
       ll[k+1] = neg_binomial_2_lpmf(n_human | lambda[i, 1], phi[1]) +
                 neg_binomial_2_lpmf(n_light | lambda[i, 2], phi[2]);
     }
     target += log_sum_exp(ll); // likelihood marginal to cause
    }

    // known, only lightning
    if(condition[i] == 3) {
     ig_counts[i, 2] ~ neg_binomial_2(lambda[i, 2], phi[2]); // lightning
    }
  }

  // Ignition location model ---------------------------------------------

  // priors
  to_vector(c_fi_raw) ~ std_normal();
  // c_fwi_raw ~ std_normal();
  c_dist_raw ~ std_normal();

  // likelihood
  for(i in 1:npoint) {
    // human- or lightning-caused (cause \in 1:2)
    if(cause[i] < 3) {
      int cc = cause[i];
      int yy = year_id[i];
      real pp = theta_point[i, cc];
      target += log(pp / (pp + theta_land[yy, cc]));
    }

    // unknown cause
    if(cause[i] == 3) {
      vector[2] like;
      int ff = fort_id[i];
      int yy = year_id[i];
      row_vector[2] props = lambda[ff, ] / sum(lambda[ff, ]);
      for(k in 1:2) {
        real pp = theta_point[i, k];
        like[k] = (pp / (pp + theta_land[yy, k])) * props[k];
      }
      target += log(sum(like)); // likelihood marginal to cause
    }
  }
}

generated quantities {
  real lambda_h = sum(lambda[1:nfort_both, 1]);
  real lambda_l = sum(lambda[1:nfort_both, 2]);
  real human_prop = lambda_h / (lambda_h + lambda_l);

  real c_vfi_h = c_fi[1, 1]; // VFI human
  real c_vfi_l = c_fi[1, 2]; // VFI lightning

  real c_tfi_h = c_fi[2, 1]; // TFI human
  real c_tfi_l = c_fi[2, 2]; // TFI lightning

  // real c_fwi_h = c_fwi[1]; // FWI human
  // real c_fwi_l = c_fwi[2]; // FWI lightning

  real c_dist_h = c_dist[1]; // dist human settlements
  real c_dist_r = c_dist[2]; // dist roads
}