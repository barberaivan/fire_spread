data {
  // ignition rate model
  int<lower=0> nfort;
  int<lower=0> nfort_both;
  int<lower=0> npix;
  int<lower=0> nlag;

  int ig_counts[nfort, 3];   // 1: human; 2: lightning; 3: unknown
  int condition[nfort];      // 1: known, either huamn or lightning;
                             // 2: some unknown; 3: only lightning

  matrix[nfort, nlag] fwi[npix];
  vector[npix] pix_weights;      // relative size of each pixel in the study area

  real prior_a_sd;
  real prior_a_mu;
  real prior_b_sd;
  real prior_phi_sd;
  real prior_ls_sd;

  // ignition location model
  int npoint;
  int nland;
  vector[nland] xland;
  vector[npoint] xpoints;
  int cause[npoint];
  int row_rate[npoint];

  real prior_c_sd;
}

parameters {
  vector[2] a_raw;
  vector<lower=0>[2] b_raw;
  vector<lower=0>[2] phi_raw;
  real<lower=0> ls_raw;

  vector[2] c_raw;
}

transformed parameters {
  // ignition rate model
  vector[2] a = a_raw * prior_a_sd + prior_a_mu;
  vector<lower=0>[2] b = b_raw * prior_b_sd;
  vector<lower=0>[2] phi = phi_raw * prior_phi_sd;
  real ls = ls_raw * prior_ls_sd;

  matrix[nfort, npix] fwi_mat;
  matrix[nfort, 2] mu;

  vector[nlag] lag_weights_improper;
  vector[nlag] lag_weights;

  // ignition location model
  vector[2] c = c_raw * prior_c_sd;
  vector[2] theta_land; // sums of unnormalized selection probabilities in the
                        // landscape

  // compute variables for both models

  for(t in 1:nlag)
    lag_weights_improper[t] = exp(-((t-1) / ls) ^ 2);
  // normalize weights
  lag_weights = lag_weights_improper / sum(lag_weights_improper);

  for(p in 1:npix) {
    fwi_mat[, p] = fwi[p] * lag_weights;
  }

  for(i in 1:2) mu[, i] = exp(a[i] + b[i] * fwi_mat) * pix_weights;
  // mu is the weighted mean of each pixel's lambda.

  for(i in 1:2) theta_land[i] = sum(exp(xland * c[i]));

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
     ig_counts[i, 1] ~ neg_binomial_2(mu[i, 1], phi[1]); // human
     ig_counts[i, 2] ~ neg_binomial_2(mu[i, 2], phi[2]); // lightning
    }

    // some unknown
    if(condition[i] == 2) {
     int nunk = ig_counts[i, 3];
     vector[nunk+1] ll;
     for(k in 0:nunk) {
       int n_human = ig_counts[i, 1] + k;
       int n_light = ig_counts[i, 2] + (nunk - k);
       ll[k+1] = neg_binomial_2_lpmf(n_human | mu[i, 1], phi[1]) +
                 neg_binomial_2_lpmf(n_light | mu[i, 2], phi[2]);
     }
     target += log_sum_exp(ll); // likelihood marginal to cause
    }

    // known, only lightning
    if(condition[i] == 3) {
     ig_counts[i, 2] ~ neg_binomial_2(mu[i, 2], phi[2]); // lightning
    }
  }

  // Ignition location model ---------------------------------------------

  // priors
  c_raw ~ std_normal();

  // likelihood
  for(i in 1:npoint) {
    // human- or lightning-caused (cause \in 1:2)
    if(cause[i] < 3) {
     int cc = cause[i];
     real pp = exp(xpoints[i] * c[cc]);
     target += log(pp / (pp + theta_land[cc]));
    }

    // unknown cause
    if(cause[i] == 3) {
      vector[2] like;
      int rr = row_rate[i];
      row_vector[2] props = mu[rr, ] / sum(mu[rr, ]);
      for(k in 1:2) {
        real pp = exp(xpoints[i] * c[k]);
        like[k] = (pp / (pp + theta_land[k])) * props[k];
      }
      target += log(sum(like)); // likelihood marginal to cause
    }
  }
}

generated quantities {
  real muh = sum(mu[1:nfort_both, 1]);
  real mul = sum(mu[1:nfort_both, 2]);
  real humanp = muh / (muh + mul);
}