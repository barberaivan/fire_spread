data {
  int<lower = 0> N;
  int<lower = 0> V; // vegetation types

  int y[N];
  int veg_num[N];   // veg type {1: forest, 2: shrubland, 3: grassland}
  vector[N] ndvi;   // unitless, in [-1, 1]
  vector[N] elev;   // m
  vector[N] north;  // unitless, in [-1, 1]
  vector[N] slope;  // degrees [0, 90]

  // prior parameters
  // For NDVI hyper-parameters
  real priorscale_a_mu;     // 10
  real priorscale_a_sigma;   // 5

  // b at log scale
  real priormean_b_mu;      // log(100)
  real priorscale_b_mu;     // 2
  real priorscale_b_sigma;   // 2

  // optim at logit scale
  real priorscale_o_mu;     // 5
  real priorscale_o_sigma;   // 3

  // prior scale for the slope of topographic variables
  real priorscale_b_topo;   // 10
}

transformed data {
  // Scaling is done to sample parameters in a convenient scale for HMC, but
  // then parameters are scaled so the raw data can be used to compute the
  // linear predictor.

  // compute means to center predictors or quadratic terms
  real elev_mean = mean(elev);
  real slope_mean = mean(slope);
  real north_mean = mean(north);

  real elev_scale = sd(elev);
  real north_scale = sd(north);
  real slope_scale = sd(slope);

  // standardize topographic variables
  vector[N] elev_z = (elev - elev_mean) / elev_scale;
  vector[N] slope_z = (slope - slope_mean) / slope_scale;
  vector[N] north_z = (north - north_mean) / north_scale;
}

parameters {
  // ndvi parameters for each vegetation type (random effects)
  vector[V] a_raw;
  vector[V] b_log_raw;
  vector[V] o_logit_raw;

  // hyper parameters, at the unconstrained scale
  real a_mu_raw;
  real<lower=0> a_sigma_raw;

  real b_mu_raw;
  real<lower=0> b_sigma_raw;

  real o_mu_raw;
  real<lower=0> o_sigma_raw;

  // topographic parameters
  real b_elev_raw;
  real b_north_raw;
  real b_slope_raw;
}

transformed parameters {
  // ndvi hyper parameters
  real a_mu = a_mu_raw * priorscale_a_mu;
  real a_sigma = a_sigma_raw * priorscale_a_sigma;

  real b_mu = b_mu_raw * priorscale_b_mu;
  real b_sigma = b_sigma_raw * priorscale_b_sigma;

  real o_mu = o_mu_raw * priorscale_o_mu;
  real o_sigma = o_sigma_raw * priorscale_o_sigma;

  // ndvi parameters
  vector[V] a = a_raw * a_sigma + a_mu;
  vector[V] b = exp(b_log_raw * b_sigma + b_mu) * (-1);
  vector[V] o = inv_logit(o_logit_raw * o_sigma + o_mu);

  // topographic parameters
  real b_elev = b_elev_raw * priorscale_b_topo;
  real b_slope = b_slope_raw * priorscale_b_topo;
  real b_north = b_north_raw * priorscale_b_topo / north_scale;
}

model {
  // priors
  a_mu_raw ~ std_normal();
  a_sigma_raw ~ std_normal();
  a_raw ~ std_normal();

  b_mu_raw ~ std_normal();
  b_sigma_raw ~ std_normal();
  b_log_raw ~ std_normal();

  o_mu_raw ~ std_normal();
  o_sigma_raw ~ std_normal();
  o_logit_raw ~ std_normal();

  b_elev_raw ~ std_normal();
  b_north_raw ~ std_normal();
  b_slope_raw ~ std_normal();

  // likelihood
  {
    vector[N] eta;

    for(n in 1:N) {
      eta[n] = a[veg_num[n]] +
               b[veg_num[n]] * (ndvi[n] - o[veg_num[n]]) ^ 2 +
               elev_z[n] * b_elev +
               north_z[n] * b_north +
               slope_z[n] * b_slope;
    }

    y ~ bernoulli_logit(eta);
  }
}

generated quantities {
  // linear slopes at original scale
  real b_elev_ori = b_elev / elev_scale;
  real b_slope_ori = b_slope / slope_scale;
  real b_north_ori = b_north / north_scale;

  // non-centered intercepts (quadratic terms are not centered)
  vector[V] intercepts = a -
                         (elev_mean * b_elev_ori +
                          slope_mean * b_slope_ori +
                          north_mean * b_north_ori);
                          // b must be scaled to the original
}