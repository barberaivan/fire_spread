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
  real priorscale_alpha;
  real priorscale_b_ndvi_log;
  real priormean_b_ndvi_log;
  real priorscale_b;
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
  vector[V] alpha_raw;
  real b_ndvi_log_raw;
  real<lower=0, upper=1> pi_ndvi_grass;
  real<lower=0, upper=1> pi_ndvi_shrub;
  vector<lower=0, upper=1>[V] optim_ndvi;
  real b_elev_raw;
  real b_north_raw;
  real b_slope_raw;
}

transformed parameters {
  vector[V] alpha = alpha_raw * priorscale_alpha;
  real b_ndvi_for = exp(b_ndvi_log_raw * priorscale_b_ndvi_log +
                        priormean_b_ndvi_log) * (-1);
  real b_ndvi_shrub = b_ndvi_for * pi_ndvi_shrub;
  real b_ndvi_grass = b_ndvi_for * pi_ndvi_grass;
  vector[V] b_ndvi;
  real b_elev = b_elev_raw * priorscale_b;
  real b_slope = b_slope_raw * priorscale_b;
  real b_north = b_north_raw * priorscale_b / north_scale;

  b_ndvi[1] = b_ndvi_for;
  b_ndvi[2] = b_ndvi_shrub;
  b_ndvi[3] = b_ndvi_grass;
}

model {
  // priors
  alpha_raw ~ std_normal();
  b_ndvi_log_raw ~ std_normal();
  b_elev_raw ~ std_normal();
  b_north_raw ~ std_normal();
  b_slope_raw ~ std_normal();

  // likelihood
  {
    vector[N] eta;

    for(n in 1:N) {
      eta[n] = alpha[veg_num[n]] +
               b_ndvi[veg_num[n]] * (ndvi[n] - optim_ndvi[veg_num[n]]) ^ 2 +
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
  vector[V] intercepts = alpha -
                         (elev_mean * b_elev_ori +
                          slope_mean * b_slope_ori +
                          north_mean * b_north_ori);
                          // b must be scaled to the original
}