// Add pp (quadratic) and distances to humans

data {
  int<lower = 0> N;
  int<lower = 0> V; // vegetation types
  int y[N];

  matrix[N, V] veg; // binary matrix, columns: [wet, subalpine, dry, shrubland, steppe]

  vector[N] ndvi;   // unitless, in [-1, 1]
  vector[N] pp;     // mm / year

  vector[N] elev;   // m
  vector[N] tpi;    // unitless, in [0, 1]
  vector[N] north;  // unitless, in [-1, 1]
  vector[N] slope;  // degrees [0, 90]

  vector[N] dist_human; // scaled
  vector[N] dist_roads; // scaled

  // prior scales
  real prior_sd_linear; // 5
  real prior_sd_quad;   // 15
  real prior_sd_int;    // 5

  // ndvi means and scales by veg type (too verbose to compute in stan)
  vector[V] ndvi_mean;
  vector[V] ndvi_scale;
}

transformed data {
  // Compute centers and scales to transform parameters from the raw to the
  // data scale. For the linear effects, scale is just sd(x), but for
  // quadratic effects, it is sd((x - mean(x)) ^ 2).
  // Scaling is done to sample parameters in a convenient scale for HMC, but
  // then parameters are scaled so the raw data can be used to compute the
  // linear predictor.

  // compute means to center predictors or quadratic terms
  real pp_mean = mean(pp);
  real elev_mean = mean(elev);
  real tpi_mean = mean(tpi);
  real slope_mean = mean(slope);

  real pp_scale = sd((pp - pp_mean) .* (pp - pp_mean));
  real elev_scale = sd((elev - elev_mean) .* (elev - elev_mean));
  real tpi_scale = sd((tpi - tpi_mean) .* (tpi - tpi_mean));
  real north_scale = sd(north);
  real slope_scale = sd(slope);

  // center slope (north has mean near 0)
  vector[N] slope_z = slope - mean(slope);
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  // intercepts
  vector[V] alpha_raw;

  // quadratic slopes (b < 0) and optimums
  vector<upper = 0>[V] b_ndvi_raw;
  real<upper = 0> b_pp_raw;
  real<upper = 0> b_elev_raw;
  real<upper = 0> b_tpi_raw;

  real b_human_raw;
  real b_roads_raw;

  vector<lower = 0, upper = 1>[V] optim_ndvi;
  real<lower = 0, upper = 1> optim_pp_raw;   // [380, 2200]
  real<lower = 0, upper = 1> optim_elev_raw; // upper = 2500
  real<lower = 0, upper = 1> optim_tpi;

  // linear slopes
  real<lower = 0> b_north_raw;
  real<lower = 0> b_slope_raw;
}

transformed parameters {
  // intercepts
  vector[V] alpha = alpha_raw * prior_sd_int;
  vector[V] intercepts; // compute later

  // quadratic slopes (b < 0) and optimums
  vector[V] b_ndvi; // need a loop, compute at the end of the block
  real b_pp = b_pp_raw * prior_sd_quad / pp_scale;
  real b_elev = b_elev_raw * prior_sd_quad / elev_scale;
  real b_tpi = b_tpi_raw * prior_sd_quad / tpi_scale;

  real optim_pp = optim_pp_raw * (2200 - 380) + 380;
  real optim_elev = optim_elev_raw * 2500;

  // linear slopes
  real b_north = b_north_raw * prior_sd_linear / north_scale;
  real b_slope = b_slope_raw * prior_sd_linear / slope_scale;

  real b_human = b_human_raw * prior_sd_linear;
  real b_roads = b_roads_raw * prior_sd_linear;

  // ndvi slopes
  for(v in 1:V)
    b_ndvi[v] = b_ndvi_raw[v] * prior_sd_quad / ndvi_scale[v];

  // non-centered intercepts
  intercepts = alpha - (slope_mean * b_slope);
                        // north already has mean 0
                        // quadratic terms not centred.

                        // distances are always treated as if they were in the
                        // original scale, because we don't care about them.

}

model {

  // priors
  alpha_raw ~ std_normal();

  b_ndvi_raw ~ std_normal();
  b_pp_raw ~ std_normal();
  b_elev_raw ~ std_normal();
  b_tpi_raw ~ std_normal();

  b_human_raw ~ std_normal();
  b_roads_raw ~ std_normal();

  b_north_raw ~ std_normal();
  b_slope_raw ~ std_normal();

  // likelihood
  {
    // linear predictor
    vector[N] eta;

    // variables for quadratic terms
    matrix[N, V] ndvi_diff_quad;
    vector[N] ndvi_diff;
    vector[N] pp_diff;
    vector[N] elev_diff;
    vector[N] tpi_diff;

    // ndvi effects, varying by vegetation type
    for(v in 1:V) {
      ndvi_diff = ndvi - optim_ndvi[v];
      ndvi_diff_quad[, v] = (ndvi_diff .* ndvi_diff) .* veg[, v];
      // this centers the term and selects the vegetation type
    }

    // remaining differences for quadratic terms
    pp_diff = pp - optim_pp;
    elev_diff = elev - optim_elev;
    tpi_diff = tpi - optim_tpi;

    // add all terms
    eta = veg * alpha +             // intercepts for centered predictors/effects
          ndvi_diff_quad * b_ndvi +
          (pp_diff .* pp_diff) * b_pp +
          (elev_diff .* elev_diff) * b_elev +
          (tpi_diff .* tpi_diff) * b_tpi +
          north * b_north +         // naturally centered
          slope_z * b_slope +       // slope_z is centered
          dist_human * b_human +    // both distances are centered
          dist_roads * b_roads;

    y ~ bernoulli_logit(eta);
  }
}