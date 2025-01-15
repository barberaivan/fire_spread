data {
  int<lower=0> N1;
  int<lower=0> N2;
  
  vector[N1] steps1;
  
  vector[N1] area1; // log area
  vector[N2] area2;
  
  vector[N1] fwi1; 
  vector[N2] fwi2;
  
  real steps_lower;
  
  // priors
  real upper_max;
  real upper_min;
  
  real prior_steps_mid_mn;
  real prior_steps_mid_sd;
  real prior_steps_b_sd;
  real prior_steps_sigma_sd;
  
  real prior_area_int_mn;
  real prior_area_int_sd;
  real prior_area_b_sd;
  real prior_area_sigma_sd;
}

transformed data {
  vector[N1] steps_log1 = log(steps1);
}

parameters {
  real steps_mid;
  real<lower=0> steps_b;
  real<lower=upper_min, upper=upper_max> steps_upper;
  real<lower=0> steps_sigma;
  
  real area_int;
  real<lower=0> area_b;
  real<lower=0> area_sigma;
  
  vector[N2] steps_logit2;
}

transformed parameters {
  real range = steps_upper - steps_lower; // steps_lower = 2
  vector[N1] steps_logit1 = logit((steps1 - steps_lower) / range);
  vector[N2] steps2 = inv_logit(steps_logit2) * range + steps_lower;
  vector[N2] steps_log2 = log(steps2);
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // priors
  steps_mid ~ normal(prior_steps_mid_mn, prior_steps_mid_sd);
  steps_b ~ normal(0, prior_steps_b_sd);
  steps_sigma ~ normal(0, prior_steps_sigma_sd);
  
  area_int ~ normal(prior_area_int_mn, prior_area_int_sd);
  area_b ~ normal(0, prior_area_b_sd);
  area_sigma ~ normal(0, prior_area_sigma_sd);
  
  // likelihoods
  steps_logit1 ~ normal(steps_b * (fwi1 - steps_mid), steps_sigma);
  steps_logit2 ~ normal(steps_b * (fwi2 - steps_mid), steps_sigma);
  
  area1 ~ normal(area_int + area_b * steps_log1, area_sigma);
  area2 ~ normal(area_int + area_b * steps_log2, area_sigma);
  
  /*
    add the log-absolute jacobian for the transform from steps1 to steps_logit1.
    transform: y = logit((x - L) / (U - L))
    where x = steps1, a = steps_upper
  
    expand logit
    y = log( ((x - L) / (U - L)) / (1 - ((x - L) / (U - L)))  )
  
    according to wolfram alpha:
    alternative form:
    y = log(x - L) - log(U - x)
  
    derivative:
    dy / dx = 1 / (x - L) + 1 / (U - x)
  */
  
  for(n in 1:N1) {
    target += log(fabs(
      1 / (steps1[n] - steps_lower) + 1 / (steps_upper - steps1[n])
    )); 
  }
}



