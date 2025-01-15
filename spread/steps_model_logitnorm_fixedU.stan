data {
  int<lower=0> N1;
  int<lower=0> N2;
  
  vector[N1] steps1;
  
  vector[N1] area1; // log area
  vector[N2] area2;
  
  vector[N1] fwi1; 
  vector[N2] fwi2;
  
  real L; // lower steps = 2
  real U; // upper steps = 2000
  
  // priors
  real prior_steps_int_sd;
  real prior_steps_b_sd;
  real prior_steps_sigma_sd;
  
  real prior_area_int_mn;
  real prior_area_int_sd;
  real prior_area_b_sd;
  real prior_area_sigma_sd;
}

transformed data {
  vector[N1] steps_log1 = log(steps1);
  vector[N1] steps_logit1 = logit((steps1 - L) / (U - L));
}

parameters {
  real steps_int;
  real<lower=0> steps_b;
  real<lower=0> steps_sigma;
  
  real area_int;
  real<lower=0> area_b;
  real<lower=0> area_sigma;
  
  vector[N2] steps_logit2;
}

transformed parameters {
  vector[N2] steps2 = inv_logit(steps_logit2) * (U - L) + L;
  vector[N2] steps_log2 = log(steps2);
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // priors
  steps_int ~ normal(0, prior_steps_int_sd);
  steps_b ~ normal(0, prior_steps_b_sd);
  steps_sigma ~ normal(0, prior_steps_sigma_sd);
  
  area_int ~ normal(prior_area_int_mn, prior_area_int_sd);
  area_b ~ normal(0, prior_area_b_sd);
  area_sigma ~ normal(0, prior_area_sigma_sd);
  
  // likelihoods
  steps_logit1 ~ normal(steps_int + steps_b * fwi1, steps_sigma);
  steps_logit2 ~ normal(steps_int + steps_b * fwi2, steps_sigma);
  
  area1 ~ normal(area_int + area_b * steps_log1, area_sigma);
  area2 ~ normal(area_int + area_b * steps_log2, area_sigma);
}