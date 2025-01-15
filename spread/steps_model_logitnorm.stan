data {
  int<lower=0> N1;
  int<lower=0> N2;

  vector[N1] steps1;
  
  vector[N1] area1; // log area
  vector[N2] area2;
  real areaL;
  
  vector[N1] fwi1; 
  vector[N2] fwi2;
  
  real L; // lower steps = 2
  
  // priors
  real Umax;
  real Umin;
  
  real prior_steps_int_mn;
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
}

parameters {
  real steps_int;
  real<lower=0> steps_b;
  real<lower=Umin, upper=Umax> U;
  real<lower=0> steps_sigma;
  
  real area_int;
  real<lower=0> area_b;
  real<lower=0> area_sigma;
  
  vector[N2] steps_logit2;
}

transformed parameters {
  vector[N1] steps_logit1 = logit((steps1 - L) / (U - L));
  vector[N2] steps2 = inv_logit(steps_logit2) * (U - L) + L;
  vector[N2] steps_log2 = log(steps2);
}

model {
  // priors
  steps_int ~ normal(prior_steps_int_mn, prior_steps_int_sd);
  steps_b ~ normal(0, prior_steps_b_sd);
  steps_sigma ~ normal(0, prior_steps_sigma_sd);
  
  area_int ~ normal(prior_area_int_mn, prior_area_int_sd);
  area_b ~ normal(0, prior_area_b_sd);
  area_sigma ~ normal(0, prior_area_sigma_sd);
  
  // likelihoods
  steps_logit1 ~ normal(steps_int + steps_b * fwi1, steps_sigma);
  steps_logit2 ~ normal(steps_int + steps_b * fwi2, steps_sigma);
  
  for(n in 1:N1) {
    area1[n] ~ normal(area_int + area_b * steps_log1[n], area_sigma)T[areaL, ];
  }
  
  for(n in 1:N2) {
    area2[n] ~ normal(area_int + area_b * steps_log2[n], area_sigma)T[areaL, ];
  }
  
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
      1 / (steps1[n] - L) + 1 / (U - steps1[n])
    )); 
  }
}