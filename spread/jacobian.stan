functions {
  real exp_logit(real x, real L, real U) {
    real exp_x = exp(x);
    real exp_unit = (exp_x - L) / (U - L);
    return logit(exp_unit);
  }

  // derivative of exp_logit
  real exp_logit_d(real x, real L, real U) {
    real exp_x = exp(x);
    return (exp_x * (L - U)) / ((exp_x - L) * (exp_x - U));
  }
}

data {
 real L;
 real U;
 real mu;
 real sigma;
}

parameters {
  real steps_log;
}

transformed parameters {
  real steps = exp(steps_log);
  real steps_unit = (steps - L) / (U - L);
  real steps_logit = logit(steps_unit);
  real steps_logit2 = exp_logit(steps_log, L, U);
}

model {
  steps_logit ~ normal(mu, sigma);
  // jacobian adjustment:
  target += log(fabs(exp_logit_d(steps_log, L, U)));
}