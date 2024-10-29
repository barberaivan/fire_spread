data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  real a;
  real s0;
  real t0;
  real d0;
}

parameters {
  real b0;
  real b1;
  real<lower=0> sigma2;
}

transformed parameters {
  vector[N] mu = b0 + b1 * x;
  real sigma = sqrt(sigma2);
}

model {
  for(n in 1:N)
    y[n] ~ normal(mu[n], sigma)T[a, ];

  b0 ~ normal(0, s0);
  b1 ~ normal(0, s0);
  sigma2 ~ inv_gamma(t0, d0);
}