data {
  int<lower=0> J;
  int<lower=0> K;
  matrix[K, J] mus;
  matrix[K, K] Ps[J];
}

parameters {
  matrix[K, J] x;
}

model {
  for(j in 1:J)
    x[, j] ~ multi_normal_prec(mus[, j], Ps[j]);
}