data {
  int<lower=0> J;
  int<lower=0> K;
  int<lower=0> Kf;
  int<lower=0> Kr;
  int<lower=0> idfix[Kf];
  int<lower=0> idran[Kr];
  matrix[K, J] mus;
  matrix[K, K] Ps[J];
}

parameters {
  // varying coefficients
  matrix[Kr, J] x;
  // constant coefficients
  vector[Kf] xf;
}

model {
  vector[K] xmerge; 
  for(j in 1:J) {
    xmerge[idfix] = xf;
    xmerge[idran] = x[, j];
    xmerge ~ multi_normal_prec(mus[, j], Ps[j]);
  }
}