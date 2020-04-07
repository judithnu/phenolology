// Simulate ordinal logistic data with a covariate and group effect

data {
  int<lower=1> N;
  int<lower=2> K;
  int<lower=1> G;

  vector[N] x; //covariate
  positive_ordered[K-1] c; //cutpoints
  vector[N] gbeta; //group effects
  real beta; //slope
}

generated quantities {
  int<lower=1, upper=K> y[N];                     // Simulated ordinals
  vector[N] gamma;

  for (n in 1:N) {
    gamma[n] = (x[n] + gbeta[n]) * beta; // Latent effect
    y[n] = ordered_logistic_rng(gamma[n], c); // ordinals
  }
}

