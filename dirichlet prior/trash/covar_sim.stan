// Simulate from a model with a gamma distribution

data {
  int<lower=1> N;
  int<lower=2> K;

  vector[N] x; //covariate
  positive_ordered[K-1] c; //cutpoints
}

transformed data {
  real beta = 0.5; // population slope effect
}

generated quantities {
  int<lower=1, upper=K> y[N];                     // Simulated ordinals
  vector[N] gamma;

  for (n in 1:N) {
    gamma[n] = x[n] * beta; // Latent effect
    y[n] = ordered_logistic_rng(gamma[n], c); // ordinals
  }
}

