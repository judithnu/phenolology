//simulate from an ordinal distribution

functions {
  vector induced_dirichlet_rng(int K, real phi) {
    vector[K - 1] c;
    vector[K] p = dirichlet_rng(rep_vector(1, K));
    
    c[1] = phi - logit(1 - p[1]);
    for (k in 2:(K - 1))
      c[k] = phi - logit(inv_logit(phi - c[k - 1]) - p[k]);
    
    return c;
  }
}

transformed data {
  int<lower=1> N = 50; // Number of observations
  int<lower=1> K = 3;  // Number of ordinal categories
  vector[N] forcing; //forcing
}

generated quantities {
  real beta = exponential_rng(2); // population effect
  real gamma = forcing * beta;                  // Latent effect
  ordered[K - 1] c = induced_dirichlet_rng(K, gamma); // (Internal) cut points
  int<lower=1, upper=K> y[N];                     // Simulated ordinals
  
  for (n in 1:N)
    y[n] = ordered_logistic_rng(gamma, c);
}

