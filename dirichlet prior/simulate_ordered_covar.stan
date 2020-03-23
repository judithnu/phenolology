// Modified from Betancourt's case study https://betanalpha.github.io/assets/case_studies/ordinal_regression.html to add a covariate

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
  int<lower=1> N = 500; // Number of observations
  int<lower=1> K = 3;  // Number of ordinal categories
  real beta = 0.5;
}

generated quantities {
  //real beta = exp(3); // covariate effect
  real x[N]; //simulate covariate
  vector[N] gamma; //simulated latent effect

  ordered[K - 1] c = induced_dirichlet_rng(K, 7); // (Internal) cut points
  int<lower=1, upper=K> y[N];                     // Simulated ordinals


  for (n in 1:N) {
    x[n] = uniform_rng(0,20); // covariate
    gamma[n] = x[n] * beta; // Latent effect
    y[n] = ordered_logistic_rng(gamma[n], c); // ordinals
  }
}

