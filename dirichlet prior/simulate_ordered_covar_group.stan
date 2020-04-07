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

data {
  int<lower=1> N;
  int<lower=1> G;
  int<lower=2> K;

  vector[N] x; //covariate
  vector[N] gbeta; //group effect
}

transformed data {
  real beta = 0.5; // population slope effect
}

generated quantities {
  int<lower=1, upper=K> y[N];                     // Simulated ordinals
  vector[N] gamma;

  ordered[K - 1] c = induced_dirichlet_rng(K, 7); // (Internal) cut points

  for (n in 1:N) {
    gamma[n] = (x[n] + gbeta[n]) * beta; // Latent effect
    y[n] = ordered_logistic_rng(gamma[n], c); // ordinals
  }
}

