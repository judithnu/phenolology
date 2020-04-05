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
  int<lower=1> G = 10; // Number of members in cluster
  real beta = 0.5; // population slope effect
  real beta_g_mu = 0.1; //group slope mean effect
  real beta_g_sd = 0.1;
}

generated quantities {
  //real beta = exp(3); // covariate effect
  real x[N]; //simulate covariate
  vector[N] gamma; //simulated latent effect
  vector[G] beta_g; //group effects
  vector[N] beta_g_dat;

  ordered[K - 1] c = induced_dirichlet_rng(K, 7); // (Internal) cut points
  int<lower=1, upper=K> y[N];                     // Simulated ordinals

  // build group effects
  for (g in 1:G) {
    beta_g[g] = normal_rng(beta_g_mu, beta_g_sd);
  }

  //beta_g_dat = rep_array(beta_g, 50);

  for (n in 1:N) {
    x[n] = uniform_rng(0,20); // covariate
    gamma[n] = (x[n]+group) * beta; // Latent effect
    y[n] = ordered_logistic_rng(gamma[n], c); // ordinals
  }
}

