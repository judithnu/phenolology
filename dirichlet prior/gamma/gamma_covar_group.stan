
data {
    int<lower=1> N;             // Number of observations
    int<lower=2> K;             // Number of ordinal categories
    int<lower=2> G;             // Number of groups

    int<lower=1, upper=K> y[N]; // Observed ordinals
    vector[N] x;                  // Covariate
    int GID[N]; //Groups

    real<lower=1> shape; // shape parameter for gamma prior on cutpoints
    real<lower=0> beta_rate; // rate parameter for exponential prior on beta
    real<lower=0> cut_rate; //
}

parameters {
    positive_ordered[K - 1] c; // (Internal) cut points
    real<lower=0> beta; // population level effect

    vector[G] betag;
    //real<lower=0> sigma_group;
}


model {
    vector[N] phi;

    // Prior model
    beta ~ exponential(3);
    c ~ gamma(10,1);
    betag ~ normal(0, 0.15);

    // Observational model

  for (i in 1:N ) {
    phi[i] = (beta * betag[GID[i]]) * x[i];
    y[i] ~ ordered_logistic(phi[i], c);
    }

}

generated quantities {
  row_vector[G] h1;
  row_vector[G] h2;

  for (g in 1:G) {
    h1[g] = c[1]/(beta + betag[g]);
    h2[g] = c[2]/(beta + betag[g]);
  }
}
