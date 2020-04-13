
data {
    int<lower=1> N;             // Number of observations
    int<lower=2> K;             // Number of ordinal categories
    int<lower=2> G;             // Number of groups

    int<lower=1, upper=K> y[N]; // Observed ordinals
    vector[N] x;                  // Covariate

    int GID[N]; //Groups
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
    //betag ~ normal(0,0.25);
    //sigma_group ~ exponential(3);
    betag ~ normal(0, 0.25);


    // Observational model

  for (i in 1:N ) {
    phi[i] = (beta * betag[GID[i]]) * x[i];
    y[i] ~ ordered_logistic(phi[i], c);
    }

}

generated quantities {
  row_vector[G] fstart;

  for (g in 1:G) {
    fstart[g] = c[1]/(beta + betag[g]);
  }
}
