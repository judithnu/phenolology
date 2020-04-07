
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
    real<lower=0> sigma_group;
}

model {
    vector[N] phi;

    // Prior model
    beta ~ exponential(3);
    c ~ gamma(10,1);
    //betag ~ normal(0,0.25);
    sigma_group ~ exponential(4);
    betag ~ normal(0, sigma_group);


    // Observational model

  for (i in 1:N ) {
    phi[i] = (beta * betag[GID[i]]) * x[i];
    y[i] ~ ordered_logistic(phi[i], c);
    }

}

// generated quantities {
    //add fstart and fend here
//     vector[N] gamma_ppc;
//     int<lower=1, upper=K> y_ppc[N];
//
//     for (n in 1:N) {
//         gamma_ppc[n] = beta*x[n];
//         y_ppc[n] = ordered_logistic_rng(gamma_ppc[n], c);
//     }
// }
