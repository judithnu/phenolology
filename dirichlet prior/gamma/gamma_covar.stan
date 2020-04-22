
data {
    int<lower=1> N;             // Number of observations
    int<lower=1> K;             // Number of ordinal categories

    int<lower=1, upper=K> y[N]; // Observed ordinals
    vector[N] x;                  // Covariate
    real<lower=1> shape; // shape parameter for gamma prior on cutpoints
    real<lower=0> rate; // rate parameter for exponential prior on beta
}

parameters {
    positive_ordered[K - 1] c; // (Internal) cut points
    real beta; // population level effect
}

model {
    vector[N] phi;

    // Prior model
    beta ~ exponential(rate);
    c ~ gamma(shape,1);


    // Observational model

  for (i in 1:N ) {
    phi[i] = beta * x[i];
    y[i] ~ ordered_logistic(phi[i], c);
    }

}

generated quantities {
    //vector[N] gamma_ppc;
    real<lower=0> h1;
    real<lower=h1> h2;
    //int<lower=1, upper=K> y_ppc[N];

    h1 = c[1]/beta;
    h2 = c[2]/beta;

    // for (n in 1:N) {
    //     gamma_ppc[n] = beta*x[n];
    //     y_ppc[n] = ordered_logistic_rng(gamma_ppc[n], c);
    // }
}
