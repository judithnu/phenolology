// simulate phenology data
// 3 phenophases (1,2,3)
// heatsum thresholds and transition rate calculated

/// generate fake data. 30 clones with 2 ramets with 10 observations each
transformed data {
    int<lower=2> K = 3; // possible number of phenophases
    int<lower=0> N = 3000; // number of phenophase observations
    int<lower=1> N_clone = 30; // categories/clones
    ordered[2] cutpoints = [11, 18]';
    // int<lower=1, upper=K> y[N]; // outcomes from the model
    // row_vector[H] heatsum[N]; //heatsums
}

generated quantities {
    row_vector[N] heatsum; // vector for heatsum
    real phenophase[N]; // Simulated phenophases
    vector[N_clone] beta; //Simulated betas
    vector[N_clone] clone; //Clone ids

    // Construct a heatsum covariate
    for (n in 1:N) {
        heatsum[n] = uniform_rng(0,800);
    }

    // Construct beta parameter
    for (i in 1:N_clone) {
        beta[i] = beta_rng(.5,5);
    }

    // Simulate data from observational model
    for (n in 1:N) {
      for (i in 1:N_clone) {
        phenophase[n] = ordered_logistic_rng(heatsum[n] * beta[i], cutpoints);
        clone[i] = i;
}}
}
