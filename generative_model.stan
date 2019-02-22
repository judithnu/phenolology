// simulate phenology data
// 3 phenophases (1,2,3)
// heatsum thresholds and transition rate calculated

/// generate fake data
transformed data {
    int<lower=2> K = 3; // possible number of phenophases
    int<lower=0> N = 100; // number of phenophase observations
    //int<lower=1> H; // number of heatsum observations
    real beta = .06; // true beta
    ordered[2] cutpoints = [11, 18]';
    // int<lower=1, upper=K> y[N]; // outcomes from the model
    // row_vector[H] heatsum[N]; //heatsums
}

generated quantities {
    row_vector[N] heatsum; // vector for heatsum
    real phenophase[N]; // Simulated phenophases

    // Construct a heatsum covariate
    for (n in 1:N) {
        heatsum[n] = uniform_rng(0,500);
    }

    // Simulate data from observational model
    for (n in 1:N)
      phenophase[n] = ordered_logistic_rng(heatsum[n] * beta, cutpoints);
}

