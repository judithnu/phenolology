// data{
// int<lower=2> after;
// int<lower=1> flowering;
// int<lower=0> before;
// int<lower=1, upper=after> y[N]
// row_vector[flowering] x[N]
// }
//
// parameters {
//     vector[flowering] beta;
//     ordered[after - 1] c; //cutpoints
// }

data{
    int<lower=2> K; // number of possible phenophases
    int<lower=0> N; // number of phenophases observations
    int<lower=1, upper=K> y[N]; //phenophase outcomes (response)

    int<lower=0> D; //number of heatsum observations
    row_vector[D] x[N]; //heatsums
}

parameters{
    vector[D] beta;
    ordered[K-1] c; //cutpoints mean
}

model {
    for (n in 1:N)
      y[n] ~ ordered_logistic(x[n] * beta, c);
}


// generated quantities {
//     // Simulate model configuration from prior model
//
//     //Simulate data from observational model
//     int y[N] = rep_array(0,N);
//     for (n in 1:N)
//       y[n] = ordered_logistic_rng(theta);
// }

// model {
//     for (n in 1:N)
//     y[n] ~ ordered_logistic(x[n] * beta, c);
// }
