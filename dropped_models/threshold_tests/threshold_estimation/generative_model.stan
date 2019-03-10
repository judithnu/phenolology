// simulate phenology data
// 3 phenophases (1,2,3)
// heatsum thresholds and transition rate calculated

/// generate fake data. 30 clones with 2 ramets with 10 observations each
transformed data {
    int<lower=2> K = 3; // possible number of phenophases
    int<lower=0> T = 50; // timesteps
    //int<lower=1> N_ind = 30; // number of individuals
    //int<lower=0> N = T*N; // number of phenophase observations
    ordered[2] cutpoints = [11, 18]';
    real<lower=0, upper=10> threshold = 4;
    real<lower=0, upper=1> beta = 0.2;
    // int<lower=1, upper=K> y[N]; // outcomes from the model
    // row_vector[H] heatsum[N]; //heatsums
}

generated quantities {
    row_vector[T] temp; //vector for temperatures
    row_vector[T] heat; //vector for heat
    row_vector[T] heatsum; // vector for heatsum
    real phenophase[T]; // Simulated phenophases
    int time[T]; //time
    // vector[N_clone] beta; //Simulated betas
    // vector[N_clone] clone; //Clone ids

    // Construct a heatsum covariate
        // temperatures
    for (t in 1:T) {
        temp[t] = uniform_rng(-5,30);
    }

    temp = sort_asc(temp);

        //daily heating
    for (t in 1:T) {
        if (temp[t] < threshold) {
          heat[t] = 0;
        }
        else if (temp[t] >= threshold) {
          heat[t] = temp[t] - threshold;
        }
    }
        //heatsums
    for (t in 1:T) {
        if (t == 1) {
          heatsum[t] = heat[t];
        }
        else if (t > 1) {
          heatsum[t] = heat[t] + heatsum[t-1];
        }
    }

    // // Construct beta parameter
    // for (i in 1:N_clone) {
    //     beta[i] = beta_rng(.5,5);
    // }

    // Simulate data from observational model
    for (t in 1:T) {
        phenophase[t] = ordered_logistic_rng(heatsum[t] * beta, cutpoints);
        time[t] = t;
    }
}

