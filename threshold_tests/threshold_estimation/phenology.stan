// ordered logistic model with threshold

functions {
    // function declarations and definitions
    real heatsumcalculator(time, timesteps, temperature, threshold) {
        for (t in 1:T) {
          if (temperature[t] < threshold) {
            heat[t] = 0;
          }
          else if (temperature[t] >= threshold) {
            heat[t] = temperature[t] - threshold;
          }
        }

    for (t in 1:T) {
        if (t == 1) {
          heatsum[t] = heat[t];
        }
        else if (t > 1) {
          heatsum[t] = heat[t] + heatsum[t-1];
        }
    }
    return heatsum;
    }
}

data{
    int<lower=2> K; // number of possible phenophases
    int<lower=0> T; // number of timesteps
    int<lower=0> N; // number of indivs
    int<lower=1> time[T];
    real<lower=-30, upper=30> temp[T]; //temperature observations
    matrix<lower=1, upper=K>[T,N] phenophase; //phenophase outcomes (response)
}

parameters{
    real<lower=0,upper=10> threshold;
    //real<lower=0, upper=1> beta[N_clone];
    ordered[K-1] c; //how many cutpoints are there
    //real<lower=0.001, upper=2> shape1;
    //real<lower=1, upper=5> shape2;
    real<lower=0, upper=1> beta;
    real heat[T];
    real heatsum[T];
}


model {
    c ~ normal(5, 10); //cutpoints prior
    beta ~ beta(1,1); //beta prior
    threshold ~ normal(5,1);

//calculate heatsum
    heatsum =

    for (n in 1:N) {
        for (t in 1:T) {
            phenophase[t,n] ~ ordered_logistic(heatsum[t] * beta, c);
    }
    }

}

