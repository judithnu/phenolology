// ordered logistic model


data{
    int<lower=2> K; // number of possible phenophases
    int<lower=0> N; // number of phenophase observations
    int<lower=1, upper=K> phenophase[N]; //phenophase outcomes (response)

    real<lower=0> heatsum[N]; //heatsums
}

parameters{
    real<lower=0, upper=1> beta;
    ordered[K-1] c; //how many cutpoints are there
}

model {
    c ~ uniform(0, 500); //cutpoints prior
    for (n in 1:N)
      phenophase[n] ~ ordered_logistic(heatsum[n] * beta, c);
}

