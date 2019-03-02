// ordered logistic model


data{
    int<lower=2> K; // number of possible phenophases
    int<lower=0> N; // number of phenophase observations
    int<lower=1> N_clone; // number of clones (groups)
    int<lower=1, upper=K> phenophase[N]; //phenophase outcomes (response)

    real<lower=0> heatsum[N]; //heatsums
}

parameters{
    real<lower=0, upper=1> beta[N_clone];
    ordered[K-1] c; //how many cutpoints are there
    real<lower=0.001, upper=2> shape1;
    real<lower=1, upper=5> shape2;
}

model {
    c ~ normal(5, 10); //cutpoints prior
    for (i in 1:N_clone)
        beta[N_clone] ~ beta(shape1, shape2); //beta prior
    for (n in 1:N)
        phenophase[n] ~ ordered_logistic(heatsum[n] * beta[N_clone], c);
}

