data {
    int<lower=1> K; //number of (real) phenophase categories (z)
    int<lower=1> V; //number of (observed) phenophase categories (y)
    int<lower=0> t; //number of phenophase observations
    int<lower=1, upper=V> w[T]; //phenophase observations
    int<lower=1, upper=K> z[T]; //phenophases
    vector<lower=0>[K] alpha; //transition prior (process model)
    vector<lower=0>[V] beta; //emission prior (obs model)
}

parameters {
    simplex[K] theta[K]; //transition probabilities (between states in process model)
    simplex[V] phi[K]; //emission probabilties (obs | hidden state)
}

model {
    for (k in 1:K)
        theta[k] ~ dirichlet(alpha); // transition matrix for process model
    for (k in 1:K)
        phi[k] ~ dirichlet(beta); // obs error?
    for (t in 1:T)
        w[t] ~ categorical(phi[z[t]]); // generate obs based on process
    for (t in 2:T)
        z[t] ~ categorical(theta[z[t-1]]); // process model
}
