data{
    int N; //number of observations
    int K; //number of possible states
    int state[N];
    vector[N] forcing;
}

parameters{
    positive_ordered[K-1] cutpoints; //cutpoints are positive and ordered because forcing units are never negative
    real<lower=0,upper=1> beta;
    real<lower=0> shape;
    real<lower=0> scale;
}

transformed parameters {
    real<lower=0> cutdiff;
    vector[K-1] h50;
    for (i in 1:2) {
        h50[i] = cutpoints[i]/beta;
    }
}

model{
    vector[N] eta;
    beta ~ beta( 0.5 , 5 );
    shape ~ gamma(7.5, 1);
    scale ~ gamma(1, 2);
    cutpoints[2]-cutpoints[1] ~ gamma(shape, scale);
    for ( i in 1:N ) {
        eta[i] = beta * forcing[i];
    }
    for ( i in 1:N ) {
        state[i] ~ ordered_logistic( eta[i] , cutpoints );
    }
}
