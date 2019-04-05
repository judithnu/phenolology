data{
    int N; //number of observations
    int K; //number of possible states
    int state[N];
    vector[N] forcing;
}

parameters{
    positive_ordered[K-1] cutpoints; //cutpoints are positive and ordered because forcing units are never negative
    real<lower=0,upper=1> beta;
    real<lower=0> alpha;
}

transformed parameters {
    //declare params
    vector[K-1] h50; // transition points on forcing scale
    real<lower=0> cdiff;
    //define params
    for ( i in 1:2 ) {
        h50[i] = cutpoints[i]/beta;
    }
    cdiff = cutpoints[2]-cutpoints[1];
}

model{
    //declarations
    vector[N] eta;
    //priors
    beta ~ beta( 0.5 , 5 );
    alpha ~ normal(10,2);
    cutpoints[1] ~ normal(5,5); // pin down first cutpoint
    //model
    cdiff ~ gamma(alpha, 1); //model difference between first and second cutpoint
    for ( i in 1:N ) {
        eta[i] = beta * forcing[i];
    }
    for ( i in 1:N ) {
        state[i] ~ ordered_logistic( eta[i] , cutpoints );
    }
}
