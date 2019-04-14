data{
    int N; //number of observations
    int K; //number of possible states
    int state[N];
    vector[N] forcing; //predictor
}

parameters{
    positive_ordered[K-1] cutpoints; //cutpoints are positive and ordered because forcing units are never negative
    real<lower=0,upper=1> beta;
    real<lower=0> alpha;
}

transformed parameters {
    //declare params
    vector[K-1] h50; // transition points on forcing scale
    real<lower=0> cdiff; //difference between cutpoints
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
        //for eta
    beta ~ beta( 0.5 , 5 );
        //for cutpoints
    cutpoints[1] ~ exponential(0.3); // pin down first cutpoint
    alpha ~ lognormal(1,1);
    cdiff ~ gamma(alpha, 1); // difference between first and second cutpoint
    //model
    for ( i in 1:N ) {
        eta[i] = beta * forcing[i];
    }
        state ~ ordered_logistic( eta , cutpoints );
}
