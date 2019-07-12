data{
    int<lower=0> N; //number of observations
    int<lower=2> K; //number of possible states

    int state[N]; //response
    vector[N] forcing; //predictor

    int<lower=0> Ngroup; //number of groups
    int<lower=1, upper=Ngroup> group[N]; //group
}

parameters{
    vector[K-1] cutpoints[Ngroup];
    vector<lower=0>[K-1] cutpoints_sigma;
    vector[K-1] cutpoints_mean;
    real<lower=0,upper=1> beta;
  //  real<lower=0> alpha;
}
//
// transformed parameters {
//     //declare params
//     // vector[K-1] h50; // transition points on forcing scale
//     // real<lower=0> cdiff;
//     // //define params
//     // for ( i in 1:2 ) {
//     //     h50[i] = cutpoints[i]/beta;
//     // }
//     //cdiff = cutpoints[2]-cutpoints[1];
// }

model{
    //declarations
    vector[N] eta;
    //priors
    cutpoints_sigma ~ exponential(1);
    cutpoints_mean ~ normal(0, 10);
    cutpoints[,1] ~ normal(cutpoints_mean[1], cutpoints_sigma[1]); //cutpoint 1 for all groups
    cutpoints[,2] ~ normal(cutpoints_mean[2], cutpoints_sigma[2]); //delta cutpoint 2 for all groups
    beta ~ beta( 0.5 , 5 );

    for (i in 1:N) {
        vector[K-1] cuts;

        //map to ordered sequence by Ordered Inverse Transform
        cuts = cutpoints[group[i]];
        cuts[2] = cuts[1] + exp(cuts[2]);

        eta[i] = beta * forcing[i];
        state[i] ~ ordered_logistic( eta[i] , cuts );
    }
  //  alpha ~ normal(10,2);
  //  cutpoints[1] ~ normal(5,5); // pin down first cutpoint
    //model
    //cdiff ~ gamma(alpha, 1); //model difference between first and second cutpoint
    // for ( i in 1:N ) {
    //     eta[i] = beta * forcing[i];
    // }
    // for ( i in 1:N ) {
    //     state[i] ~ ordered_logistic( eta[i] , cuts );
    // }
}


