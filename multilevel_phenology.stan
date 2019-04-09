data{
    int N; // number of observations

    int K; //number of phenophases
    int Nsex; //number of sexes
    int Nclone; //number of clones
    int Nprovenance; //number of provenances
    int Nsite; //number of sites

    int SiteID[N];
    int CloneID[N];
    int ProvenanceID[N];
    int SexID[N];

    vector[N] forcing;
    int state[N];

    real L; // lower bound for cutpoints difference
}

parameters{
    vector<lower=0,upper=1>[Nsex] beta;
    positive_ordered[K-1] cutpoints;
    vector[Nsex] a_sex;
    vector[Nclone] a_clone;
    vector[Nprovenance] a_provenance;
    vector[Nsite] a_site;

    real<lower=0> sigma_sex;
    real<lower=0> sigma_clone;
    real<lower=0> sigma_provenance;
    real<lower=0> sigma_site;

}

transformed parameters{
    //declare
    real<lower=1>cdiff; // don't allow cutpoint collapse

    //define
    cdiff = cutpoints[2]-cutpoints[1];
}

model{
    vector[N] phi;

    //sample priors

    sigma_site ~ exponential( 1.5 );
    sigma_provenance ~ exponential( 1.5 );
    sigma_clone ~ exponential( 1.5 );

    a_sex ~ normal( 0, sigma_sex);
    a_site ~ normal( 0 , sigma_site );
    a_provenance ~ normal( 0 , sigma_provenance );
    a_clone ~ normal( 0 , sigma_clone );

    cutpoints[1] ~ gamma( 7.5 , 1 );
    cdiff ~ gamma(5,1) T[L,];

    beta ~ beta( 0.5 , 5 );

    //model
    for ( i in 1:N ) {
        phi[i] = beta[SexID[i]] * forcing[i] + a_sex[SexID[i]] + a_provenance[ProvenanceID[i]] + a_clone[CloneID[i]] + a_site[SiteID[i]];
    }
    for ( i in 1:N ) state[i] ~ ordered_logistic( phi[i] , cutpoints ); //consider vectorizing
}


