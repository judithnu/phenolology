data{
    //how many
    int N;
    int K;
    int Nsite;
    int Nprovenance;
    int Nclone;
    int Nyear;
    // data
    int state[N];
    vector[N] forcing;
    int CloneID[N];
    int YearID[N];
    int SiteID[N];
    int ProvenanceID[N];
}
parameters{
    positive_ordered[K-1] kappa;
    vector[Nsite] b_site;
    vector[Nprovenance] b_prov;
    vector[Nclone] b_clone;
    vector[Nyear] b_year;

    real<lower=0> sigma_site;
    real<lower=0> sigma_prov;
    real<lower=0> sigma_clone;
    real<lower=0> sigma_year;

}
model{
    vector[N] phi;
    sigma_year ~ exponential( 2 );
    sigma_clone ~ exponential( 2 );
    sigma_prov ~ exponential( 2 );
    sigma_site ~ exponential( 2 );
    b_year ~ normal( 0 , sigma_year );
    b_clone ~ normal( 0 , sigma_clone );
    b_prov ~ normal( 0 , sigma_prov );
    b_site ~ normal( 0 , sigma_site );
    kappa ~ gamma( 7.5 , 1 );
    for ( i in 1:N ) {
        phi[i] = forcing[i] * (b_site[SiteID[i]] + b_prov[ProvenanceID[i]] + b_clone[CloneID[i]] + b_year[YearID[i]]);
    }
    for ( i in 1:N ) state[i] ~ ordered_logistic( phi[i] , kappa );
}
generated quantities{
    vector[42] beta;
    beta = b_site[SiteID] + b_prov[ProvenanceID];
}

