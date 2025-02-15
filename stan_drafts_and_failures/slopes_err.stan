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
  vector<lower=0>[N] forcing_true;
    positive_ordered[K-1] kappa;
    real<lower=0> beta;
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
    beta ~ exponential(2);
    sigma_year ~ exponential( 3 );
    sigma_clone ~ exponential( 3 );
    sigma_prov ~ exponential( 3 );
    sigma_site ~ exponential( 3 );
    b_year ~ normal( 0 , sigma_year );
    b_clone ~ normal( 0 , sigma_clone );
    b_prov ~ normal( 0 , sigma_prov );
    b_site ~ normal( 0 , sigma_site );
    kappa ~ gamma( 7.5 , 1 );
    for ( i in 1:N ) {
      forcing_true[i] ~ normal(300,60);
      forcing[i] ~ normal(forcing_true[i],5);
        phi[i] = forcing[i] * (beta + b_site[SiteID[i]] + b_prov[ProvenanceID[i]] + b_clone[CloneID[i]] + b_year[YearID[i]]);
    }
    for ( i in 1:N ) state[i] ~ ordered_logistic( phi[i] , kappa );
}

