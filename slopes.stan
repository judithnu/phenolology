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
    int SiteID[N];
    int ProvenanceID[N];
    int CloneID[N];
    int YearID[N];
}
parameters{
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

transformed parameters{
  real b_site_mean;
  real b_prov_mean;
  real b_clone_mean;
  real b_year_mean;

  b_site_mean = mean(b_site);
  b_prov_mean = mean(b_prov);
  b_clone_mean = mean(b_clone);
  b_year_mean = mean(b_year);
}

model{
    vector[N] phi;
    beta ~ exponential(2);
    kappa ~ gamma( 7.5 , 1 );

    // adaptive priors on effects
    sigma_site ~ exponential( 6 );
    sigma_prov ~ exponential( 6 );
    sigma_clone ~ exponential( 6 );
    sigma_year ~ exponential( 6 );
    b_site ~ normal( 0 , sigma_site );
    b_prov ~ normal( 0 , sigma_prov );
    b_clone ~ normal( 0 , sigma_clone );
    b_year ~ normal( 0 , sigma_year );

    for ( i in 1:N ) {
        phi[i] = forcing[i] * (beta + b_site[SiteID[i]] + b_prov[ProvenanceID[i]] + b_clone[CloneID[i]] + b_year[YearID[i]]);
    }
    for ( i in 1:N ) state[i] ~ ordered_logistic( phi[i] , kappa );
}

