data{
    //how many
    int N;
    int K;
    int Nsite;
    int Nprovenance;
    int Nclone;
    int Ntree;
    int Nyear;
    // data
    int state[N];
    vector[N] forcing;
    int SiteID[N];
    int ProvenanceID[N];
    int CloneID[N];
    int TreeID[N];
    int YearID[N];
}
parameters{
    positive_ordered[K-1] kappa;
    real<lower=0> beta;
    vector[Nsite] b_site;
    vector[Nprovenance] b_prov;
    vector[Nclone] b_clone;
    vector[Ntree] b_tree;
    vector[Nyear] b_year;
    real<lower=0> sigma_site;
    real<lower=0> sigma_prov;
    real<lower=0> sigma_clone;
    real<lower=0> sigma_tree;
    real<lower=0> sigma_year;
}

transformed parameters{
  real b_clone_mean;
  real b_tree_mean;

  b_clone_mean = mean(b_clone);
  b_tree_mean = mean(b_tree);
}

model{
    vector[N] phi;
    //main effects
    kappa ~ gamma( 7.5 , 1 );
    beta ~ exponential(2);
    //site,prov,clone,tree,year effects
    sigma_site ~ exponential( 4 );
    sigma_prov ~ exponential( 4 );
    sigma_clone ~ exponential( 4 );
    sigma_tree ~ exponential( 4 ) ;
    sigma_year ~ exponential( 4 );

    b_site ~ normal( 0 , sigma_site );
    b_prov ~ normal( 0 , sigma_prov );
    b_clone ~ normal( 0 , sigma_clone );
    b_tree ~ normal(0, sigma_tree);
    b_year ~ normal( 0 , sigma_year );
    //eta
    for ( i in 1:N ) {
        phi[i] = forcing[i] * (beta + b_site[SiteID[i]] + b_prov[ProvenanceID[i]] + b_clone[CloneID[i]] + b_year[YearID[i]]);
    }
    //likelihood
    for ( i in 1:N ) state[i] ~ ordered_logistic( phi[i] , kappa );
}

