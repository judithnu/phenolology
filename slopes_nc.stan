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
    vector[Nprovenance] z_prov;
    vector[Nclone] b_clone;
    vector[Nyear] z_year;
    real<lower=0> sigma_site;
    real<lower=0> sigma_prov;
    real<lower=0> sigma_clone;
    real<lower=0> sigma_year;
}


model{
    vector[N] phi;
    beta ~ exponential(2);
    kappa ~ gamma( 20 , 1 );

    // adaptive priors on effects
    sigma_site ~ exponential( 5 );
    sigma_prov ~ exponential( 5 );
    sigma_clone ~ exponential( 5 );
    sigma_year ~ exponential( 5 );
    b_site ~ normal( 0 , sigma_site );
    z_prov ~ normal( 0 , 1 );
    b_clone ~ normal( 0 , sigma_clone );
    z_year ~ normal( 0 , 1 );

    // model

    for ( i in 1:N ) {
        phi[i] = forcing[i] * (beta + b_site[SiteID[i]] + z_prov[ProvenanceID[i]]*sigma_prov + b_clone[CloneID[i]] + z_year[YearID[i]]*sigma_year);
    }
    for ( i in 1:N ) state[i] ~ ordered_logistic( phi[i] , kappa );
}

generated quantities{

  //centered vars
  vector[Nprovenance] b_prov;
  vector[Nyear] b_year;

  //mean effects
  real b_site_mean;
  real b_clone_mean;
  real b_year_mean;
  real b_prov_mean;

  //ppc y_rep
  vector[N] phi;
  vector[N] state_rep;

  // recalculate b parameters that were un-centered
  b_prov = z_prov*sigma_prov;
  b_year = z_year*sigma_year;

  // calculate mean effects across groups
  b_site_mean = mean(b_site);
  b_prov_mean = mean(b_prov);
  b_clone_mean = mean(b_clone);
  b_year_mean = mean(b_year);

  // simulate data for model testing
  for ( i in 1:N ) {
        phi[i] = forcing[i] * (beta + b_site[SiteID[i]] + b_prov[ProvenanceID[i]] + b_clone[CloneID[i]] + b_year[YearID[i]]);
    }

  for (i in 1:N) {
    state_rep[i] = ordered_logistic_rng(phi[i], kappa);
  }
}


