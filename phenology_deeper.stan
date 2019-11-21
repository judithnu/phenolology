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

  real mu_site;
  real mu_prov;
  real mu_clone;
  real mu_year;

  real<lower=0> sigma_site;
  real<lower=0> sigma_prov;
  real<lower=0> sigma_clone;
  real<lower=0> sigma_year;
}


model{
  vector[N] phi;
  beta ~ exponential(2);
  kappa ~ gamma( 20, 1 ) ;

  // adaptive priors on effects
  sigma_site ~ exponential( 5 );
  sigma_prov ~ exponential( 5 );
  sigma_clone ~ exponential( 5 );
  sigma_year ~ exponential( 5 );

  mu_site ~ normal(0,1);
  mu_prov ~ normal(0,1);
  mu_clone ~ normal(0,1);
  mu_year ~ normal(0,1);

  b_site ~ normal( mu_site , sigma_site );
  b_prov ~ normal( mu_prov , sigma_prov );
  b_clone ~ normal( mu_clone , sigma_clone );
  b_year ~ normal( mu_year, sigma_year );


  // model

  for ( i in 1:N ) {
    phi[i] = forcing[i] * ( beta + b_site[SiteID[i]] + b_prov[ProvenanceID[i]] + b_clone[CloneID[i]] + b_year[YearID[i]] );
  }
  for ( i in 1:N ) state[i] ~ ordered_logistic( phi[i] , kappa );
}

generated quantities{
  //DECLARE

  //mean effects
  real b_site_mean;
  real b_clone_mean;
  real b_year_mean;
  real b_prov_mean;

  //ppc y_rep : uncomment to generate ppc yrep
  vector[N] phi;
  vector[N] state_rep;

  //DEFINE

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

