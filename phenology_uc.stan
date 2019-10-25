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
  
  vector[Nsite] z_site;
  vector[Nprovenance] z_prov;
  vector[Nclone] z_clone;
  vector[Nyear] z_year; //uncentered
  
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
  
  z_site ~ normal( 0 , 1 );
  z_prov ~ normal( 0 , 1 );
  z_clone ~ normal( 0 , 1 );
  z_year ~ normal( 0 , 1 ); 
  
  // model
  
  for ( i in 1:N ) {
    phi[i] = forcing[i] * ( beta + 
    z_site[SiteID[i]]*sigma_site + z_prov[ProvenanceID[i]]*sigma_prov + 
    z_clone[CloneID[i]]*sigma_clone + z_year[YearID[i]]*sigma_year );
  }
  for ( i in 1:N ) state[i] ~ ordered_logistic( phi[i] , kappa );
}

generated quantities{
  //DECLARE
  //re-centered vars
  vector[Nsite] b_site;
  vector[Nprovenance] b_prov;
  vector[Nclone] b_clone;
  vector[Nyear] b_year;
  
  //mean effects
  real b_site_mean;
  real b_clone_mean;
  real b_year_mean;
  real b_prov_mean;

  //ppc y_rep : uncomment to generate ppc yrep
  // vector[N] phi;
  // vector[N] state_rep;
  
  //DEFINE
  
  // recalculate b parameters that were un-centered
  b_site = z_site*sigma_site;
  b_prov = z_prov*sigma_prov;
  b_clone = z_clone*sigma_clone;
  b_year = z_year*sigma_year;
  
  // calculate mean effects across groups
  b_site_mean = mean(b_site);
  b_prov_mean = mean(b_prov);
  b_clone_mean = mean(b_clone);
  b_year_mean = mean(b_year);
  
  // simulate data for model testing
  //   for ( i in 1:N ) {
    //     phi[i] = forcing[i] * (beta + b_site[SiteID[i]] + b_prov[ProvenanceID[i]] + b_clone[CloneID[i]] + b_year[YearID[i]]);
    //   }
  //   
    //   for (i in 1:N) {
      //     state_rep[i] = ordered_logistic_rng(phi[i], kappa);
      //   }
}

