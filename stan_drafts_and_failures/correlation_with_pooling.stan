data{
  int Phenophase_Derived[22416];
  vector[22416] forcing_accum;
  int YearID[22416];
  int CloneID[22416];
  int SiteID[22416];
  int ProvenanceID[22416];
  int SexID[22416];
}
parameters{
  real<lower=0,upper=1> b;
  ordered[2] kappa;
  vector[7] a_site;
  vector[260] a_clone;
  vector[15] a_year;
  vector[2] b_sex;
  vector[2] a_sex;
  real as;
  real bs;
  vector<lower=0>[2] sigma_sex;
  vector[6] b_prov;
  vector[6] a_prov;
  real ap;
  real bp;
  vector<lower=0>[2] sigma_prov;
  corr_matrix[2] Rhos;
  corr_matrix[2] Rhop;
  real<lower=0> site_sigma;
  real<lower=0> clone_sigma;
  real<lower=0> year_sigma;
}
model{
  vector[22416] phi;
  vector[22416] alphas;
  vector[22416] beta;
  year_sigma ~ exponential( 2 );
  clone_sigma ~ exponential( 2 );
  site_sigma ~ exponential( 2 );
  Rhop ~ lkj_corr( 2 );
  Rhos ~ lkj_corr( 2 );
  sigma_prov ~ exponential( 1 );
  bp ~ normal( 0 , 0.1 );
  ap ~ normal( 0 , 1 );
  {
    vector[2] YY[6];
    vector[2] MU;
    MU = [ ap , bp ]';
    for ( j in 1:6 ) YY[j] = [ a_prov[j] , b_prov[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rhop , sigma_prov) );
  }
  sigma_sex ~ exponential( 1 );
  bs ~ normal( 0 , 0.1 );
  as ~ normal( 0 , 1 );
  {
    vector[2] YY[2];
    vector[2] MU;
    MU = [ as , bs ]';
    for ( j in 1:2 ) YY[j] = [ a_sex[j] , b_sex[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rhos , sigma_sex) );
  }
  a_year ~ normal( 0 , year_sigma );
  a_clone ~ normal( 0 , clone_sigma );
  a_site ~ normal( 0 , site_sigma );
  kappa ~ gamma(7.5,1);
  b ~ beta( 0.5 , 5 );
  for ( i in 1:22416 ) {
    beta[i] = b + b_sex[SexID[i]] + b_prov[ProvenanceID[i]];
  }
  for ( i in 1:22416 ) {
    alphas[i] = a_sex[SexID[i]] + a_site[SiteID[i]] + a_prov[ProvenanceID[i]] + a_clone[CloneID[i]] + a_year[YearID[i]];
  }
  for ( i in 1:22416 ) {
    phi[i] = alphas[i] + beta[i] * forcing_accum[i];
  }
  for ( i in 1:22416 ) Phenophase_Derived[i] ~ ordered_logistic( phi[i] , kappa );
}
