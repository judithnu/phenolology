data{
    int Phenophase_Derived[10891];
    vector[10891] sum_forcing;
    int YearID[10891];
    int CloneID[10891];
    int ProvenanceID[10891];
    int SiteID[10891];
}
parameters{
    real<lower=0> beta; //progress cannot go backward
    positive_ordered[2] kappa; //finishing flowering cannot occur before flowering
    vector[260] a_clone;
    vector[15] a_year;
    vector[7] b_site;
    vector[7] a_site;
    real as;
    real bs;
    vector<lower=0>[2] sigma_site;
    vector[6] b_prov;
    vector[6] a_prov;
    real ap;
    real bp;
    vector<lower=0>[2] sigma_prov;
    corr_matrix[2] Rhos;
    corr_matrix[2] Rhop;
    real<lower=0> clone_sigma;
    real<lower=0> year_sigma;
}

model{
    vector[10891] phi;
    //fixed priors 
    kappa ~ gamma( 7.5 , 2 ); //positive and over entire distribution of forcing temperatures
    beta ~ exponential(2); //must be positive
    //adaptive priors for intercepts
    year_sigma ~ exponential( 2 );
    clone_sigma ~ exponential( 2 );
    a_year ~ normal( 0 , year_sigma );
    a_clone ~ normal( 0 , clone_sigma );
    //correlation between site slope and intercept - provenance
    Rhop ~ lkj_corr( 2 );
    sigma_prov ~ exponential( 1.5 );
    bp ~ normal(0, .25);
    ap ~ normal( 0, 1 );
    { 
    vector[2] YY[6];
    vector[2] MU;
    MU = [ ap , bp ]';
    for ( j in 1:6 ) YY[j] = [ a_prov[j] , b_prov[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rhop , sigma_prov) );
    }
    //correlation between site slope and intercept - site
    Rhos ~ lkj_corr( 2 );
    sigma_site ~ exponential( 1.5 );
    bs ~ normal(0, .25) ;
    as ~ normal( 0 , 1 );
    {
    vector[2] YY[7];
    vector[2] MU;
    MU = [ as , bs ]';
    for ( j in 1:7 ) YY[j] = [ a_site[j] , b_site[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rhos , sigma_site) );
    }
    for ( i in 1:10891 ) {
        phi[i] = a_site[SiteID[i]] + a_prov[ProvenanceID[i]] + a_clone[CloneID[i]] + a_year[YearID[i]] + (beta+ b_site[SiteID[i]] + b_prov[ProvenanceID[i]]) * sum_forcing[i];
    }
    for ( i in 1:10891 ) Phenophase_Derived[i] ~ ordered_logistic( phi[i] , kappa );
}

generated quantities {
  int<lower=0, upper=3> state_exp[10891];
  real phi_exp[10891];
  
  for (i in 1:10891) {
    phi_exp[i] = a_site[SiteID[i]] + a_prov[ProvenanceID[i]] + a_clone[CloneID[i]] + a_year[YearID[i]] + (beta+ b_site[SiteID[i]] + b_prov[ProvenanceID[i]]) * sum_forcing[i];
    state_exp[i] = ordered_logistic_rng( phi_exp[i], kappa);
  }
}
