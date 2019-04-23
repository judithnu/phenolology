data{
    int Phenophase_Derived[10963];
    vector[10963] sum_forcing;
    int CloneID[10963];
    int YearID[10963];
    int SiteID[10963];
    int ProvenanceID[10963];
}
parameters{
    vector[6] a_prov;
    vector[7] a_site;
    vector[15] a_year;
    vector[259] a_clone;
    vector[7] beta_site;
    vector[6] beta_prov;
    real<lower=0> beta;
    positive_ordered[2] kappa;
    real<lower=0> prov_sigma;
    real<lower=0> site_sigma;
    real<lower=0> year_sigma;
    real<lower=0> clone_sigma;
    real<lower=0> site_nu;
    real<lower=0> prov_nu;

}
transformed parameters{
    vector[10963] alpha_tot;
    vector[10963] beta_tot;
    vector[10963] h1;
    vector[10963] h2;
    real kappa_diff;

    kappa_diff = kappa[2] - kappa[1];

    //inflection points

    for (i in 1:10963) {
        h1[i] = (kappa[1] + alpha_tot[i])/beta_tot[i];
    }
    for (i in 1:10963) {
        h1[i] = (kappa[2] + alpha_tot[i])/beta_tot[i];
    }

    //total alphas
    for ( i in 1:10963 ) {
        alpha_tot[i] = a_prov[ProvenanceID[i]] + a_site[SiteID[i]] + a_year[YearID[i]] + a_clone[CloneID[i]];
    }

    //total betas

    for ( i in 1:10963 ) {
        beta_tot[i] = beta + beta_site[SiteID[i]] + beta_prov[ProvenanceID[i]];
    }

}
model{
    vector[10963] phi;

    //fixed priors
    prov_nu ~ exponential( 2 );
    site_nu ~ exponential( 2 );

    clone_sigma ~ exponential( 1.5 );
    year_sigma ~ exponential( 1.5 );
    site_sigma ~ exponential( 1.5 );
    prov_sigma ~ exponential( 1.5 );

    beta ~ beta( 2, 2 );
    //alpha_tot ~ normal(0,1);
    kappa_diff ~ gamma(7.5 , 1 ); //prior on difference between cutpoints
    kappa[1] ~ gamma(7.5 , 1 ); //nail down the first cutpoint


    // adaptive priors
    beta_prov ~ uniform( -1, 1 );
    beta_site ~ uniform( -1, 1 );
    a_clone ~ normal( 0 , clone_sigma );
    a_year ~ normal( 0 , year_sigma );
    a_site ~ normal( 0 , site_sigma );
    a_prov ~ normal( 0 , prov_sigma );

    //model
    for ( i in 1:10963 ) {
        phi[i] = a_prov[ProvenanceID[i]] + a_site[SiteID[i]] + a_year[YearID[i]] + a_clone[CloneID[i]] + (beta + beta_site[SiteID[i]] + beta_prov[ProvenanceID[i]]) * sum_forcing[i];
    }
    for ( i in 1:10963 ) Phenophase_Derived[i] ~ ordered_logistic( phi[i] , kappa );
}


