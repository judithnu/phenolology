//male intercept only model
data{
    int Phenophase_Derived[10963];
    vector[10963] forcing_accum;
    int CloneID[10963];
    int YearID[10963];
    int SiteID[10963];
    int ProvenanceID[10963];
}
parameters{
    real<lower=0,upper=1> beta;
    positive_ordered[2] kappa;
    vector[6] a_prov;
    vector[7] a_site;
    vector[15] a_year;
    vector[259] a_clone;
    real<lower=0> prov_sigma;
    real<lower=0> site_sigma;
    real<lower=0> year_sigma;
    real<lower=0> clone_sigma;
}
model{
    vector[10963] phi;
    clone_sigma ~ exponential( 1.5 );
    year_sigma ~ exponential( 1.5 );
    site_sigma ~ exponential( 1.5 );
    prov_sigma ~ exponential( 1.5 );
    a_clone ~ normal( 0 , clone_sigma );
    a_year ~ normal( 0 , year_sigma );
    a_site ~ normal( 0 , site_sigma );
    a_prov ~ normal( 0 , prov_sigma );
    kappa ~ gamma( 7.5, 1 );
    beta ~ beta( 0.5 , 5 );
    for ( i in 1:10963 ) {
        phi[i] = a_prov[ProvenanceID[i]] + a_site[SiteID[i]] + a_year[YearID[i]] + a_clone[CloneID[i]] + beta * forcing_accum[i];
    }
    for ( i in 1:10963 ) Phenophase_Derived[i] ~ ordered_logistic( phi[i] , kappa );
}
