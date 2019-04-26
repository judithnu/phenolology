data{
    int Phenophase_Derived[10891];
    int OrchardID[10891];
    int YearID[10891];
    int CloneID[10891];
    int ProvenanceID[10891];
    int SiteID[10891];
    vector[10891] sum_forcing;
}
parameters{
    positive_ordered[2] kappa;
    real<lower=0> beta;
    vector[7] b_site;
    vector[6] b_prov;
    vector[260] b_clone;
    vector[15] b_year;
    vector[17] b_orch;
    real<lower=0> sigma_site;
    real<lower=0> sigma_prov;
    real<lower=0> sigma_clone;
    real<lower=0> sigma_year;
    real<lower=0> sigma_orch;
}
model{
    vector[10891] phi;
    sigma_orch ~ exponential( 2 );
    sigma_year ~ exponential( 2 );
    sigma_clone ~ exponential( 2 );
    sigma_prov ~ exponential( 2 );
    sigma_site ~ exponential( 2 );
    b_orch ~ normal( 0 , sigma_orch );
    b_year ~ normal( 0 , sigma_year );
    b_clone ~ normal( 0 , sigma_clone );
    b_prov ~ normal( 0 , sigma_prov );
    b_site ~ normal( 0 , sigma_site );
    beta ~ exponential( 3 );
    kappa ~ gamma( 7.5 , 2 );
    for ( i in 1:10891 ) {
        phi[i] = sum_forcing[i] * (beta + b_site[SiteID[i]] + b_prov[ProvenanceID[i]] + b_clone[CloneID[i]] + b_year[YearID[i]] + b_orch[OrchardID[i]]);
    }
    for ( i in 1:10891 ) Phenophase_Derived[i] ~ ordered_logistic( phi[i] , kappa );
}


