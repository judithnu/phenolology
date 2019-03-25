data{
    int OrchardID[22392];
    int Orchard[22392];
    int Index[22392];
    int Year[22392];
    int Phenophase_Derived[22392];
    int TreeID_new[22392];
    int SiteID[22392];
    int CloneID[22392];
    int ProvenanceID[22392];
    vector[22392] forcing_accum;
    int SexID[22392];
}
parameters{
    vector<lower=0,upper=1>[2] beta;
    positive_ordered[2] cutpoints;
    vector[259] a_clone;
    vector[6] a_provenance;
    vector[7] a_site;
    vector[854] a_tree;
    real<lower=0> sigma_clone;
    real<lower=0> sigma_provenance;
    real<lower=0> sigma_site;
    real<lower=0> sigma_tree;
}
model{
    vector[22392] phi;
    sigma_tree ~ exponential( 1.5 );
    sigma_site ~ exponential( 1.5 );
    sigma_provenance ~ exponential( 1.5 );
    sigma_clone ~ exponential( 1.5 );
    a_tree ~ normal( 0 , sigma_tree );
    a_site ~ normal( 0 , sigma_site );
    a_provenance ~ normal( 0 , sigma_provenance );
    a_clone ~ normal( 0 , sigma_clone );
    cutpoints ~ gamma( 7.5 , 1 );
    beta ~ beta( 0.5 , 5 );
    for ( i in 1:22392 ) {
        phi[i] = beta[SexID[i]] * forcing_accum[i] + a_provenance[ProvenanceID[i]] + a_clone[CloneID[i]] + a_site[SiteID[i]] + a_tree[TreeID_new[i]];
    }
    for ( i in 1:22392 ) Phenophase_Derived[i] ~ ordered_logistic( phi[i] , cutpoints );
}


