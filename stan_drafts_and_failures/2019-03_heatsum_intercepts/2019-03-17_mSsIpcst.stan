data{
    int OrchardID[20741];
    int Orchard[20741];
    int Index[20741];
    int Year[20741];
    int DoY[20741];
    int Phenophase_Derived[20741];
    int TreeID_new[20741];
    int SiteID[20741];
    int CloneID[20741];
    int ProvenanceID[20741];
    vector[20741] Heatsum;
    int SexID[20741];
}
parameters{
    vector<lower=0,upper=1>[2] beta;
    ordered[2] cutpoints;
    vector[260] a_clone;
    vector[6] a_provenance;
    vector[7] a_site;
    vector[858] a_tree;
    real<lower=0> sigma_clone;
    real<lower=0> sigma_provenance;
    real<lower=0> sigma_site;
    real<lower=0> sigma_tree;
}
model{
    vector[20741] phi;
    sigma_tree ~ exponential( 1.5 );
    sigma_site ~ exponential( 1.5 );
    sigma_provenance ~ exponential( 1.5 );
    sigma_clone ~ exponential( 1.5 );
    a_tree ~ normal( 0 , sigma_tree );
    a_site ~ normal( 0 , sigma_site );
    a_provenance ~ normal( 0 , sigma_provenance );
    a_clone ~ normal( 0 , sigma_clone );
    cutpoints ~ normal( 197 , 176 );
    beta ~ beta( 0.5 , 5 );
    for ( i in 1:20741 ) {
        phi[i] = beta[SexID[i]] * Heatsum[i] + a_provenance[ProvenanceID[i]] + a_clone[CloneID[i]] + a_site[SiteID[i]] + a_tree[TreeID_new[i]];
    }
    for ( i in 1:20741 ) Phenophase_Derived[i] ~ ordered_logistic( phi[i] , cutpoints );
}


