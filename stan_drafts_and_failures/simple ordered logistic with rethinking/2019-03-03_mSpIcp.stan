data{
    int SiteID[9992];
    int OrchardID[9992];
    int Orchard[9992];
    vector[9992] Heat;
    int Clone[9992];
    int Index[9992];
    int DoY[9992];
    int Year[9992];
    int Phenophase_Derived[9992];
    int CloneID[9992];
    vector[9992] Heatsum;
    int ProvenanceID[9992];
}
parameters{
    vector<lower=0,upper=1>[4] b_provenance;
    ordered[2] cutpoints;
    vector[108] a_clone;
    vector[4] a_provenance;
    real<lower=0> sigma_clone;
    real<lower=0> sigma_provenance;
}
model{
    vector[9992] phi;
    sigma_provenance ~ cauchy( 0 , 1 );
    sigma_clone ~ cauchy( 0 , 1 );
    a_provenance ~ normal( 0 , sigma_provenance );
    a_clone ~ normal( 0 , sigma_clone );
    cutpoints ~ normal( 197 , 176 );
    b_provenance ~ beta( 0.5 , 5 );
    for ( i in 1:9992 ) {
        phi[i] = b_provenance[ProvenanceID[i]] * Heatsum[i] + a_clone[CloneID[i]] + a_provenance[ProvenanceID[i]];
    }
    for ( i in 1:9992 ) Phenophase_Derived[i] ~ ordered_logistic( phi[i] , cutpoints );
}


