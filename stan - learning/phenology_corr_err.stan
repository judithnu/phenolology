data{
    int Phenophase_Derived[10891];
    int YearID[10891];
    int CloneID[10891];
    int ProvenanceID[10891];
    int SiteID[10891];
    vector<lower=0>[10891] sum_forcing;
}
parameters{
    vector<lower=0>[10891] sum_forcing_true;
    real<lower=0> beta;
    positive_ordered[2] kappa;
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
    real<lower=0> sigma_err;
}
model{
    vector[10891] phi;
    year_sigma ~ exponential( 2 );
    clone_sigma ~ exponential( 2 );
    Rhop ~ lkj_corr( 2 );
    Rhos ~ lkj_corr( 2 );
    sigma_prov ~ exponential( 1.5 );
    bp ~ normal( 0 , 0.5 );
    ap ~ normal( 0 , 1 );
    {
    vector[2] YY[6];
    vector[2] MU;
    MU = [ ap , bp ]';
    for ( j in 1:6 ) YY[j] = [ a_prov[j] , b_prov[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rhop , sigma_prov) );
    }
    sigma_site ~ exponential( 1.5 );
    bs ~ normal( 0 , 0.5 );
    as ~ normal( 0 , 1 );
    {
    vector[2] YY[7];
    vector[2] MU;
    MU = [ as , bs ]';
    for ( j in 1:7 ) YY[j] = [ a_site[j] , b_site[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rhos , sigma_site) );
    }
    a_year ~ normal( 0 , year_sigma );
    a_clone ~ normal( 0 , clone_sigma );
    kappa ~ gamma( 7.5 , 1 );
    beta ~ exponential( 2 );

    for ( i in 1:10891 ) {
        sum_forcing_true[i] ~ gamma( 7.5, 2 );
        sum_forcing[i] ~ normal(sum_forcing_true[i],.25);
        phi[i] = a_site[SiteID[i]] + a_prov[ProvenanceID[i]] + a_clone[CloneID[i]] + a_year[YearID[i]] + (beta + b_site[SiteID[i]] + b_prov[ProvenanceID[i]]) * sum_forcing_true[i];
    }
    for ( i in 1:10891 ) Phenophase_Derived[i] ~ ordered_logistic( phi[i] , kappa );
}


