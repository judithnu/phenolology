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
    int CloneID[N];
    int YearID[N];
    int SiteID[N];
    int ProvenanceID[N];
}
parameters{
    //intercepts
    vector[Nprovenance] a_prov;
    vector[Nsite] a_site;
    vector[Nyear] a_year;
    vector[Nclone] a_clone;
    //slope effects
    vector[Nsite] beta_site;
    vector[Nprovenance] beta_prov;
    real<lower=0, upper=1> beta;
    //cutpoints
    positive_ordered[K-1] kappa;
    //variance for adaptive priors
    real<lower=0> prov_sigma;
    real<lower=0> site_sigma;
    real<lower=0> year_sigma;
    real<lower=0> clone_sigma;
    real<lower=0> site_nu;
    real<lower=0> prov_nu;
}
transformed parameters{
    real kappa_diff; //difference between cutpoints
    vector[N] alpha_tot; //all intercept effects
    vector[N] beta_tot; //total slope
    matrix[N,K-1] h50; //inflection points
    vector[N] h50diff; //difference between inflection points

    kappa_diff = kappa[K-1] - kappa[K-2];
    alpha_tot = a_prov[ProvenanceID] + a_site[SiteID] + a_year[YearID] + a_clone[CloneID];
    beta_tot = beta + beta_site[SiteID] + beta_prov[ProvenanceID];

    //inflection points
    for (i in 1:N) {
        h50[i,K-2] = (kappa[K-2] + alpha_tot[i])/beta_tot[i];
}
    for (i in 1:N) {
        h50[i,K-1] = (kappa[K-1] + alpha_tot[i])/beta_tot[i];
    }
    h50diff = h50[,K-1] - h50[,K-2];
}


model{
    vector[N] phi; //outcomes
    //hyper priors
    prov_nu ~ exponential( 6 );
    site_nu ~ exponential( 6 );
    clone_sigma ~ exponential( 1.5 );
    year_sigma ~ exponential( 1.5 );
    site_sigma ~ exponential( 1.5 );
    prov_sigma ~ exponential( 1.5 );
    //fixed priors
    kappa_diff ~ gamma( 7.5 , 1 );
    kappa[K-2] ~ gamma( 5 , 1 );
    beta ~ beta(.5, 5 );
    //adaptive priors
    beta_prov ~ normal( 0 , prov_nu );
    beta_site ~ normal( 0 , site_nu );
    a_clone ~ normal( 0 , clone_sigma );
    a_year ~ normal( 0 , year_sigma );
    a_site ~ normal( 0 , site_sigma );
    a_prov ~ normal( 0 , prov_sigma );
    //linear model
    for ( i in 1:N ) {
        phi[i] = a_prov[ProvenanceID[i]] + a_site[SiteID[i]] + a_year[YearID[i]] + a_clone[CloneID[i]] + (beta + beta_site[SiteID[i]] + beta_prov[ProvenanceID[i]]) * forcing[i];
    }
    //likelihood
    for ( i in 1:N ) state[i] ~ ordered_logistic( phi[i] , kappa );
}


