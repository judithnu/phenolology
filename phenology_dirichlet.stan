functions {
  //from Michael Betancourt https://betanalpha.github.io/assets/case_studies/ordinal_regression.html
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}

data{
  //how many
  int<lower=1> N;
  int<lower=1> K;
  int Nsite;
  int Nprovenance;
  int Nclone;
  int Nyear;
  
  // data
  int<lower=1, upper=K> state[N];
  vector[N] forcing;
  int SiteID[N];
  int ProvenanceID[N];
  int CloneID[N];
  int YearID[N];
}

parameters{
  positive_ordered[K-1] kappa; //(Internal) cut points
  real<lower=0> beta; //population slope
  
  vector[Nsite] b_site;
  vector[Nprovenance] b_prov;
  vector[Nclone] b_clone;
  vector[Nyear] b_year;
  
  real<lower=0> sigma_site;
  real<lower=0> sigma_prov;
  real<lower=0> sigma_clone;
  real<lower=0> sigma_year;
}


model{
  vector[N] phi;
  beta ~ exponential(2);
  //kappa ~ gamma( 20, 1 ) ;
  kappa ~ induced_dirichlet(rep_vector(1,K),0);
  
  // adaptive priors on effects
  sigma_site ~ exponential( 5 );
  sigma_prov ~ exponential( 5 );
  sigma_clone ~ exponential( 5 );
  sigma_year ~ exponential( 5 );
  
  b_site ~ normal( 0 , sigma_site );
  b_prov ~ normal( 0 , sigma_prov );
  b_clone ~ normal( 0 , sigma_clone );
  b_year ~ normal( 0, sigma_year );

  
  // model
  
  for ( i in 1:N ) {
    phi[i] = forcing[i] * ( beta + b_site[SiteID[i]] + b_prov[ProvenanceID[i]] + b_clone[CloneID[i]] + b_year[YearID[i]] );
  } //linear model
  for ( i in 1:N ) state[i] ~ ordered_logistic( phi[i] , kappa ); //obs model
}

// generated quantities{
//   //DECLARE
//   
//   //mean effects
//   real b_site_mean;
//   real b_clone_mean;
//   real b_year_mean;
//   real b_prov_mean;
// 
//   //ppc y_rep : uncomment to generate ppc yrep
//   vector[N] phi;
//   vector[N] state_rep;
//   
//   //DEFINE
//   
//   // calculate mean effects across groups
//   b_site_mean = mean(b_site);
//   b_prov_mean = mean(b_prov);
//   b_clone_mean = mean(b_clone);
//   b_year_mean = mean(b_year);
//   
//   // simulate data for model testing
//     for ( i in 1:N ) {
//       phi[i] = forcing[i] * (beta + b_site[SiteID[i]] + b_prov[ProvenanceID[i]] + b_clone[CloneID[i]] + b_year[YearID[i]]);
//     }
// 
//     for (i in 1:N) {
//       state_rep[i] = ordered_logistic_rng(phi[i], kappa);
//     }
// }
// 
