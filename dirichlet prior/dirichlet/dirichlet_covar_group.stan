//Betancourt's ordinal logistic with a dirichlet prior
functions {
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

data {
    int<lower=1> N;             // Number of observations
    int<lower=2> K;             // Number of ordinal categories
    int<lower=2> G;             // Number of groups

    int<lower=1, upper=K> y[N]; // Observed ordinals
    vector[N] x;                  // Covariate

    int GID[N]; //Groups
}

parameters {
    positive_ordered[K - 1] c; // (Internal) cut points
    real beta; // population level effect

    vector[G] betag;
    real<lower=0> sigma_group;
}

model {
    vector[N] phi;

    // Prior model
    beta ~ exponential(3);
    c ~ induced_dirichlet(rep_vector(1, K), 0);
    //betag ~ normal(0,0.25);
    sigma_group ~ exponential(4);
    betag ~ normal(0, sigma_group);


    // Observational model

  for (i in 1:N ) {
    phi[i] = (beta * betag[GID[i]]) * x[i];
    y[i] ~ ordered_logistic(phi[i], c);
    }

}


// generated quantities {
//     vector[N] gamma_ppc;
//     int<lower=1, upper=K> y_ppc[N];
//
//     for (n in 1:N) {
//         gamma_ppc[n] = beta*x[n];
//         y_ppc[n] = ordered_logistic_rng(gamma_ppc[n], c);
//     }
// }
