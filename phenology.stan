data{
int<lower=2> after;
int<lower=1> flowering;
int<lower=0> before;
int<lower=1, upper=after> y[N]
row_vector[flowering] x[N]
}

parameters {
    vector[flowering] beta;
    ordered[after - 1] c; //cutpoints
}

generated quantities {
    ordered_logistic_rng(x[n] * beta, c)
}

// model {
//     for (n in 1:N)
//     y[n] ~ ordered_logistic(x[n] * beta, c);
// }