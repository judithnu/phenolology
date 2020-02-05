#meta for dirichlet_prior

This folder contains files where I figure out the dirichlet prior in an ordered logistic model.

The following files directly follow [Mike Betancourt's case study](https://betanalpha.github.io/assets/case_studies/ordinal_regression.html)

* `dirichlet_learn.R`
* `simulate_ordered.stan`
* `ordered_logistic_induced.stan`

Next, I attempted adding a simple covariate $\beta * x$ where both $\beta$ and $x$ are drawn from the standard normal.

* `dirichlet_covar.R`
* `simulated_ordered_covar.stan`