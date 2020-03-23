# Take dirichlet_learn.R's simple example and add a covariate


library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())


simu <- stan(file='dirichlet prior/simulate_ordered_covar.stan', iter=1, chains=1,
             algorithm="Fixed_param")

simu_params <- extract(simu)

input_data <- list("N" = 500, "K" = 3, "x"= array(simu_params$x[1,]), "y" = array(simu_params$y[1,]))
round(unique(simu_params$c), 2)

table(input_data$y)
plot(input_data$x, input_data$y)

#Fit the simulated data using the model with a dirichlet prior

# need to make a version of this file that has a covariate model
fit <- stan(file='dirichlet prior/ordered_logistic_induced_with_covar.stan', data=input_data, chains=1)

params <- data.frame(extract(fit))
summary(params)
