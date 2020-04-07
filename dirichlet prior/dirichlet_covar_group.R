# Take dirichlet_learn.R's simple example and add a covariate


library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())

# fauxdat

N = 5000
G = 7
K = 3

# covariate
x <- rnorm(N, mean = 12, sd=2)


# group effects
gbeta_mu <- 0
gbeta_sd <- 0.25
gbeta_vec <- rnorm(7, mean=gbeta_mu, sd=gbeta_sd)
gbeta <- sample(gbeta_vec, size=N, replace=TRUE)
gid <- as.numeric(as.factor(gbeta))


# simulate group level effects

input_data_for_simulation <- list("N" = N, "G"= G, "K"=3, "x" = x, "gbeta" = gbeta)

simu <- stan(file='dirichlet prior/simulate_ordered_covar_group.stan', iter=1, chains=1,
             algorithm="Fixed_param", data=input_data_for_simulation)

simu_params <- extract(simu)

input_data_for_model <- list("N" = N, "K" = K, "G"= G, "x" = x, "GID"=gid, "y" = array(simu_params$y[1,]))
round(unique(simu_params$c), 2)

table(input_data_for_model$y)
plot(input_data_for_model$x, input_data_for_model$y)

#Fit the simulated data using the model with a dirichlet prior

# need to make a version of this file that has a covariate model
fit <- stan(file='dirichlet prior/ordered_logistic_induced_with_covar_and_group.stan', data=input_data_for_model, chains=1)

params <- data.frame(extract(fit))
summary(params)

hist(params$c.1/(params$beta+params$gbeta.1))
hist(params$c.2/(params$beta+params$gbeta.1))
