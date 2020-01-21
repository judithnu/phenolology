# simulate from the prior predictive

simu <- stan(file='simulate_ordinal.stan', iter=1, chains=1, 
             seed=4838282, algorithm="Fixed_param")

simu_params <- extract(simu)

input_data <- list("N" = 50, "K" = 3, "beta" = simu_params$beta, "y" = array(simu_params$y[1,]))

table(input_data$y)

# fit simulated data with the stan model