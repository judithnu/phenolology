############################################################
# Initial setup
############################################################
library(ggplot2)
library(tidyr)
library(dplyr)

library(rstan)
library(bayesplot)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

util <- new.env()
source('~/Documents/classes/STANworkshop/material/1 - workflow/stan_utility.R', local=util)

############################################################
# 1. Conceptual Analysis -------------------------
############################################################
# I want to be able to predict when trees will flower and for how long. Flowering is controlled by heatsum, but may also have a genetic component. The relationship between heatsum and flowering is logistic. In, a logistic model of transition between phenophases with logit(-k(xi-h)) where xi is the heatsum for a given record, h describes the heatsum at which there is a 50% chance an individual tree has begun transitioned between phases (or the point at which 50% of the trees in a population have transitioned between phases). k describes how quickly trees transition from one to another. Trees can only go from 1 to 2 to 3 - they cannot go backwards. An ordered logistic distribution keeps transitions in order.

# Model is usually written as logit(ak - b*xi) where a is the "intercept" or cutpoint for each transition. to convert it to h, divide b by a.

############################################################
# 2. Define Observations ----------------------------
############################################################
phenology_and_heatsum_data <- read.csv('~/Documents/research_phenolology/data/phenology_heatsum.csv', stringsAsFactors = FALSE)

pd <- phenology_and_heatsum_data

ggplot(pd, aes(x=Heatsum, y=Phenophase_Derived, colour=Sex)) +
    geom_jitter(shape = 3, alpha = 0.5) +
    facet_wrap(.~Year)


hist(pd$Phenophase_Derived, xlab = "Phenophase")
# discrete proportion of each phenophase value

pd_narm <- pd[!is.na(pd$Phenophase_Derived),] #drop nas
pr_pd <- table(pd_narm$Phenophase_Derived)/nrow(pd_narm)

# cumsum converts to cumulative proportions
cum_pr_pd <- cumsum(pr_pd)

# plot
plot(c(1:3), cum_pr_pd, type = "b", xlab = "Phenophase", ylab = "cumulative proportion", ylim=c(0,1))

# cumulative proportion by heatsum
pr_by_heat <- table(pd$Phenophase_Derived, pd$Heatsum)


# This graph for lab meeting
ggplot(pd, aes(x = Heatsum, color = as.factor(Phenophase_Derived) )) +
    stat_ecdf() +
    scale_colour_viridis_d() +
    ggtitle("Cumulative proportion of each phenophase as heatsum increases")

ggplot(pd, aes(x = DoY, color = as.factor(Phenophase_Derived))) +
    stat_ecdf() +
    scale_colour_viridis_d() +
    ggtitle("Cumulative proportion of each phenophase as day of year increases")

# apply link function

logit <- function(x) log(x/(1-x))
lco <- logit(cum_pr_pd) #log cumulative odds
# logistic <- function(x) exp(x)/(1 + exp(x))
# logistic(cum_pr_pd)
# logistic(pr_pd)

############################################################
# 3. Identify Relevant Summary Statistics
############################################################
# h and k are relevant summary statistics.


############################################################
# 4. Build a Generative Model (Create fake data)
############################################################
writeLines(readLines("~/Documents/research_phenolology/phenology.stan"))
gen_model <- "~/Documents/research_phenolology/generative_model.stan"

generated_data_model <- stan(file=gen_model, iter=1,
            chains=1, algorithm="Fixed_param")

sim_data <- data.frame(heatsum = t(extract(generated_data_model)$heatsum),
           phenophase = t(extract(generated_data_model)$phenophase))

ggplot(sim_data, aes(x = heatsum, y=phenophase)) +
    geom_point() +
    ggtitle('Simulated heatsum and phenophase data')

############################################################
# 5. Analyze the Generative Ensemble (Build a model with your fake data and see if you can get the parameters back)
############################################################

R <- 100 # 1000 draws from the Bayesian joint distribution
N <- nrow(sim_data)
K <- length(unique(sim_data$phenophase))
heatsum <- sim_data$heatsum
phenophase <- sim_data$phenophase

stan_rdump(c("N", "K", "heatsum", "phenophase"), file="simulated_data_from_gen_model.data.R")




############################################################
# 5a. Analyze the Simulated Data
############################################################

sim_data_for_stan <- read_rdump("simulated_data_from_gen_model.data.R")

sim_fit <- stan(file='phenology.stan', data=sim_data_for_stan,
            iter=10000, warmup=1000, chains=5)

util$check_all_diagnostics(sim_fit)

params <- extract(sim_fit)
posterior_cp <- as.array(sim_fit) #extract posterior draws
lp_cp <- log_posterior(sim_fit)

np_cp <- nuts_params(sim_fit) #nuts params
head(np_cp)

color_scheme_set("darkgray")
mcmc_parcoord(posterior_cp, np = np_cp) # parallel coordinates plot. show one line per iteration, connecting the parameter values at this iteration, with divergences in red. let's you see global patterns in divergences

mcmc_pairs(posterior_cp, np = np_cp, pars = c("beta", "c[1]", "c[2]")) # show univariate histograms and bivariate scatter plots for selected parameters and is especially useful in identifying collinearity between variables (which manifests as narrow bivariate plots) as well as the presence of multiplicative non-identifiabilities (bananas). Each bivariate plot occurs twice and contains half the chains - so you can compare if chains produce similar results

scatter_c1_cp <- mcmc_scatter( # investigate a single bivariate plot to investigate it more closely. plot centered parameterization
    posterior_cp,
    pars = c("c[1]", "beta"),
    transformations = list(beta = "log"),
    np = np_cp
)

scatter_c1_cp

mcmc_trace(posterior_cp, np=np_cp) #trace plot. time series plot of markov chains. use window argument to zoom in on suspicious sections.

color_scheme_set("red")
mcmc_nuts_divergence(np_cp, lp_cp) # understand how divergences interact with the model globally. Identify light tails and incomplete exploration of target distribution. use chain argument to overlay the plot for a particular Markov chain on the plot

color_scheme_set("red")
mcmc_nuts_energy(np_cp)
            # energy plot. shows overlaid histograms of the (centered) marginal energy distribution piE and the first-differenced distribution pi_deltaE. id overly heavy tails (also challenging for sampling). the energy diagnostic for HMC (and the related energy-based Bayesian fraction of missing info) quantifies the heaviness of the tails of the posterior. Ideally the two histograms will look the same

# MCMC diagnostics: convergence?

# Rhat: potential scale reduction statistic
# compare a chain's behavior to other randomly intialized chains. Split R_hat measures ratio of the average variance of draws within each chain to the variance of the pooled draws across chains. If all chains at equilibrium, 1. If they haven't converged, > 1.

rhats <- rhat(sim_fit)
color_scheme_set("brightblue")
mcmc_rhat(rhats) +
    yaxis_text(hjust = 1)

# Effective sample size
# estimate of the number of independent draws from the posterior dist of the estimand of interest. n_eff in stan is based on ability of draws to estimate the true mean value of the param. because draws are not independent if there is autocorrelation between draws, neff is usually smaller than total N. the larger the ration of n_eff to N, the better. ratios depend not just on the model but on the algorithm used to fit the model

ratios_cp <- neff_ratio(sim_fit)
print(ratios_cp)

mcmc_neff(ratios_cp, size = 2) +
    yaxis_text(hjust = 1)

# Autocorrelation
#n_eff/N decreases as autocorrelation becomes more extreme. Visualize autocorrelation using mcmc_acf or mcmc_acf_bar. Postive autocorrelation is bad because it means the chain gets stuck. Ideally, it drops quickly to zero as lag increasses. negative autocorrelation indicates fast convergence of sample mean towards true

mcmc_acf(posterior_cp, lags = 10)
############

N <- nrow(pd)
