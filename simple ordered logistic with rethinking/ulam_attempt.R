library(rethinking)
library(dplyr)
library(bayesplot)

util <- new.env()
source('~/Documents/classes/STANworkshop/material/1 - workflow/stan_utility.R', local=util)

options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

phenology_data <- read.csv("~/phenolology/data/stan_input/phenology_heatsum.csv", stringsAsFactors = FALSE, header = TRUE)
SPU_dat <- read.csv("~/phd/data/OrchardInfo/LodgepoleSPUs.csv", header=TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select(SPU_Name, Orchard)

phendf <- phenology_data
phendf <- dplyr::left_join(phenology_data, SPU_dat)
#Create Clone IDs
phendf$CloneID <- group_indices(phendf, Clone)
#Create OrchardIDs
phendf$OrchardID <- group_indices(phendf, Orchard)
#Create ProvenanceIDs
phendf$ProvenanceID <- group_indices(phendf, SPU_Name)
#Create SiteIDs
phendf$SiteID <- group_indices(phendf, Site)
#Create SexIDs
phendf$SexID <- group_indices(phendf, Sex)


#drop "non-numeric" cols to make stan quit complaining
phendf <- select(phendf, -Site, -SPU_Name, -Sex)

# mdf <- subset(phendf, Sex == "MALE") %>%
#     select(-Sex)
#
# fdf <- subset(phendf, Sex == "FEMALE")

# read in data for priors
heatsum_priors_dat <- read.csv('~/phd/data/PhenologyAndPollenCounts/orchard_heatsums_WalshWebber2008.csv')
# calculate priors

simple_beta_prior = 0.5
pre_pollination_summary <- heatsum_priors_dat %>%
    dplyr::filter(period == 'pre-pollination') %>%
    dplyr::summarise(mean=mean(Tsum_day, na.rm=TRUE), sd=sd(Tsum_day, na.rm=TRUE))
pre_pollination_summary <- pre_pollination_summary*simple_beta_prior

cutpoint_mean <- pre_pollination_summary$mean
cutpoint_sd <- pre_pollination_summary$sd

#varying intercepts on clones and Provenance
mIcp <- ulam(
    alist(
        #likelihood
        Phenophase_Derived ~ ordered_logistic(phi, cutpoints),
        #model
        phi <- beta*Heatsum + a_provenance[ProvenanceID] + a_clone[CloneID],
        #priors
        beta ~ dbeta(.5,5),
        cutpoints ~ dnorm(197, 176),
        a_clone[CloneID] ~ dnorm(0, sigma_clone),
        a_provenance[ProvenanceID] ~ dnorm(0, sigma_provenance),
        #hyperpriors
        sigma_clone ~ exponential(1.5),
        sigma_provenance ~ exponential(1.5)
    ),
    data = mdf,
    start= list(a1=100*.2, a2=250*.2),
    warmup = 1000, iter=5000, chains = 10, cores = 10
)

mSsIcp <- ulam(
  alist(
    #likelihood
    Phenophase_Derived ~ ordered_logistic(phi, cutpoints),
    #model
    phi <- beta[SexID]*Heatsum + a_provenance[ProvenanceID] + a_clone[CloneID],
    #priors
    beta[SexID] ~ dbeta(.5,5),
    cutpoints ~ dnorm(197, 176),
    a_clone[CloneID] ~ dnorm(0, 1),
    a_provenance[ProvenanceID] ~ dnorm(0, 1)
    #hyperpriors
    #sigma_clone ~ exponential(1.5),
    #sigma_provenance ~ exponential(1.5)
  ),
  data = phendf,
  start= list(a1=100*.2, a2=250*.2),
  warmup = 1000, iter=5000, chains = 10, cores = 10
)

fit <- mSsIcp@stanfit
precis(fit, depth=2)
plot(fit)
plot(m3)
traceplot(m3)

saveRDS(fit@stanfit, file = "~/phenolology/2019-03-04_mSsIcp.rds")

fit <- readRDS("~/Documents/research_phenolology/2019-03-04_mSsIcp.rds")


saveRDS(sfm3, file="2019-03-03_mSpIcp.rds")
stancode(sfm3)
write(stancode(sfm3), file="2019-03-03_mIcp.stan")
m3fit <- readRDS("~/Documents/research_phenolology/2019-03-03_mIcp.rds")

util$check_all_diagnostics(m2fit)

# VISUAL MCMC DIAGNOSTICS
############################################################
#params <- extract(fit)
posterior_cp <- as.array(fit) #extract posterior draws
lp_cp <- log_posterior(fit)

param_names <- attributes(posterior_cp)$dimnames$parameters

np_cp <- nuts_params(fit) #nuts params

head(np_cp)

color_scheme_set("darkgray")
mcmc_parcoord(posterior_cp, np = np_cp, pars = param_names[c(1:4, 115:120)]) # parallel coordinates plot. show one line per iteration, connecting the parameter values at this iteration, with divergences in red. let's you see global patterns in divergences
mcmc_parcoord(posterior_cp, np = np_cp, pars = param_names[c(7:114, 119:120)])


mcmc_pairs(posterior_cp, np = np_cp, pars = c("b_provenance[1]", "cutpoints[1]", "a_clone[1]", "a_provenance[1]")) # show univariate histograms and bivariate scatter plots for selected parameters and is especially useful in identifying collinearity between variables (which manifests as narrow bivariate plots) as well as the presence of multiplicative non-identifiabilities (bananas). Each bivariate plot occurs twice and contains half the chains - so you can compare if chains produce similar results

scatter_c1_cp <- mcmc_scatter( # investigate a single bivariate plot to investigate it more closely. plot centered parameterization
    posterior_cp,
    pars = c("c[1]", "beta"),
    transformations = list(beta = "log"),
    np = np_cp
)

scatter_c1_cp

mcmc_trace(posterior_cp, np=np_cp) #trace plot. time series plot of markov chains. use window argument to zoom in on suspicious sections.

color_scheme_set("red")
mcmc_nuts_divergence(np_cp, lp_cp, chain = 2) # understand how divergences interact with the model globally. Identify light tails and incomplete exploration of target distribution. use chain argument to overlay the plot for a particular Markov chain on the plot

color_scheme_set("red")
mcmc_nuts_energy(np_cp)
# energy plot. shows overlaid histograms of the (centered) marginal energy distribution piE and the first-differenced distribution pi_deltaE. id overly heavy tails (also challenging for sampling). the energy diagnostic for HMC (and the related energy-based Bayesian fraction of missing info) quantifies the heaviness of the tails of the posterior. Ideally the two histograms will look the same

# MCMC diagnostics: convergence?

# Rhat: potential scale reduction statistic
# compare a chain's behavior to other randomly intialized chains. Split R_hat measures ratio of the average variance of draws within each chain to the variance of the pooled draws across chains. If all chains at equilibrium, 1. If they haven't converged, > 1.

rhats <- rhat(fit)
color_scheme_set("brightblue")
mcmc_rhat(rhats) +
    yaxis_text(hjust = 1)

# Effective sample size
# estimate of the number of independent draws from the posterior dist of the estimand of interest. n_eff in stan is based on ability of draws to estimate the true mean value of the param. because draws are not independent if there is autocorrelation between draws, neff is usually smaller than total N. the larger the ration of n_eff to N, the better. ratios depend not just on the model but on the algorithm used to fit the model

ratios_cp <- neff_ratio(m2fit)
print(ratios_cp)

mcmc_neff(ratios_cp, size = 2) +
    yaxis_text(hjust = 1)

# Autocorrelation
#n_eff/N decreases as autocorrelation becomes more extreme. Visualize autocorrelation using mcmc_acf or mcmc_acf_bar. Postive autocorrelation is bad because it means the chain gets stuck. Ideally, it drops quickly to zero as lag increasses. negative autocorrelation indicates fast convergence of sample mean towards true

mcmc_acf(posterior_cp, lags = 10)


# PARAMETER ESTIMATES
############################################################
# Plotting MCMC draws

posterior <- as.array(m3)
dim(sim_fit)

dimnames(posterior)

### Posterior uncertainty intervals ("credible intervals")

# central posterior uncertainty intervals. by default shows 50% intervals (thick) and 90% intervals (thinner outer). Use point_est to show or hide point estimates
color_scheme_set("red")
mcmc_intervals(posterior, pars = param_names[c(1:4)], point_est = "none")
mcmc_intervals(posterior, regex_pars = "a_clone")
mcmc_intervals(posterior, regex_pars = "cutpoints")
mcmc_intervals(posterior, regex_pars = "a_provenance", point_est = "none")

#NB
mcmc_intervals(posterior_cp, regex_pars = c("beta")) + ggtitle("Female and male transition speed")
mcmc_intervals(posterior_cp, regex_pars = c("provenance")) + ggtitle("Provenance intercepts")
mcmc_intervals(posterior_cp, regex_pars = c("clone")) + ggtitle("Clone intercepts")

# show uncertainty intervals as shaded areas under estimated posterior density curves

mcmc_areas(
    posterior,
    pars = c("c[1]", "c[2]", "beta"),
    prob = 0.8, # 80% intervals
    prob_outer = 0.99, # 99% intervals
    point_est = "mean"
)

#PRES
mcmc_areas(posterior_cp,
           regex_pars = "beta")

mcmc_areas_ridges(posterior_cp, regex_pars=c("clone"))

### Univariate marginal posterior distributions
# Look at histograms or kernel density estimates of marginal posterior distributions, chains together or separate

# plot marginal posterior distributions combining all chains. Tranformations argument can be helpful
color_scheme_set("green")
mcmc_hist(posterior, regex_pars="a_provenance")

color_scheme_set("brightblue")
mcmc_hist_by_chain(posterior, pars = c("c[1]", "c[2]", "beta"))

# plot densities rather than histograms
color_scheme_set("purple")
mcmc_dens(posterior, regex_pars = c("a_provenance"))
# by chain
mcmc_dens_overlay(posterior, pars = c("c[1]", "c[2]", "beta"))

# plot violins

color_scheme_set("teal")
mcmc_violin(posterior, regex_pars = "a_provenance", probs = c(0.1, 0.5, 0.9))

### Bivariate plots

# just a scatterplot of two params
color_scheme_set("gray")
mcmc_scatter(posterior, pars = c("c[1]", "c[2]"))
mcmc_scatter(posterior, pars = c("c[1]", "beta"))
mcmc_scatter(posterior, pars = c("c[2]", "beta"))

# hexagonal binner helps avoid overplotting
mcmc_hex(posterior, pars = c("c[1]", "c[2]"))
mcmc_hex(posterior, pars = c("c[1]", "beta"))
mcmc_hex(posterior, pars = c("c[2]", "beta"))

# with more than two params. use hex plots instead with off_diag_fun = "hex". half markov chains above and half below diagonal
color_scheme_set("pink")
mcmc_pairs(posterior, pars = c("c[1]", "c[2]", "beta"), off_diag_args = list(size = 1.5))

### Trace plots
# no divergence info (see mcmc diagnostics graphs)

color_scheme_set("blue")
mcmc_trace(posterior, pars = c("c[1]", "c[2]", "beta"))

# with mixed colour scheme

color_scheme_set("mix-blue-red")
mcmc_trace(posterior, pars = c("c[1]", "c[2]", "beta"),
           facet_args = list(ncol=1, strip.position = "left")) #facet_arg passes parameters to facet_wrap in ggplot2. ncol puts traceplots in a column instead of side by side and strip.position moves the facet labels

# with viridis colour scheme
color_scheme_set("viridis")
mcmc_trace(posterior, pars = c("c[1]", "c[2]", "beta"))

# use points instead of lines and highlight a single chain
mcmc_trace_highlight(posterior, regex_pars = "a_provenance", highlight = 3)



# ------------------
