# model examination
library(rstan)
library(shinystan)
library(bayesplot)
#library(ggplot2)

# MODEL DATA #####################


ffit.stan <- readRDS("2019-10-28phenologyFEMALE.rds")
mfit.stan <- readRDS("2019-10-28phenologyMALE.rds")



fshiny <- as.shinystan(ffit.stan)
launch_shinystan(fshiny)

mshiny <- as.shinystan(mfit.stan)
launch_shinystan(mshiny)

# female model #########################################
fsum <- rstan::summary(ffit.stan)$summary
fsum <- as.data.frame(fsum)
farray <- as.array(ffit.stan) # diagnostics
#fpostdf <- as.data.frame(ffit.stan)
#fpost_re <- extract.samples(ffit.stan) # requ
#post <- as.array(ffit.stan)
flp <- log_posterior(ffit.stan)
#param_names <- attributes(postdf)$dimnames$parameters
#fparam_names <- colnames(fpostdf)
fnp <- nuts_params(ffit.stan) #nuts params
frhats <- rhat(ffit.stan)
fratios <- neff_ratio(ffit.stan)

### Male model ############

msum <- rstan::summary(mfit.stan)$summary
msum <- as.data.frame(msum)
marray <- as.array(mfit.stan) # diagnostics
#postlist <- rstan::extract(ffit.stan)
#mpostdf <- as.data.frame(mfit.stan)
#mpost_re <- extract.samples(mfit.stan)
#post <- as.array(ffit.stan)
mlp <- log_posterior(mfit.stan)
#param_names <- attributes(postdf)$dimnames$parameters
#mparam_names <- colnames(mpostdf)
mnp <- nuts_params(mfit.stan) #nuts params
mrhats <- rhat(mfit.stan)
mratios <- neff_ratio(mfit.stan)

# Diagnostics #######################
sex = "MALE"

if (sex=="FEMALE") {
    sarray <- farray
    np <- fnp
    # param_names <- fparam_names
    lp <- flp
    rhats <- frhats
    ratios <- fratios
    color_scheme_set("purple")
}

if (sex=="MALE") {
    sarray <- marray
    np <- mnp
    #param_names <- mparam_names
    lp <- mlp
    rhats <- mrhats
    ratios <- mratios
    color_scheme_set("yellow")
}

# Trace and rank plots
mcmc_trace(sarray, regex_pars = "site") + ggtitle(paste(sex, "site"))
mcmc_trace(sarray, regex_pars = "prov") + ggtitle(paste(sex, "prov"))
mcmc_trace(sarray, regex_pars = "sigma") + ggtitle(paste(sex, "sigma"))
mcmc_trace(sarray, regex_pars = "year") + ggtitle(paste(sex, "year"))
mcmc_trace(sarray, regex_pars = "kappa") + ggtitle(paste(sex, "kappa"))

## rank plots

mcmc_rank_hist(ffit.stan, regex_pars = "year")
mcmc_rank_overlay(ffit.stan, regex_pars = c("year", "prov"))
mcmc_rank_overlay(ffit.stan, regex_pars = "site")
mcmc_rank_overlay(ffit.stan, regex_pars = "clone\\[.0")
mcmc_rank_overlay(ffit.stan, regex_pars = "sigma")
mcmc_rank_overlay(ffit.stan, regex_pars = "mean")

# Divergences

divergences <- filter(np, Parameter=="divergent__" & Value==1)
print("divergences = ")
nrow(divergences)

#color_scheme_set("darkgray")
# mcmc_parcoord(sarray, np = np, pars = param_names[c(3:16, 276:294)]) # parallel coordinates plot. show one line per iteration, connecting the parameter values at this iteration, with divergences in red. let's you see global patterns in divergences
# mcmc_parcoord(sarray, np=np, regex_pars = c("clone"))

mcmc_pairs(sarray, np = np, regex_pars = c("kappa")) # show univariate histograms and bivariate scatter plots for selected parameters and is especially useful in identifying collinearity between variables (which manifests as narrow bivariate plots) as well as the presence of multiplicative non-identifiabilities (bananas). Each bivariate plot occurs twice and contains half the chains - so you can compare if chains produce similar results
mcmc_pairs(sarray, np = np, regex_pars = c("site"))
mcmc_pairs(sarray, np = np, regex_pars="prov")
mcmc_pairs(sarray, np = np, regex_pars = "sigma")
mcmc_pairs(sarray, np = np, pars = c("b_year[1]", "b_year[2]", "b_year[3]", "b_year[4]", "sigma_year"))
mcmc_pairs(sarray, np=np, pars=c("b_site[1]", "b_clone[1]", "b_prov[1]"))


mcmc_nuts_divergence(np, lp)
mcmc_nuts_divergence(np, lp, chain = 1) # understand how divergences interact with the model globally. Identify light tails and incomplete exploration of target distribution. use chain argument to overlay the plot for a particular Markov chain on the plot
mcmc_nuts_divergence(np, lp, chain = 2)
mcmc_nuts_divergence(np, lp, chain = 3)
mcmc_nuts_divergence(np, lp, chain = 4)
mcmc_nuts_divergence(np, lp, chain = 5)
mcmc_nuts_divergence(np, lp, chain = 6)
mcmc_nuts_divergence(np, lp, chain = 7)
mcmc_nuts_divergence(np, lp, chain = 8)

### ENERGY

mcmc_nuts_energy(np) + ggtitle(paste(sex, "energy plot"))
# energy plot. shows overlaid histograms of the marginal energy distribution piE and the first-differenced distribution pi_deltaE. id overly heavy tails (also challenging for sampling). the energy diagnostic for HMC (and the related energy-based Bayesian fraction of missing info) quantifies the heaviness of the tails of the posterior. Ideally the two histograms will look the same

#Look at the pairs plot to see which primitive parameters are correlated with the energy__ margin. There should be a negative relationship between lp__ and energy__ in the pairs plot, which is not a concern because lp__ is the logarithm of the posterior kernel rather than a primitive parameter.

energy <- dplyr::filter(np, Parameter== "energy__")


# Rhat: potential scale reduction statistic
# compare a chain's behavior to other randomly intialized chains. Split R_hat measures ratio of the average variance of draws within each chain to the variance of the pooled draws across chains. If all chains at equilibrium, 1. If they haven't converged, > 1.

print("rhats > 1 for")
names(which(rhats > 1.01))
which(is.na(rhats))

color_scheme_set("purple")
mcmc_rhat(rhats) +
    yaxis_text(hjust = 1)

# Effective sample size
# estimate of the number of independent draws from the posterior dist of the estimand of interest. n_eff in stan is based on ability of draws to estimate the true mean value of the param. because draws are not independent if there is autocorrelation between draws, neff is usually smaller than total N. the larger the ration of n_eff to N, the better. ratios depend not just on the model but on the algorithm used to fit the model

print("Effective sample size < 0.5 for")
names(which(ratios < 0.5))
badratios <- ratios[which(ratios<0.5)]

print("Effective sample size < 0.1 for")
names(which(ratios < 0.1))
vbadratios <- ratios[which(ratios < 0.1)]

mcmc_neff(badratios, size = 1.5) +
    yaxis_text(hjust = 1)

mcmc_neff_hist(neff_ratio(ffit.stan)) # THIS IS FEMALE SPECIFIC, CHANGE TO CORRECT SEX

# calculate bulk andtail ess


ness <- data.frame(parameter = attr(sarray, which="dimnames")[3], tail=NA, bulk=NA, rhat=NA)
for (i in 1:nrow(ness)) {
    ness$tail[i] <- ess_tail(sarray[,,i])
    ness$bulk[i] <- ess_bulk(sarray[,,i]) # rank normalized
    ness$rhat[i] <- Rhat(sarray[,,i])
}

ness <- ness[!str_detect(ness$parameters, "state_rep"),]

head(dplyr::arrange(ness, tail))
dplyr::arrange(ness, -tail)
dplyr::arrange(ness, bulk)
dplyr::arrange(ness, -bulk)
dplyr::arrange(ness, rhat) %>%
    filter(rhat > 1.01)

# Autocorrelation
#n_eff/N decreases as autocorrelation becomes more extreme. Visualize autocorrelation using mcmc_acf or mcmc_acf_bar. Postive autocorrelation is bad because it means the chain gets stuck. Ideally, it drops quickly to zero as lag increasses. negative autocorrelation indicates fast convergence of sample mean towards true

mcmc_acf(sarray, lags = 10, regex_pars = c("kappa"))
mcmc_acf(sarray, lags=10, regex_pars = c('b_site'))
mcmc_acf(sarray, lags=10, regex_pars=c('b_prov'))
mcmc_acf(sarray, lags=10, regex_pars=c('b_year'))

# Parameter Estimation ##################

#compare m and f
#modify for plotting
fsum <- fsum[rownames(fsum) %in% rownames(msum),]
msum <- msum[rownames(msum) %in% rownames(fsum),]
fsum <- fsum[!rownames(fsum) %in% c("lp__", "kappa[1]", "kappa[2]", "beta"),]
msum <- msum[!rownames(msum) %in% c("lp__", "kappa[1]", "kappa[2]", "beta"),]
plot(fsum[,1], msum[,1], pch=rownames(fsum))
abline(0,1)

mvf <- lm(fsum$mean ~ msum$mean)

# central posterior uncertainty intervals. by default shows 50% intervals (thick) and 90% intervals (thinner outer). Use point_est to show or hide point estimates
color_scheme_set("red")
fint_bsp <- mcmc_intervals(farray, regex_pars = c("b_site", "b_prov"))
mint_bsp <- mcmc_intervals(marray, regex_pars = c("b_site", "b_prov"))

compare_fm(fint_bsp, mint_bsp)

farea_k <- mcmc_areas(farray, regex_pars = c("kappa"))
marea_k <- mcmc_areas(marray, regex_pars = c("kappa"))

compare_fm(farea_k, marea_k)

fint_clone <- mcmc_intervals(farray, regex_pars=c("clone")) + ggtitle("female")
mint_clone <- mcmc_intervals(marray, regex_pars=c("clone")) + ggtitle("male")

#NB
fint_beta <- mcmc_intervals(farray, regex_pars = c("beta"))
mint_beta <- mcmc_intervals(marray, regex_pars = c("beta"))

compare_fm(fint_beta, mint_beta)
mcmc_intervals(fpostdf, regex_pars = c("prov", "site", "sigma")) + ggtitle("Scaled ristos effects")
mcmc_intervals(fpostdf, regex_pars = c("clone")) + ggtitle("Clone intercepts")




