library(rethinking)
library(dplyr)
library(bayesplot)
#library(magrittr)

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


#drop "non-numeric" cols to make stan quit complaining
phendf <- select(phendf, -Site, -SPU_Name)

mdf <- subset(phendf, Sex == "MALE") %>%
    select(-Sex)

fdf <- subset(phendf, Sex == "FEMALE")

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

#varying intercepts on clones and Provenance, varying provenance slopes
mSpIcp <- ulam(
    alist(
        #likelihood
        Phenophase_Derived ~ ordered_logistic(phi, cutpoints),
        #model
        phi <- b_provenance[ProvenanceID]*Heatsum + a_clone[CloneID] + a_provenance[ProvenanceID],
        #priors
        b_provenance[ProvenanceID] ~ dbeta(.5,5),
        cutpoints ~ dnorm(197, 176),
        a_clone[CloneID] ~ dnorm(0, sigma_clone),
        a_provenance[ProvenanceID] ~ dnorm(0, sigma_provenance),
        #hyperpriors
        sigma_clone ~ dcauchy(0,1),
        sigma_provenance ~ dcauchy(0,1)
    ),
    data = mdf,
    start= list(a1=100*.2, a2=250*.2),
    warmup = 2000, iter=10000, chains = 10, cores = 10
)

m2 <- mSpIcp
precis(m2, depth=2)
plot(m2)
traceplot(m2)

sfm2 <- m2@stanfit
saveRDS(sfm2, file="2019-03-03_mSpIcp.rds")
stancode(sfm2)
write(stancode(sfm2), file="2019-03-03_mSpIcp.stan")
