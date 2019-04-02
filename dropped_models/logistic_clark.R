library(rethinking)
library(dplyr)
library(bayesplot)

util <- new.env()
source('~/Documents/classes/STANworkshop/material/1 - workflow/stan_utility.R', local=util)

options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

phenology_data <- read.csv("data/stan_input/phenology_heatsum.csv", stringsAsFactors = FALSE, header = TRUE)
SPU_dat <- read.csv("~/Documents/research_phd/data/OrchardInfo/LodgepoleSPUs.csv", header=TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select(SPU_Name, Orchard)

phendf <- phenology_data %>%
    na.omit()
phendf <- dplyr::left_join(phenology_data, SPU_dat)

phendf$SexID <- group_indices(phendf, Sex)

phendf <- select(phendf, -Site, -SPU_Name, -Sex, -TreeID, -Phenophase, -Date, -Clone) %>%
    mutate(Phenophase_scaled = case_when(Phenophase_Derived == 3 ~ 1,
                                         Phenophase_Derived == 2 ~ 0.5,
                                         Phenophase_Derived == 1 ~ 0))

fit <- map2stan(
    alist(
        #likelihood
        Phenophase_Derived ~ dmultinom(size=1, prob = theta),
        theta = c()
        pp_continuous = 3/(1 + exp((-k[SexID]*(Heatsum-h[SexID])))),
        k[SexID] ~ dbeta(.5,5),
        h[SexID] ~ dnorm(190,100)),
    data = phendf,
    warmup = 5, iter=10, chains = 1, cores = 1
)

fit <- fit <- ulam(
    alist(
        #likelihood
        Phenophase_Derived ~ categorical(theta),
        theta ~ beta(1,2)),

    data = phendf,
    warmup = 5, iter=10, chains = 1, cores = 1
)

