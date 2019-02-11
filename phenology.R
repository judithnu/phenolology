############################################################
# Initial setup
############################################################
library(ggplot2)
library(tidyr)
library(dplyr)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

util <- new.env()
source('stan_utility.R', local=util)

############################################################
# model
############################################################

phenology_data <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/data_formatted/derived_phenophase.csv', header = TRUE, stringsAsFactors = FALSE)

pd <- phenology_data

ggplot(pd, aes(x=DoY, y=Phenophase_Derived, colour=Sex)) +
    geom_jitter(shape = 3, alpha = 0.5) +
    facet_wrap(.~Year)

hist(pd$Phenophase_Derived, xlab = "Phenophase")
# discrete proportion of each phenophase value

pd_narm <- pd[!is.na(pd$Phenophase_Derived),] #drop nas
pr_pd <- table(pd_narm$Phenophase_Derived)/nrow(pd_narm)

# cumsum converts to cumulative proportions
cum_pr_pd <- cumsum(pr_pd)

# plot
plot(c(0:3), cum_pr_pd, type = "b", xlab = "Phenophase", ylab = "cumulative proportion", ylim=c(0,1))

# apply link function

logit <- function(x) log(x/(1-x))
lco <- logit(cum_pr_pd) #log cumulative odds
# logistic <- function(x) exp(x)/(1 + exp(x))
# logistic(cum_pr_pd)
# logistic(pr_pd)
