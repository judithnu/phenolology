############################################################
# Initial setup
############################################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

util <- new.env()
source('stan_utility.R', local=util)

############################################################
# model
############################################################

phenology_data <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/Combined_Wagner_Walsh_pollination_phenology.csv')

# Transform data so Phenophase is in 3 stages appropriate for model
# Model formatting

# Identify recorded start and end dates


