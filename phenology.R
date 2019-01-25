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

phenology_data <- read.csv('~/Documents/research_phd/analysis/censorship/PGTIScensorship.csv')
