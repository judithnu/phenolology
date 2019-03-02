# Analysis for lab meeting - Prince George data

library(rstan)
library(bayesplot)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

util <- new.env()
source('~/Documents/classes/STANworkshop/material/1 - workflow/stan_utility.R', local=util)

# read in data ------------------------
phenologyheatsum_data <- read.csv('data/phenology_heatsum.csv')

phd <- phenologyheatsum_data

mphd <- subset(phd, Sex == "MALE") # separate out male and female datasets
fphd <- subset(phd, Sex == "FEMALE")

K = 3
V = 3
t = nrow(mphd)
w[t] = mphd$Phenophase_Derived
