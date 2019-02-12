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
# 1. Conceptual Analysis -------------------------
############################################################
# I want to be able to predict when trees will flower and for how long. Flowering is controlled by heatsum, but may also have a genetic component. The relationship between heatsum and flowering is logistic. In, a logistic model of transition between phenophases with logit(-k(xi-h)) where xi is the heatsum for a given record, h describes the heatsum at which there is a 50% chance an individual tree has begun transitioned between phases (or the point at which 50% of the trees in a population have transitioned between phases). k describes how quickly trees transition from one to another. Trees can only go from 1 to 2 to 3 - they cannot go backwards. An ordered logistic distribution keeps transitions in order.

# Model is usually written as logit(ak - b*xi) where a is the "intercept" or cutpoint for each transition. to convert it to h, divide b by a.

############################################################
# 2. Define Observations ----------------------------
############################################################
phenology_data <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/data_formatted_and_derived//derived_phenophase.csv', header = TRUE, stringsAsFactors = FALSE)

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

############################################################
# 3. Identify Relevant Summary Statistics
############################################################
# h and k are relevant summary statistics.


############################################################
# 4. Build a Generative Model
############################################################
writeLines(readLines("~/Documents/research_phenolology/phenology.stan"))
gen_model <- "~/Documents/research_phenolology/phenology.stan"
############################################################
# 5. Analyze the Generative Ensemble
############################################################
R <- 100 # 1000 draws from the Bayesian joint distribution
N <- 3
#theta <- c(1,1,1)/sum(c(1,1,1))

############################################################
# 5a. Analyze the Prior Predictive Distribution
############################################################
simu_data <- list("N" = N, "theta" = theta)

fit <- stan(file=gen_model, data=simu_data,
            iter=R, warmup=0, chains=1, refresh=R,
            seed=4838282, algorithm="Fixed_param")

simu_ys <- extract(fit)$y
hist(simu_ys)

############

N <- nrow(pd)
