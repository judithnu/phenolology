# Run phenology model

# Set sex and forcing type ################
# Choose sex and forcing type
#sex <- "FEMALE"
sex <- "MALE"

forcingtype = "scaled_ristos"

test = FALSE


# Dependencies and options ##################
# library(rethinking)
library(tidyverse)
library(rstan)

source('phenology_functions.R')

# rstan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# functions ################

# write data file that stan can actually use
prepforstan <- function(df, file) {
  N <- nrow(df)
  K <- length(unique(df$Phenophase_Derived))
  Nsite <- length(unique(df$SiteID))
  Nprovenance <- length(unique(df$ProvenanceID))
  Norchard <- length(unique(df$OrchardID))
  Nclone <- length(unique(df$CloneID))
  Ntree <- length(unique(df$TreeID))
  Nyear <- length(unique(df$YearID))
  
  SiteID <- df$SiteID
  ProvenanceID <- df$ProvenanceID
  OrchardID <- df$OrchardID
  CloneID <- df$CloneID
  TreeID <- df$TreeID
  YearID <- df$YearID
  
  forcing <- df$sum_forcing
  state <- df$Phenophase_Derived
  
  rstan::stan_rdump(c("N", "K", "Nsite","Nprovenance", "Nclone", "Nyear", "Ntree", "SiteID", "ProvenanceID", "CloneID", "YearID", "TreeID", "forcing", "state"), file)
}

# Read in and process data
phendf <- read_data(slim=FALSE)
#
# msf <- mean(phendf$sum_forcing)
# sdf <- sd(phendf$sum_forcing)
#
# phendf$normalized <- (phendf$sum_forcing - msf)/sdf
#
# phendf <- rename(phendf, sum_forcing_old = sum_forcing, sum_forcing=normalized)

# filter for sex of interest
df <- filter(phendf, Sex == sex)
df <- stanindexer(df)

# write out and read in data structured for stan
prepforstan(df, paste("data/stan_input/", sex, ".rdump", sep=""))

rdump <- read_rdump(paste("data/stan_input/", sex, ".rdump", sep=""))

# Draft stan code using rethinking ####
# slopedraft <- ulam(
#     alist(
#         Phenophase_Derived ~ ordered_logistic(phi, kappa),
#         phi <- sum_forcing * (beta + b_site[SiteID] + b_prov[ProvenanceID] + b_clone[CloneID] + b_year[YearID] + b_orch[OrchardID]),
#         #fixed priors
#         kappa ~ normal(7.5,2), #make gamma in stan
#         beta ~ exponential(3),
#         #adaptive priors
#         b_site[SiteID] ~ normal(0, sigma_site),
#         b_prov[ProvenanceID] ~ normal(0, sigma_prov),
#         b_clone[CloneID] ~ normal(0, sigma_clone),
#         b_year[YearID] ~ normal(0, sigma_year),
#         b_orch[OrchardID] ~ normal(0, sigma_orch),
#         #hyperpriors
#         sigma_site ~ exponential(2),
#         sigma_prov ~ exponential(2),
#         sigma_clone ~ exponential(2),
#         sigma_year ~ exponential(2),
#         sigma_orch ~ exponential(2)
#     ),
#     data=fdf, sample=FALSE, declare_all_data = FALSE)
#
# # write(stancode(fit_draft), file="slopes.stan")

################

# Fit model  #############
if (test==TRUE)
{
  test <- stan("phenology_deeper_noncentered.stan",
               model_name = paste("test", Sys.Date(), sex, forcingtype, sep="_"),
               data = rdump,
               pars = c( "z_year", "phi"), include=FALSE,
               chains = 1, cores = 1, warmup = 20, iter = 25,
               save_dso=FALSE
  ) # quick check for small problems
} else {
  fit <- stan("phenology_deeper_noncentered.stan",
              model_name = paste(Sys.Date(), sex, "phenology", forcingtype, sep="_"),
              data = rdump,
              pars = c( "phi"), include=FALSE,
              chains = 4, cores = 4, warmup = 1000, iter = 4500,
              control = list(max_treedepth = 12)#, adapt_delta = .9),
  )
  
  saveRDS(fit, file = paste(Sys.Date(), "phenology", sex, ".rds", sep=''))
}

