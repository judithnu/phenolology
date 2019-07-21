# Run phenology model

# Set sex and forcing type ################
# Choose sex and forcing type
sex <- "FEMALE"
#sex <- "MALE"

forcingtype <- "scaled_ristos"


# Dependencies and options ##################
# library(rethinking)
library(tidyverse)
library(rstan)

# rstan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# Functions #################

# Stan can only take consecutive integers for factors, so turn factors into consecutive integers
stanindexer <- function(df) {
  df$CloneID <- group_indices(df, Clone)
  df$OrchardID <- group_indices(df, Orchard)
  df$ProvenanceID <- group_indices(df, SPU_Name)
  df$SiteID <- group_indices(df, Site)
  df$YearID <- group_indices(df, Year)
  df$TreeID <- group_indices(df, TreeUnique)
  return(df)
}

# Drop non-integer columns. Necessary if using rethinking to draft new stan code
# stancleaner <- function(df) {
#     df <- dplyr::select(df, -Site, -SPU_Name, -Sex, -TreeID, -Phenophase, -Date, -Clone, -forcing_type)
#     return(df)
# }

# write file that stan can actually use
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

  rstan::stan_rdump(c("N", "K", "Nsite","Nprovenance", "Nclone", "Nyear", "SiteID", "ProvenanceID", "CloneID", "YearID", "forcing", "state"), file)
}


# Read in data ##################
## phenology
phenology_data <- read.csv("data/phenology_heatsum.csv",
                           stringsAsFactors = FALSE, header = TRUE
) %>%
  filter(forcing_type == forcingtype)

if(forcingtype == "gdd") { #scale growing degree days
  phenology_data$sum_forcing <- phenology_data$sum_forcing/10
}

## sample for less domination by PGTIS

foo <- phenology_data %>%
  filter(Site=="PGTIS" & Year %in% c(2006, 2007, 2010, 2011, 2005)) 
           
uclone <- unique(foo$Clone) %>%
  base::sample(size=15)

pgtis_sampled <- foo %>%
  filter(Clone %in% uclone)

nopgtis <- filter(phenology_data, Site != "PGTIS")

phenology_data <- full_join(pgtis_sampled, nopgtis)

## provenance
SPU_dat <- read.csv("../phd/data/OrchardInfo/LodgepoleSPUs.csv",
                    header=TRUE, stringsAsFactors = FALSE) %>%
  dplyr::select(SPU_Name, Orchard)



# Data Processing ##################
# join provenance and phenology data

phendf <- phenology_data %>%
  na.omit() %>%
  left_join(SPU_dat) %>%
  unique()

# filter for sex of interest
df <- filter(phendf, Sex == sex)
df <- stanindexer(df)

# write out and read in data structured for stan
prepforstan(df, paste(sex, ".rdump", sep=""))

rdump <- read_rdump(paste(sex, ".rdump", sep=""))

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


# Fit model  #############
test <- stan("slopes.stan",
             model_name = paste(sex, "slopes gdd"),
             data = rdump,
             chains = 1, cores = 1, warmup = 20, iter = 25
) # quick check for small problems

fit <- stan("slopes.stan",
            model_name = paste(sex, "slopes gdd"),
            data = rdump,
            chains = 10, cores = 10, warmup = 1000, iter = 1300,
            control = list(max_treedepth = 15, adapt_delta = .9)
)

saveRDS(fit, file = paste(sex, "slopes_gdd.rds", sep=''))

