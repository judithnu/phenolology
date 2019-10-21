# Average Predictive comparisons

source('phenology_functions.R')

library(tidyverse)
library(gtools)

fmod <- readRDS("slopes_nc_scaled_ristos_FEMALE2019-10-04climatena.rds") %>%
    as.data.frame() %>%
  select(-contains("state_rep")) %>%
  sample_frac(0.2)

mmod <- readRDS("slopes_nc_scaled_ristos_MALE2019-10-04climatena.rds") %>%
    as.data.frame() %>%
  select(-contains("state_rep")) %>%
  sample_frac(0.2)


# original data (this code should match relevant bits in run_stan)
phendf <- read_data(slim = FALSE)

# identify combinations of effects that actually occur

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")


udf <- unique_grouper(fdf, mdf)

fpardf <- build_par_df(mcmcdf = fmod, datdf = fdf, sex = "FEMALE")
mpardf <- build_par_df(mcmcdf = mmod, datdf = mdf, sex = "MALE")

pardf <- rbind(fpardf, mpardf)


# calculate the start of forcing for all model configurations for the first "observation"
udf$SiteID[1]


site <- select(fmod, matches(paste("b_site", fdf$SiteID[1], sep="\\[")))
prov <- select(fmod, matches(paste("b_prov", fdf$ProvenanceID[1], sep="\\[")))
clone <- select(fmod, matches(paste("b_clone", fdf$CloneID[1], sep="\\[")))
beta <- select(fmod, beta)

betas <- beta + site + prov + clone
kappa1 <- select(fmod, `kappa[1]`)
kappa2 <- select(fmod, `kappa[2]`)

calcstageforcing <- function(p=0.2, beta=betas, kappa) {
  prob <- (logit(p) + kappa)/beta
  return(prob)
}

start <- calcstageforcing(kappa=kappa1)
