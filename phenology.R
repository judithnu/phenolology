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
# mutate(Phenophase_Stages = case_when(str_detect(Phenophase, "receptive20") ~ 2,
#                                      str_detect(Phenophase, "pollenshed20") ~ 2,
#                                      str_detect(Phenophase, "receptive80") ~ 3,
#                                      str_detect(Phenophase, "pollenshed80") ~ 3)) #change phenophase to 2 for 20% active and 3 for 80% done
## Get a list of dates trees were definitely observed at each Orchard on each Site in each Year

## Add stage 1 using collection dates
### STOPPED HERE BROKEN BELOW. MERGE SEEMS TO NOT ADD NA PHENOPHASES, WHICH I THINK IT SHOULD


# cw_temp$SiteOrchardYear_Index <- group_indices(cw_temp, Site, Orchard, Year)  # create index of site, orchard, and year
# observation_dates <- cw_temp %>%
#     select(SiteOrchardYear_Index, Date) %>%
#     rename(ObservationDate = Date)
# unique()
#
# temp <- full_join(observation_dates, cw_temp)
# # So now I need a list of dates that each ramet in a given year+site+orchard+sex has a record at
# cw_temp$SiteOrchardYearSexCloneXY_Index <- group_indices(cw_temp, Site, Orchard, Year, Clone, X, Y)
# phenophase_dates <- cw_temp %>%
#     select(SiteOrchardYearSexCloneXY_Index, Date) %>%
#     unique()
