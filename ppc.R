#### posterior predictive checks
## for forcing units across the range observed in my data (length X), simulate 30 phenological states for each model configuration (Y).
## This will result in X x 30 x Y states

forcingtype = 'scaled_ristos'



library(tidyverse)
library(rstan)
library(bayesplot)
library(rethinking)


#model
fmod <- readRDS("slopes_nc_scaled_ristos_FEMALE2019-08-27_climatena.rds") #%>%
    as.data.frame()

mmod <- readRDS("slopes_nc_scaled_ristos_MALE2019-08-27_climatena.rds") %>%
    as.data.frame()

# original data (this code should match relevant bits in run_stan)
phenology_data <- read.csv("data/phenology_heatsum.csv",
                           stringsAsFactors = FALSE, header = TRUE
) %>%
    filter(forcing_type == forcingtype)

## provenance
SPU_dat <- read.csv("../research_phd/data/OrchardInfo/LodgepoleSPUs.csv",
                    header=TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select(SPU_Name, Orchard)



# Data Processing ##################
# join provenance and phenology data

phendf <- phenology_data %>%
    na.omit() %>%
    left_join(SPU_dat) %>%
    unique() %>%
    mutate(Phenophase = as.factor(Phenophase_Derived)) %>%
    select(-Phenophase_Derived)


# identify combinations of effects that actually occur

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")


udf <- unique_grouper(fdf, mdf)
rm(fdf, mdf)

rordlogit()
dordlogit()

# sequence from 0 to x forcing units, seq length 30
# choose 20% of the unique groups
# choose 500 draws for each udf
# calculate phi for each draw for each udf
# simulate 30 states for each draw of each udf
 30*500*(2955*.1)*30

# maybe do more draws for fewer groups?

 #male
mdf <- df # from run_stan

*7027

sampling(fmod, )
