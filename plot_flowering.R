library(rstan)
library(dplyr)

# functions ############
build_par_df <- function(mcmcdf, datdf = udf, sex) {

    udf <- filter(datdf, Sex==sex)
    mcmcdf$var <- rownames(mcmcdf))

    singledimpars <- mcmcdf %>%
        dplyr::select(quart, beta, sigma_site, sigma_prov, sigma_clone, sigma_year, contains("kappa"), contains("mean")) %>%
        rename(kappa1 = `kappa[1]`) %>%
        rename(kappa2 = `kappa[2]`)
    siteb <- tidypar(mcmcdf, "b_site", "SiteID")
    provb <- tidypar(mcmcdf, "b_prov", "ProvenanceID")
    cloneb <- tidypar(mcmcdf, "b_clone", "CloneID")
    yearb <- tidypar(mcmcdf, "b_year", "YearID")

    sitemerge <- left_join(udf, siteb) %>%
        select(-name)
    provmerge <- left_join(udf, provb) %>%
        select(-name)
    clonemerge <- left_join(udf, cloneb) %>%
        select(-name)
    yearmerge <- left_join(udf, yearb) %>%
        select(-name)

    pardf <- left_join(udf, clonemerge) %>%
        left_join(provmerge) %>%
        left_join(sitemerge) %>%
        left_join(yearmerge) %>%
        left_join(singledimpars)
}

# data ##############

fmod <- readRDS("FEMALE_slopes_scaled.rds")

mmod <- readRDS("MALE_slopes_scaled.rds")

fsum <- rstan::summary(fmod)$summary %>%
    data.frame() %>%
    dplyr::select(-n_eff, -Rhat, -se_mean, sd)
fsum <- t(fsum)
msum <- rstan::summary(mmod)$summary %>%
    data.frame() %>%
    dplyr::select(-n_eff, -Rhat, -se_mean, sd)
msum <- t(msum)


phenology_data <- read.csv("data/phenology_heatsum.csv",
                           stringsAsFactors = FALSE, header = TRUE
) %>%
    filter(forcing_type == forcingtype)

## provenance
SPU_dat <- read.csv("../research_phd/data/OrchardInfo/LodgepoleSPUs.csv",
                    header=TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select(SPU_Name, Orchard)

phendf <- phenology_data %>%
    na.omit() %>%
    left_join(SPU_dat) %>%
    unique()


# identify combinations of effects that actually occur

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")


udf <- unique_grouper(fdf, mdf)

fpardf <- build_par_df(mcmcdf = fsum, datdf = udf, sex = "FEMALE")
mpardf <- build_par_df(mcmcdf = mmod, datdf = udf, sex = "MALE")
pardf <- rbind(fpardf, mpardf) %>%
    group_by(IndSexGroup)