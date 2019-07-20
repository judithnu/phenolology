#why don't effects cause differences

library(bayesplot)

# fmod

ggplot(fmod, aes(x=`b_prov[1]`, color="prov1"))+
    geom_density() +
    geom_density(aes(x=`b_prov[2]`, color="prov2")) +
    ggtitle("fmod")

ggplot(fpardf, aes(x=b_prov, color=SPU_Name)) +
    geom_density() +
    ggtitle("fpardf")


#investigate build_par_df

mcmcdf <- fmod
datdf <- udf
sex="FEMALE"


build_par_df <- function(mcmcdf, datdf = udf, sex) {

    udf <- filter(datdf, Sex==sex)
    mcmcdf$iter <- 1:nrow(mcmcdf)

    singledimpars <- mcmcdf %>%
        dplyr::select(iter, beta, sigma_site, sigma_prov, sigma_clone, sigma_year, contains("kappa"), contains("mean")) %>%
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
ggplot(pardf_trans, aes(x=fstart, color="fstart")) +
    geom_density() +
    geom_density(aes(x=fstart1_pop, color="fstart_pop"))

ggplot(pardf_trans, aes(x=b_prov, color=SPU_Name))+
    geom_density() +
    facet_grid(ProvenanceID ~ Sex)

mcmc_areas(fmod, regex_pars = "b_prov")
mcmc_areas(mmod, regex_pars = "b_prov")



hist(pardf_trans$b_site)
