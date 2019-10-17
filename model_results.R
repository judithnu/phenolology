#### model results

# Create datasets holding parameter values - fstart, fend, and f50s

library(dplyr)
library(tidyr)
library(stringr)
library(rstan)
library(gtools)

source('phenology_functions.R')

# FUNCTIONS ####################################


#take a stan model dataframe and create a tidy dataframe for a given parameter. takes a dataframe, a string with the parameter name, and a string describing the param ID (e.g. par = "Site" and id="SiteID")
tidypar <- function(stanframe, param, id) {
    #wide to long format for parameters
    par <- stanframe %>%
        dplyr::select(contains(param), iter) %>%
        tidyr::gather(key = key, value = value, contains("b_")) %>%
        dplyr::mutate(id = str_extract(key, "[0-9]{1,}"))
    colnames(par) <- c("iter", "name", param, id)
    par[,ncol(par)] <- as.integer(par[,ncol(par)])
    return(par)
}


# Build a dataframe where each row is a unique combination of parameters for each unique combination of sex, site, clone, year, and provenance that occurs in the data (IndSexGroup). #pardf is a dataframe with N rows of parameters for each sum_forcing, site, provenance, clone, and year combination that appear in the data, where N = number of draws from the posterior distribution
build_par_df <- function(mcmcdf, datdf = udf, sex) {

    udf <- dplyr::filter(datdf, Sex==sex)
    mcmcdf$iter <- 1:nrow(mcmcdf)

    singledimpars <- mcmcdf %>%
        dplyr::select(iter, beta, sigma_site, sigma_prov, sigma_clone, sigma_year, contains("kappa"), contains("mean")) %>%
        rename(kappa1 = `kappa[1]`) %>%
        rename(kappa2 = `kappa[2]`)

    siteb <- tidypar(mcmcdf, "b_site", "SiteID")
    provb <- tidypar(mcmcdf, "b_prov", "ProvenanceID")
    cloneb <- tidypar(mcmcdf, "b_clone", "CloneID")
    yearb <- tidypar(mcmcdf, "b_year", "YearID")

    clonemerge <- left_join(udf, cloneb) %>%
        select(-name)
    pardf <- left_join(udf, clonemerge)
    rm(clonemerge, cloneb)
    sitemerge <- left_join(udf, siteb) %>%
        select(-name)
    pardf <- left_join(pardf, sitemerge)
    rm(sitemerge, siteb)
    provmerge <- left_join(udf, provb) %>%
        select(-name)
    pardf <- left_join(pardf, provmerge)
    rm(provmerge, provb)
    yearmerge <- left_join(udf, yearb) %>%
        select(-name)
    pardf <- left_join(pardf, yearmerge)
    rm(yearmerge, yearb)
    pardf <- left_join(pardf, singledimpars)

    # pardf <- left_join(udf, clonemerge) %>%
    #     left_join(provmerge) %>%
    #     left_join(sitemerge) %>%
    #     left_join(yearmerge) %>%
    #     left_join(singledimpars)
}

calc_flowering_prob <- function(forcingaccum, beta, kappa1, kappa2) {
    prob = logistic2(x=forcingaccum, b=beta, c=kappa1) - logistic2(forcingaccum, b=beta, c=kappa2)
    return(prob)
}

# This function calculates the forcing units required for a tree to have a 20% chance of being in a state > 1.
# calcstage2start <- function(p=0.2, beta=betas, kappa1=kappa1, kappa2=kappa2) {
#     x <- (exp(kappa2)*p - exp(kappa2) + exp(kappa1)*p + exp(kappa1))^2 - 4*(p^2)*exp(kappa1+kappa2)
#     y <- -exp(kappa2)*p + exp(kappa2) - exp(kappa1)*p - exp(kappa1)
#     f <- log((-sqrt(x) + y)/(2*p))/beta
#
#     return(f)
# }

calcstageforcing <- function(p=0.2, beta=betas, kappa) {
    prob <- (logit(p) + kappa)/beta
    return(prob)
}

# This function calculates the forcing units required to reach p flowering left
# calcstage2end <- function(p=0.2, beta=betas, kappa1=kappa1, kappa2=kappa2) {
#     x <- (exp(kappa2)*p - exp(kappa2) + exp(kappa1)*p + exp(kappa1))^2 - 4*(p^2)*exp(kappa1+kappa2)
#     y <- -exp(kappa2)*p + exp(kappa2) - exp(kappa1)*p - exp(kappa1)
#     f <- log((sqrt(x) + y)/(2*p))/beta
#
#     return(f)
# }
# this function calculates the number of forcing units at which the probability of having passed out of stage 2 is 80%
# calcstage2end <- function(p=0.8, beta=betas, kappa1=kappa1) {
#     end <- (logit(p) + kappa1)/beta
#     return(end)
# }

# MODEL AND ORIGINAL DATA ############

#model

fmod <- readRDS("slopes_nc_scaled_ristos_FEMALE2019-10-04climatena.rds") %>%
    as.data.frame() %>%
    sample_frac(0.2)

mmod <- readRDS("slopes_nc_scaled_ristos_MALE2019-10-04climatena.rds") %>%
    as.data.frame() %>%
    sample_frac(0.2)

# original data (this code should match relevant bits in run_stan)
phendf <- read_data()

# identify combinations of effects that actually occur

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")


udf <- unique_grouper(fdf, mdf)


#pardf is a dataframe with N rows of parameters for each sum_forcing, site, provenance, clone, and year combination that appear in the data, where N = number of draws from the posterior distribution

fpardf <- build_par_df(mcmcdf = fmod, datdf = fdf, sex = "FEMALE")
mpardf <- build_par_df(mcmcdf = mmod, datdf = mdf, sex = "MALE")

pardf <- rbind(fpardf, mpardf)

rm(fpardf)
rm(mpardf)


# Calculate transformed parameters ################
# add transformed parameters fstart and fend with and without effects
pardf_trans <- pardf %>%
    mutate(betas = b_clone + b_prov + b_site + b_year + beta) %>%
mutate(fbegin_all = calcstageforcing(p=0.2, beta=betas, kappa=kappa1)) %>%
    mutate(fend_all = calcstageforcing(p=0.8, beta=betas, kappa=kappa2)) %>%
    mutate(fbegin_noprov = calcstageforcing(p=0.2, beta=betas-b_prov, kappa=kappa1)) %>%
    mutate(fend_noprov = calcstageforcing(p=0.8, beta=betas-b_prov, kappa=kappa2)) %>%
    mutate(fbegin_nosite = calcstageforcing(p=0.2, beta=betas-b_site, kappa=kappa1)) %>%
    mutate(fend_nosite = calcstageforcing(p=0.8, beta=betas-b_site, kappa=kappa2)) %>%
    mutate(fbegin_pop = calcstageforcing(p=0.2, beta=beta, kappa=kappa1)) %>%
    mutate(fend_pop = calcstageforcing(p=0.8, beta=beta, kappa=kappa2))

pardf_trans <- data.table::as.data.table(pardf_trans)
rm(pardf)

#pardf has all parameters and transformed parameters

tpars <- pardf_trans %>%
    dplyr::select(-forcing_type, -forcing) %>%
    dplyr::select(iter, contains("ID"), Sex, Site, SPU_Name, Clone, Year, contains("Ind"), starts_with("f")) %>%
    gather(key="param", value="sum_forcing", starts_with("f")) %>%
    mutate(side = str_extract(param, ("begin|end"))) %>%
    mutate(effect = str_extract(param, "all|pop|noprov|nosite"))


write.csv(tpars, "transformed_parameters.csv", row.names = FALSE)

# Get data for flowering periods
fbloom <- dplyr::select(fdf, SiteID, ProvenanceID, CloneID, YearID, Site, SPU_Name, Clone, Year, Phenophase, Sex, sum_forcing) %>%
    filter(Phenophase==2) %>%
    rename(FUs = sum_forcing)

mbloom <- dplyr::select(mdf, SiteID, ProvenanceID, CloneID, YearID, Site, SPU_Name, Clone, Year, Phenophase, Sex, sum_forcing) %>%
    filter(Phenophase==2) %>%
    rename(FUs = sum_forcing)

bdat <- rbind(fbloom, mbloom)# %>%
#full_join(tpars)
# bdat has forcing units (FUs) that the model estimates start and end to occur at as well as the actual forcing units (sum_forcing) that flowering was recorded at.

# Plot raw data

#plot data #######################

s2 <- phendf %>%
    filter(Phenophase==2) %>%
    group_by(Index) %>%
    summarise(s2 = min(DoY))

s3 <- phendf %>%
    filter(Phenophase==3) %>%
    group_by(Index) %>%
    summarise(s3 = min(DoY))

phendf <- left_join(phendf, s2) %>%
    left_join(s3)

# This plot shows the cumulative distribution of accumulated forcing units for phenological states 2 and 3 for the collected data
ggplot(filter(phendf, DoY == s2 | DoY==s3), aes(x=sum_forcing, color=Sex, linetype=Phenophase)) +
    stat_ecdf() +
    facet_wrap(Site ~ SPU_Name) +
    scale_color_viridis_d() +
    theme_bw(base_size=18) +
    ggtitle("Cumulative phase 2 and 3 flowering phenology data") +
    theme(legend.position = "top") +
    ylab(NULL) +
    xlab("Accumulated forcing units")


ggplot(filter(phendf, Phenophase==2 & Site %in% c("PGTIS", "Kalamalka", "PRT", "Vernon")), aes(x=sum_forcing, color=Sex)) +
    stat_ecdf(size=1.1) +
    facet_grid(SPU_Name ~ Site) +
    scale_color_viridis_d() +
    theme_bw(base_size=18) +
    ggtitle("Flowering - Provenance comparison")

ggplot(filter(phendf, Phenophase==2 & SPU_Name %in% c("Bulkley Valley Low", "Central Plateau Low", "Prince George Low")), aes(x=sum_forcing, color=Sex)) +
    stat_ecdf() +
    facet_grid(Site ~ SPU_Name) +
    scale_color_viridis_d() +
    theme_bw(base_size=20) +
    ggtitle("Flowering - Site comparison")

# Plot transformed parameters #############

# halfparams <- filter(tpars, str_detect(param, "half"))
#
# ggplot(halfparams, aes(x=sum_forcing, fill=param)) +
#     geom_density(alpha=0.5) +
#     scale_fill_viridis_d()

ggplot(tpars, aes(x=sum_forcing, color=effect, linetype=side)) +
    geom_density() +
    facet_grid(SPU_Name ~ Sex) +
    scale_color_viridis_d() +
    theme_bw()


# ph1 <- filter(tpars, param %in% c("fhalf1", "fhalf2")) %>%
#     mutate(transition = case_when(param=="fhalf1" ~ 1,
#                                   param=="fhalf2" ~ 2))
# ph2 <- filter(tpars, param %in% c("fhalf1_noprov", "fhalf2_noprov")) %>%
#     rename(no_prov=param, sum_forcing_noprov=sum_forcing) %>%
#     mutate(transition = case_when(no_prov=="fhalf1_noprov" ~ 1,
#                                   no_prov=="fhalf2_poprov" ~ 2))

provcomp <- filter(tpars, effect %in% c("all", "no_prov"))

ggplot(provcomp, aes(x=sum_forcing, fill=effect, linetype=side)) +
    geom_density(alpha=0.5) +
    stat_ecdf(data=bdat, aes(x=FUs), inherit.aes=FALSE) +
    scale_fill_viridis_d(option="B") +
    facet_grid(SPU_Name ~ Sex) +
    ggtitle("Start and end by provenance") +
    xlab("Forcing units") +
    theme_bw(base_size=18) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(legend.position= "top")

temp <- filter(provcomp, effect=="all")
ggplot(provcomp, aes(x=sum_forcing, fill=effect, linetype=side)) +
    geom_density(alpha=0.5) +
    stat_ecdf(data=bdat, aes(x=FUs), inherit.aes=FALSE) +
    scale_fill_viridis_d(option="B") +
    #facet_grid(SPU_Name ~ Sex) +
    ggtitle("Start and end") +
    xlab("Forcing units") +
    theme_bw(base_size=18) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(legend.position= "top")

sitecomp <- filter(tpars, effect %in% c("all", "no_site")) %>%
    filter(Sex=="FEMALE")

ggplot(sitecomp, aes(x=sum_forcing, fill=effect, linetype=side)) +
    geom_density(alpha=0.5) +
    #stat_ecdf(data=bdat, aes(x=FUs), inherit.aes=FALSE) +
    scale_fill_viridis_d(option="A") +
    facet_grid(Site ~ .) +
    ggtitle("Start and end by site") +
    xlab("Forcing units") +
    theme_bw(base_size=18) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(legend.position= "top")




# calculate HPDIs ###############

hpd_lower = function(x, prob) rethinking::HPDI(x, prob)[1]
hpd_upper = function(x, prob) rethinking::HPDI(x, prob)[2]

full <- group_by(tpars, param, Sex) %>%
    summarize(x=hpd_lower(sum_forcing, prob=.99), xend=hpd_upper(sum_forcing, prob=.99)) %>%
    mutate(interval=".99")

fifty <- group_by(tpars, param, Sex) %>%
    summarize(x=hpd_lower(sum_forcing, prob=.5), xend=hpd_upper(sum_forcing, prob=.5)) %>%
    mutate(interval="0.5")

hpd_df <- full_join(full, fifty) %>%
    gather(key=end, value=value, x, xend) %>%
    arrange(interval)


# How to do an interval plot
ggplot(filter(hpd_df, interval==".99"), aes(x=value, y=0, color=Sex, size=interval)) +
    geom_line() +
    geom_line(data = filter(hpd_df, interval=="0.5"), aes(x=value, y=0)) +
    facet_grid(rows=vars(param, Sex)) +
    theme_bw() +
    theme(strip.text.y = element_text(angle = 0)) +
    scale_color_viridis_d()

# This is a good plot
ggplot(bdat, aes(x=sum_forcing, fill=Sex, linetype=param)) +
    geom_density(alpha=0.5) +
    geom_density(aes(x=FUs, alpha=0.5, linetype=NULL))+
    theme_bw() +
    facet_grid(SPU_Name ~ .) +
    theme(strip.text.y = element_text(angle = 0)) +
    scale_fill_viridis_d()+
    ggtitle("Modeled start and end of flowering with flowering period data", subtitle = "by provenance. fstart = 20% begun, fend = 80% done")

# This is a good plot 2
ggplot(bdat, aes(x=sum_forcing, fill=Sex, linetype=param)) +
    geom_density(alpha=0.5) +
    stat_ecdf(aes(x=FUs, color=Sex), alpha=.8)+
    theme_bw() +
    facet_grid(Site ~ .) +
    theme(strip.text.y = element_text(angle = 0)) +
    scale_fill_viridis_d()+
    scale_color_viridis_d() +
    ggtitle("Modeled start and end of flowering with cumulative flowering period data", subtitle = "by site. fstart = 20% begun, fend = 80% done")

ggplot(fdat, aes(x=sum_forcing, fill=param)) +
    geom_density( alpha=0.5) +
    geom_density(aes(x=FUs, fill=NULL))+
    theme_bw() +
    facet_grid(Site ~ .) +
    theme(strip.text.y = element_text(angle = 0)) +
    scale_fill_viridis_d()+
    ggtitle("Modeled start and end of receptivity with receptive period data ", subtitle = "by site. fstart = 20% begun, fend = 80% done")

ggplot(fdat, aes(x=sum_forcing, fill=param)) +
    geom_density( alpha=0.5) +
    geom_density(aes(x=FUs, fill=NULL))+
    theme_bw() +
    facet_wrap(Year ~ ., scales="free") +
    theme(strip.text.y = element_text(angle = 0)) +
    scale_fill_viridis_d()+
    ggtitle("Modeled start and end of receptivity with receptive period data ", subtitle = "by year. fstart = 20% begun, fend = 80% done")

ggplot(filter(tpars_long, param %in% c("fstart", "fend")), aes(x=FUs, fill=SPU_Name, linetype=param)) +
    geom_density(alpha=.2) +
    theme_bw()

ggplot(filter(tpars_long, param %in% c("fstart", "fend")), aes(x=FUs, fill=Site, linetype=param)) +
    geom_density(alpha=.2) +
    theme_bw()

ggplot(filter(tpars_long, param %in% c("fstart", "fend")), aes(x=FUs, group=as.factor(CloneID))) +
    geom_density(alpha=.2) +
    theme_bw() +
    facet_grid(param ~ .) +
    theme(
        plot.background = element_blank(),
        #panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.4),
        axis.ticks = element_line(size = 0.3),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(0.9)),
        strip.placement = "outside",
        # strip.background = element_rect(fill = "gray95", color = NA),
        panel.spacing = unit(1.5, "lines"),
        legend.position = "right",
        legend.background = element_blank(),
        legend.text = element_text(size = 13),
        legend.text.align = 0,
        legend.key = element_blank()
    )

tpars_wide <- select(tpars_wide, starts_with("f"))
mcmc_areas(tpars_wide)

# Compare provenances

full_prov <- group_by(tpars, Sex, SPU_Name, param) %>%
    summarize(x=hpd_lower(sum_forcing, prob=.99), xend=hpd_upper(sum_forcing, prob=.99)) %>%
    mutate(interval=".99")

fifty_prov <- group_by(tpars, Sex, SPU_Name, param) %>%
    summarize(x=hpd_lower(sum_forcing, prob=.5), xend=hpd_upper(sum_forcing, prob=.5)) %>%
    mutate(interval="0.5")

hpd_df_prov <- full_join(full_prov, fifty_prov) %>%
    gather(key=end, value=value, x, xend) %>%
    arrange(interval)

ggplot(filter(hpd_df_prov, interval==".99"), aes(x=value, y=SPU_Name, color=Sex, size=interval)) +
    geom_line(alpha=0.5) +
    geom_line(data = filter(hpd_df_prov, interval=="0.5"), aes(x=value, y=SPU_Name)) +
    facet_grid(rows=vars(param, Sex) ) +
    theme_bw() +
    scale_color_viridis_d()

# Calculate prob of flowering (S=2) over forcing units ############

small_pars <- pardf_trans %>%
    group_by(IndSexGroup, Year) %>%
    sample_n(30, replace=TRUE) %>%
    mutate(beta_tot = beta + b_site + b_prov + b_clone + b_year)
small_pars$rows <- rownames(small_pars)
fus <- seq(from=0, to=30, by=1.5)

fucalcstructure <- expand.grid(small_pars$rows, fus)
colnames(fucalcstructure) <- c("rows", "fus")

small_pars <- full_join(small_pars, fucalcstructure)

small_pars$stage2prob <- logistic2(small_pars$fus, small_pars$beta, small_pars$kappa1) -
    logistic2(small_pars$fus, small_pars$beta, small_pars$kappa2)

ggplot(small_pars, aes(x=fus, y=stage2prob, color=Sex, group=rows)) +
    geom_line(alpha=.05) +
    facet_wrap(Sex ~ .) +
    scale_color_viridis_d() +
    xlab("Accumulated forcing units") +
    theme_bw(base_size = 20)


fbloomd <- dplyr::select(fdf, SiteID, ProvenanceID, CloneID, YearID, Site, SPU_Name, Clone, Year, Phenophase_Derived, Sex, DoY) %>%
    filter(Phenophase_Derived==2) %>%
    rename(DoY_obs = DoY)

mbloomd <- dplyr::select(mdf, SiteID, ProvenanceID, CloneID, YearID, Site, SPU_Name, Clone, Year, Phenophase_Derived, Sex, DoY) %>%
    filter(Phenophase_Derived==2) %>%
    rename(DoY_obs = DoY)

bdatd <- rbind(fbloomd, mbloomd)# %>%
full_join(tpars)

# Calculate flowering period ############

sum_forcing <- seq(from=4, to=35, by=.5)
names(sum_forcing) <- sum_forcing


flowerprob <- map_dfc(sum_forcing, calc_flowering_prob, beta=pardf_trans$betas, kappa1=pardf_trans$kappa1, kappa2=pardf_trans$kappa2)

fprob <- data.frame(pardf, flowerprob) %>%
    gather(key=sum_forcing, value=flowerprob, starts_with("X")) %>%
    separate(sum_forcing, into=c("drop", "sum_forcing"), sep=1) %>%
    select(-drop)
fprob$index <- group_indices(fprob, IndSexGroup, iter)

ggplot(fprob, aes(x=sum_forcing, y=flowerprob, group=index, color=Sex)) +
    geom_line(alpha=0.3) +
    facet_wrap("SPU_Name")


