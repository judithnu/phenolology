#### posterior predictive checks

samples = 200
forcingtype = 'scaled_ristos'

# Create datasets holding parameter values - fstart, fend, and f50s

library(tidyverse)
library(rstan)
library(gtools)


# FUNCTIONS ####################################

# Add grouping columns that stan will tolerate.
stanindexer <- function(df) {
    df$CloneID <- group_indices(df, Clone)
    df$OrchardID <- group_indices(df, Orchard)
    df$ProvenanceID <- group_indices(df, SPU_Name)
    df$SiteID <- group_indices(df, Site)
    df$YearID <- group_indices(df, Year)
    df$TreeID <- group_indices(df,TreeUnique)
    return(df)
}

logistic2 <- function(x,b,c) {
    1/(1 + exp(-(b * (x-(c/b)))))
}

pdflogistic <- function(x,b,c) {
    num <- b * exp(-b*(x-(c/b)))
    denom <- (1 + exp(-(b * (x-(c/b)))))^2
    val <- num/denom
    return(val)
}

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

# Filter the phenology dataframe by sex. Add grouping vars a la stan and extract unique combinations of grouping variables. df is a dataframe and sex is "MALE" or "FEMALE". Stan requires grouping vars to be consecutive integers, so indexes will differ/represent different underlying groups for any groups that are not identical across sexes
splitsex <- function(df = phendf, sex) {
    df <- filter(phendf, Sex == sex)
    # add grouping columns to match with stan output
    df <- stanindexer(df)
}

# join male and female sex
unique_grouper <- function(df1 = fdf, df2 = mdf) {
    df <- full_join(df1, df2)
    #extract unique groups
    udf <- df %>%
        dplyr::select(Sex, Site, SiteID, SPU_Name, ProvenanceID, Clone, CloneID, Year, YearID) %>%
        distinct()
# Group by individual clone and group by individuals WITH sex
    udf$IndGroup <- group_indices(udf, SiteID, ProvenanceID, YearID, CloneID)
    udf$IndSexGroup <- group_indices(udf, Sex, IndGroup)
    return(udf)
}

# Build a dataframe where each row is a unique combination of parameters for each unique combination of sex, site, clone, year, and provenance that occurs in the data (IndSexGroup). #pardf is a dataframe with N rows of parameters for each sum_forcing, site, provenance, clone, and year combination that appear in the data, where N = number of draws from the posterior distribution
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

calc_flowering_prob <- function(forcingaccum, beta, kappa1, kappa2) {
    prob = logistic2(x=forcingaccum, b=beta, c=kappa1) - logistic2(forcingaccum, b=beta, c=kappa2)
    return(prob)
}

# This function calculates the forcing units required to reach p flowering started
calcstage2start <- function(p=0.2, beta=betas, kappa1=kappa1, kappa2=kappa2) {
    x <- (exp(kappa2)*p - exp(kappa2) + exp(kappa1)*p + exp(kappa1))^2 - 4*(p^2)*exp(kappa1+kappa2)
    y <- -exp(kappa2)*p + exp(kappa2) - exp(kappa1)*p - exp(kappa1)
    f <- log((-sqrt(x) + y)/(2*p))/beta

    return(f)
}

# This function calculates the forcing units required to reach p flowering left
calcstage2end <- function(p=0.2, beta=betas, kappa1=kappa1, kappa2=kappa2) {
    x <- (exp(kappa2)*p - exp(kappa2) + exp(kappa1)*p + exp(kappa1))^2 - 4*(p^2)*exp(kappa1+kappa2)
    y <- -exp(kappa2)*p + exp(kappa2) - exp(kappa1)*p - exp(kappa1)
    f <- log((sqrt(x) + y)/(2*p))/beta

    return(f)
}

# MODEL AND ORIGINAL DATA ############

#model
fmod <- readRDS("FEMALE_slopes_scaled.rds") %>%
    as.data.frame()

mmod <- readRDS("MALE_slopes_scaled.rds") %>%
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
    unique()


# identify combinations of effects that actually occur

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")


udf <- unique_grouper(fdf, mdf)

#pardf is a dataframe with N rows of parameters for each sum_forcing, site, provenance, clone, and year combination that appear in the data, where N = number of draws from the posterior distribution

fpardf <- build_par_df(mcmcdf = fmod, datdf = udf, sex = "FEMALE")
mpardf <- build_par_df(mcmcdf = mmod, datdf = udf, sex = "MALE")
pardf <- rbind(fpardf, mpardf) %>%
    group_by(IndSexGroup) %>%
    sample_n(samples) #downsample

# Calculate transformed parameters ################
# add transformed parameters fstart, fend, and f50s
pardf_trans <- pardf %>%
    mutate(betas = b_clone + b_prov + b_site + b_year + beta) %>%
    # mutate(fstart = (logit(.2) + kappa1)/betas) %>%
    # mutate(fend = (logit(.8) + kappa2)/betas) %>%
    # mutate(fhalf1 = kappa1/betas) %>%
    # mutate(fhalf2 = kappa2/betas) %>%
    #mutate(fhalf1_pop = kappa1/beta) %>%
    #mutate(fhalf2_pop = kappa2/beta) %>%
    #mutate(fstart1_pop = (logit(.2) + kappa1)/beta) %>%
    # mutate(fhalf1_noprov = kappa1/(betas-b_prov)) %>%
    # mutate(fhalf2_noprov = kappa2/(betas-b_prov)) %>%
    # mutate(fhalf1_nosite = kappa1/(betas-b_site)) %>%
    # mutate(fstart_noprov = (logit(.2) + kappa1)/(betas-b_prov))
mutate(fbegin_all = calcstage2start(p=0.2, beta=betas, kappa1=kappa1, kappa2=kappa2)) %>%
    mutate(fend_all = calcstage2end(p=0.2, beta=betas, kappa1=kappa1, kappa2=kappa2)) %>%
    mutate(fbegin_noprov = calcstage2start(p=0.2, beta=betas-b_prov, kappa1=kappa1, kappa2=kappa2)) %>%
    mutate(fend_noprov = calcstage2end(p=0.2, beta=betas-b_prov, kappa1=kappa1, kappa2=kappa2)) %>%
    mutate(fbegin_nosite = calcstage2start(p=0.2, beta=betas-b_site, kappa1=kappa1, kappa2=kappa2)) %>%
    mutate(fend_nosite = calcstage2end(p=0.2, beta=betas-b_site, kappa1=kappa1, kappa2=kappa2)) %>%
    mutate(fbegin_pop = calcstage2start(p=0.2, beta=beta, kappa1=kappa1, kappa2=kappa2)) %>%
    mutate(fend_pop = calcstage2end(p=0.2, beta=beta, kappa1=kappa1, kappa2=kappa2))

#pardf has all parameters and transformed parameters

tpars <- dplyr::select(pardf_trans, iter, contains("ID"), Sex, Site, SPU_Name, Clone, Year, contains("Ind"), starts_with("f")) %>%
    gather(key="param", value="sum_forcing", starts_with("f")) %>%
    mutate(side = case_when(str_detect(param, "begin") ~ "begin",
                            str_detect(param, "end") ~ "end")) %>%
    mutate(effect = case_when(str_detect(param, "all") ~ "all",
                              str_detect(param, "pop") ~ "pop",
                              str_detect(param, "noprov") ~ "no_prov",
                              str_detect(param, "nosite") ~ "no_site"))

#tpars$param <- factor(tpars$param)
#tpars$param = factor(tpars$param,levels(tpars$param)[c(1,3,2,4)])

#write.csv(tpars, "transformed_parameters.csv", row.names = FALSE)

# Plot raw data

#plot data #######################
ggplot(phendf, aes(x=sum_forcing, color=Sex)) +
    stat_ecdf(size=1.1) +
    facet_wrap(SPU_Name ~ Site) +
    scale_color_viridis_d() +
    theme_bw(base_size=20) +
    ggtitle("Flowering phenology data")

ggplot(filter(phendf, Phenophase_Derived==2 & Site %in% c("PGTIS", "Kalamalka", "PRT", "Vernon")), aes(x=sum_forcing, color=Sex)) +
    stat_ecdf(size=1.1) +
    facet_grid(SPU_Name ~ Site) +
    scale_color_viridis_d() +
    theme_bw(base_size=20) +
    ggtitle("Flowering - Provenance comparison")

ggplot(filter(phendf, Phenophase_Derived==2 & SPU_Name %in% c("Bulkley Valley Low", "Central Plateau Low", "Prince George Low")), aes(x=sum_forcing, color=Sex)) +
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

sitecomp <- filter(tpars, effect %in% c("all", "no_site"))

ggplot(sitecomp, aes(x=sum_forcing, fill=effect, linetype=side)) +
    geom_density(alpha=0.5) +
    stat_ecdf(data=bdat, aes(x=FUs), inherit.aes=FALSE) +
    scale_fill_viridis_d(option="A") +
    facet_grid(Site ~ Sex) +
    ggtitle("Start and end by site") +
    xlab("Forcing units") +
    theme_bw(base_size=18) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(legend.position= "top")

#present
# provhalf <- full_join(ph1, ph2)
# ggplot(provhalf, aes(x=sum_forcing, fill="withprov", group=param)) +
#     geom_density() +
#     geom_density(aes(x=sum_forcing_noprov, fill="no_prov", linetype=param, group=no_prov), alpha=0.5) +
#     scale_fill_viridis_d(option="B") +
#     facet_grid(SPU_Name ~ Sex) +
#     ggtitle("Provenance effect on 50% forcing unit requirements") +
#     xlab("Forcing units") +
#     theme_bw(base_size=18) +
#     theme(strip.text.y = element_text(angle = 0)) +
#     theme(legend.position= "top")



hpd_lower = function(x, prob) rethinking::HPDI(x, prob)[1]
hpd_upper = function(x, prob) rethinking::HPDI(x, prob)[2]

# calculate HPDIs

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


# Get data for flowering periods
fbloom <- dplyr::select(fdf, SiteID, ProvenanceID, CloneID, YearID, Site, SPU_Name, Clone, Year, Phenophase_Derived, Sex, sum_forcing) %>%
    filter(Phenophase_Derived==2) %>%
    rename(FUs = sum_forcing)

mbloom <- dplyr::select(mdf, SiteID, ProvenanceID, CloneID, YearID, Site, SPU_Name, Clone, Year, Phenophase_Derived, Sex, sum_forcing) %>%
    filter(Phenophase_Derived==2) %>%
    rename(FUs = sum_forcing)

bdat <- rbind(fbloom, mbloom)# %>%
    full_join(tpars)
# bdat has forcing units (FUs) that the model estimates start and end to occur at as well as the actual forcing units (sum_forcing) that flowering was recorded at.

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


