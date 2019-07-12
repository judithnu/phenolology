#### posterior predictive checks

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
    df$TreeUnique <- group_indices(df,TreeID)
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

tidypar <- function(stanframe, param, id) {

    #take a stan model dataframe and create a tidy dataframe for a given parameter. takes a dataframe, a string with the parameter name, and a string describing the param ID (e.g. par = "Site" and id="SiteID")

    #wide to long format for parameters
    par <- stanframe %>% dplyr::select(contains(param), iter) %>%
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
    mcmcdf$iter <- 1:nrow(fmod)

    singledimpars <- mcmcdf %>%
        dplyr::select(iter, beta, sigma_site, sigma_prov, sigma_clone, sigma_year, contains("kappa")) %>%
        rename(kappa1 = `kappa[1]`) %>%
        rename(kappa2 = `kappa[2]`)

    kappa <- dplyr::select(mcmcdf, contains("kappa"))
    siteb <- tidypar(mcmcdf, "b_site", "SiteID")
    provb <- tidypar(mcmcdf, "b_prov", "ProvenanceID")
    cloneb <- tidypar(mcmcdf, "b_clone", "CloneID")
    yearb <- tidypar(mcmcdf, "b_year", "YearID")

    clonemerge <- left_join(udf, cloneb)
    provmerge <- left_join(udf, provb)
    sitemerge <- left_join(udf, siteb)
    yearmerge <- left_join(udf, yearb)

    pardf <- data.frame(udf,
                        iter = clonemerge$iter,
                        b_clone = clonemerge$b_clone,
                        b_prov = provmerge$b_prov,
                        b_site = sitemerge$b_site,
                        b_year = yearmerge$b_year) %>%
        left_join(singledimpars)
}


# MODEL AND ORIGINAL DATA ############

#model
fmod <- readRDS("FEMALE_slopes_scaled.rds") %>%
    as.data.frame()

mmod <- readRDS("MALE_slopes_scaled.rds") %>%
    as.data.frame()

#data
phenology_data <- read.csv("data/phenology_heatsum.csv",
                           stringsAsFactors = FALSE,
                           header = TRUE) %>%
    filter(forcing_type=="scaled_ristos")

clim <- read.csv("data/all_clim_PCIC.csv",
                 stringsAsFactors = FALSE, header=TRUE) %>%
    filter(forcing_type=="scaled_ristos")

## provenance
SPU_dat <- read.csv("../research_phd/data/OrchardInfo/LodgepoleSPUs.csv",
                    header=TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select(SPU_Name, Orchard)
SPU_dat$SPU_Name <- str_split(SPU_dat$SPU_Name, " \\(", simplify=TRUE)[,1]

# Data Processing ##################
# join provenance and phenology data

phendf <- dplyr::left_join(phenology_data, SPU_dat) %>%
    unique()

# identify combinations of effects that actually occur

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")

udf <- unique_grouper(fdf, mdf)

#pardf is a dataframe with N rows of parameters for each sum_forcing, site, provenance, clone, and year combination that appear in the data, where N = number of draws from the posterior distribution

fpardf <- build_par_df(mcmcdf = fmod, datdf = udf, sex = "FEMALE")
mpardf <- build_par_df(mcmcdf = mmod, datdf = udf, sex = "MALE")
pardf <- rbind(fpardf, mpardf)

# add transformed parameters fstart, fend, and f50s
pardf_trans <- pardf %>%
    mutate(betas = b_clone + b_prov + b_site + b_year + beta) %>%
    mutate(fstart = (logit(.2) + kappa1)/betas) %>%
    mutate(fend = (logit(.8) + kappa2)/betas) %>%
    mutate(fhalf1 = kappa1/betas) %>%
    mutate(fhalf2 = kappa2/betas)

#pardf has all parameters and transformed parameters

tpars <- dplyr::select(pardf_trans, iter, contains("ID"), Sex, Site, SPU_Name, Clone, Year, contains("Ind"), starts_with("f")) %>%
    gather(key="param", value="sum_forcing", starts_with("f"))

tpars$param <- factor(tpars$param)
tpars$param = factor(tpars$param,levels(tpars$param)[c(1,3,2,4)])



# Get data for flowering periods
fbloom <- dplyr::select(fdf, SiteID, ProvenanceID, CloneID, YearID, Site, SPU_Name, Clone, Year, Phenophase_Derived, Sex, sum_forcing) %>%
    filter(Phenophase_Derived==2) %>%
    rename(FUs = sum_forcing)

mbloom <- dplyr::select(mdf, SiteID, ProvenanceID, CloneID, YearID, Site, SPU_Name, Clone, Year, Phenophase_Derived, Sex, sum_forcing) %>%
    filter(Phenophase_Derived==2) %>%
    rename(FUs = sum_forcing)

bdat <- rbind(fbloom, mbloom) %>%
    full_join(tpars) %>%
    filter(param %in% c("fstart", "fend"))
# bdat has forcing units (FUs) that the model estimates start and end to occur at as well as the actual forcing units (sum_forcing) that flowering was recorded at.

# Plot transformed parameters #############

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
