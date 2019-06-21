#### posterior predictive checks

# Create datasets holding parameter values - fstart, fend, and f50s


library(tidyverse)
library(rstan)
library(rethinking)

# FUNCTIONS ####################################
stanindexer <- function(df) {
    df$CloneID <- group_indices(df, Clone)
    df$OrchardID <- group_indices(df, Orchard)
    df$ProvenanceID <- group_indices(df, SPU_Name)
    df$SiteID <- group_indices(df, Site)
    df$YearID <- group_indices(df, Year)
    df$Tree <- group_indices(df,TreeID)
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

# MODEL AND ORIGINAL DATA ############
#model
fmod <- readRDS("female_slopes.rds") %>%
    as.data.frame()
#mfit.stan <- readRDS("male_slopes.rds")

#mmod <- as.data.frame(mfit.stan)

#far <- as.array(ffit.stan)

#data
phenology_data <- read.csv("data/phenology_heatsum.csv",
                           stringsAsFactors = FALSE, header = TRUE) %>%
    filter(forcing_type=="ristos")

clim <- read.csv("data/all_clim_PCIC.csv", stringsAsFactors = FALSE, header=TRUE) %>%
    filter(forcing_type=="ristos") %>%
    filter(DoY %in% c(90:190))

## provenance
SPU_dat <- read.csv("../research_phd/data/OrchardInfo/LodgepoleSPUs.csv",
                    header=TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select(SPU_Name, Orchard) %>%
    mutate(SPU_Name = str_)

# Data Processing ##################
# join provenance and phenology data

phendf <- phenology_data %>%
    na.omit()
phendf <- dplyr::left_join(phenology_data, SPU_dat) %>%
    unique()

# separate phenology into male and female dataframes and turn factors into integers
fdf <- filter(phendf, Sex == "FEMALE")
fdf <- stanindexer(fdf)
fdf$groupID <- group_indices(fdf, SiteID, ProvenanceID, YearID, CloneID)
fdf$recordID <- group_indices(fdf, groupID, Date )
# mdf <- filter(phendf, Sex == "MALE")
# mdf <- stanindexer(mdf)

# identify combinations of effects & predictor (sum_forcing) that actually occur
# ufdf <- fdf %>%
#     select(groupID, recordID, Site, SiteID, SPU_Name, ProvenanceID, Clone, CloneID, Year, YearID, sum_forcing, DoY) %>%
#     distinct()

#identify combinations of effects that actually occur
ufdf <- fdf %>%
    dplyr::select(groupID, Site, SiteID, SPU_Name, ProvenanceID, Clone, CloneID, Year, YearID) %>%
    distinct()

# Build a dataframe where each row is a unique combination of parameters for each unique combination of forcing, site, clone, year, and provenance #########
#split by param

fmod$iter <- 1:nrow(fmod)
# fmod_sampled <- fmod[draws,]

singledimpars <- fmod %>%
    dplyr::select(iter, beta, sigma_site, sigma_prov, sigma_clone, sigma_year, contains("kappa")) %>%
    rename(kappa1 = `kappa[1]`) %>%
    rename(kappa2 = `kappa[2]`)



kappa <- dplyr::select(fmod, contains("kappa"))
siteb <- tidypar(fmod, "b_site", "SiteID")
provb <- tidypar(fmod, "b_prov", "ProvenanceID")
cloneb <- tidypar(fmod, "b_clone", "CloneID")
yearb <- tidypar(fmod, "b_year", "YearID")

clonemerge <- left_join(ufdf, cloneb)
provmerge <- left_join(ufdf, provb)
sitemerge <- left_join(ufdf, siteb)
yearmerge <- left_join(ufdf, yearb)

pardf <- data.frame(ufdf,
                       iter = clonemerge$iter,
                       b_clone = clonemerge$b_clone,
                       b_prov = provmerge$b_prov,
                       b_site = sitemerge$b_site,
                       b_year = yearmerge$b_year) %>%
    left_join(singledimpars)

#pardf is a dataframe with N rows of parameters for each sum_forcing, site, provenance, clone, and year combination that appear in the data, where N = number of draws from the posterior distribution

# add transformed parameters fstart, fend, and f50s
pardf_trans <- pardf %>%
    mutate(betas = b_clone + b_prov + b_site + b_year + beta) %>%
    mutate(fstart = (logit(.2) + kappa1)/betas) %>%
    mutate(fend = (logit(.8) + kappa2)/betas) %>%
    mutate(fhalf1 = kappa1/betas) %>%
    mutate(fhalf2 = kappa2/betas)

tpars_wide <- dplyr::select(pardf_trans, contains("ID"), Site, SPU_Name, Clone, Year, starts_with("f"))

tpars_long <- tpars_wide %>%
    gather(key="param", value="FUs", starts_with("f"))

provpaster <- function(var) {
    paste(var, ProvenanceID, sep="")
}
tpars_wide_prov <- tpars_wide %>%
    mutate(fstartprov = paste("fstart", ProvenanceID, sep=""))

# Plot transformed parameters #############

hpd_lower = function(x, prob) rethinking::HPDI(x, prob=.9)[1]
hpd_upper = function(x, prob) rethinking::HPDI(x, prob=.9)[2]


full <- group_by(tpars, param) %>%
    summarize(x=hpd_lower(FUs, prob=1), xend=hpd_upper(FUs, prob=1)) %>%
    mutate(interval="full")

fifty <- group_by(tpars, param) %>%
    summarize(x=hpd_lower(FUs, prob=1), xend=hpd_upper(FUs, prob=1)) %>%
    mutate(interval="fifty")

hpd_df <- full_join(full, fifty) %>%
    gather(key=end, value=value, x, xend) %>%
    arrange(interval)

ggplot(hpd_df, aes(x=value, y=param, size=interval)) +
    geom_line()


ggplot(filter(tpars, param %in% c("fstart", "fend")), aes(x=FUs,y=SPU_Name, color=param)) +
    geom_point() +
    theme_bw()

ggplot(sfull, aes(x))



ggplot(dplyr)
ggplot(dplyr::filter(tpars, !param=="forcingdiff"), aes(x=param, y=value, color=Provenance)) +
    geom_violin(draw_quantiles = .5) +
    scale_color_viridis_d()

ggplot(filter(tpars, !param=="forcingdiff"), aes(x=as.factor(CloneID), y=value)) +
    geom_line(alpha=0.3) +
    facet_wrap(param ~ .)

summary <- tpars %>%
    group_by(param, CloneID) %>%
    summarise(mean(value))

# Now simulate state predictions ##################

# calculate slope and phi in columns
pred_df <- pardf %>%
    mutate(betas = b_clone + b_prov + b_site + b_year + beta) %>%
    mutate(phi = betas * sum_forcing) %>%
    mutate(t1 = logistic2(sum_forcing, betas, kappa1)) %>%
    mutate(t2 = logistic2(sum_forcing, betas, kappa2)) %>%
    mutate(p2 = pdflogistic(sum_forcing, betas, kappa1)) %>%
    mutate(p3 = pdflogistic(sum_forcing, betas, kappa2)) %>%
    mutate(f501 = kappa1/betas) %>%
    mutate(f502 = kappa2/betas) %>%
distinct()



#simulate state (1 per real obs)
state_pred <- c()
for (i in 1:nrow(pred_df)) {
    state_pred[i] <- rethinking::rordlogit(1,
                                      phi = pred_df$phi[i],
                                      a = c(pred_df$kappa1[i], pred_df$kappa2[i]))
}

pred_df <- data.frame(pred_df, state_pred=state_pred)

#merge in real data

real_state <- fdf %>%
    dplyr::select(recordID, Phenophase_Derived)

realpred <- full_join(pred_df, real_state) %>%
    select(-contains("b_")) %>%
    select(-contains("beta")) %>%
    select(-contains("kappa")) %>%
    select(-contains("sigma")) %>%
   gather(key=outcome, value=state, state_pred, Phenophase_Derived)

mrealpred <- realpred
frealpred <- realpred

ggplot(realpred, aes(x=sum_forcing, fill=as.factor(outcome))) +
    geom_histogram(alpha=0.6) +
    facet_grid(state ~ .) +
    scale_fill_viridis_d() +
    ggtitle("Male state data and predictions") +
    theme_bw()

m2rp <- filter(mrealpred, state==2)
f2rp <- filter(frealpred, state==2 & outcome=="state_pred")
ggplot(f2rp, aes(x=f501, fill=SPU_Name, linetype=outcome)) +
    geom_density(alpha=0.4) +
    #facet_wrap("Year") +
    scale_fill_viridis_d() +
    ggtitle("Female 50% first transition")

ggplot(f2rp, aes(x=DoY, fill=SPU_Name, linetype=outcome)) +
    geom_histogram(alpha=0.4) +
    facet_wrap(Year ~ Site, scales="free") +
    scale_fill_viridis_d() +
    ggtitle("Female 50% first transition")

ggplot(f2rp, aes(x=f502, fill=SPU_Name, linetype=outcome)) +
    geom_density(alpha=0.1) +
    #facet_wrap("Year") +
    scale_fill_viridis_d() +
    ggtitle("Female 50% second transition")

# calculate first and last day flowering. NEED MORE SIMULATIONS FOR THIS TO WORK.

#flowering <- filter(realpred, state=="2")


# hpd_lower = function(x) HPDI(x, prob=.9)[1]
# hpd_upper = function(x) HPDI(x, prob=.9)[2]
# df_hpd <- group_by(pred, state) %>%
#     summarise(x=hpd_lower(sum_forcing), xend=hpd_upper(sum_forcing))
#
# # calculate 90% HPDI
# flowering <- filter(pred, state=="2")
#
#
# ggplot(realpred, aes(x=sum_forcing, group=as.factor(draw), color=outcome)) +
#     stat_ecdf(alpha=0.5) +
#     facet_grid(state ~ .)  +
#     #geom_segment(data=df_hpd, aes(x=x, xend=xend, y=0, yend=0, size=3))
#     geom_vline(data=df_hpd, aes(xintercept=x)) +
#     geom_vline(data=df_hpd, aes(xintercept=xend))
#
# ggplot(realpred, aes(x=sum_forcing, fill=outcome)) +
#     geom_density(adjust=3, alpha=0.7) +
#     facet_grid(state ~ .) +
#     ggtitle("female strobili predicted and measured states with 90% HPDI") +
#     geom_vline(data = df_hpd, aes(xintercept=x)) +
#     geom_vline(data = df_hpd, aes(xintercept=xend)) +
#     scale_fill_viridis_d()

## calculate transitions #####################


# ts <- pardf %>%
#     mutate(betas = b_clone + b_prov + b_site + b_year + beta) %>%
#     mutate(t1 = logistic2(sum_forcing, betas, kappa1)) %>%
#     mutate(t2 = logistic2(sum_forcing, betas, kappa2))


plot1 <- ggplot(ts, aes(x=sum_forcing, y=t1, color=draw)) +
    geom_point(alpha=0.5)

# Calculate transitions

fmodsum <- rstan::summary(ffit.stan)
fmodsum <- as.data.frame(fmodsum$summary)[,c(1,4:8)]


foop <- as.data.frame(t(fmodsum ))
foop$estimates <- rownames(foop)

tidypar2 <- function(stanframe, param, id) {
    #take a stan model dataframe and create a tidy dataframe for a given parameter. takes a dataframe, a string with the parameter name, and a string describing the param ID (e.g. par = "Site" and id="SiteID")
    par <- stanframe %>% select(contains(param)) %>%
        mutate(estimate = rownames(stanframe)) %>%
        gather(key = key, value = value, contains("b_")) %>%
        mutate(id = str_extract(key, "[0-9]{1,}"))
    colnames(par) <- c("estimate","name", param, id)
    par[,4] <- as.integer(par[,4])
    return(par)
}

singledimparest <- foop %>%
    mutate(estimate=rownames(foop)) %>%
    select(estimate, beta, sigma_site, sigma_prov, sigma_clone, sigma_year, contains("kappa")) %>%
    rename(kappa1 = `kappa[1]`) %>%
    rename(kappa2 = `kappa[2]`)

siteest <- tidypar2(foop, "b_site", "SiteID")
provest <- tidypar2(foop, "b_prov", "ProvenanceID")
cloneest <- tidypar2(foop, "b_clone", "CloneID")
yearest <- tidypar2(foop, "b_year", "YearID")

ufdf <- ufdf %>% select(-sum_forcing, -recordID) %>% distinct()

#only include years, clones, sites, and provs actually in dataset
cloneestmerge <- left_join(ufdf, cloneest)
provestmerge <- left_join(ufdf, provest)
siteestmerge <- left_join(ufdf, siteest)
yearestmerge <- left_join(ufdf, yearest)

parestdf <- data.frame(ufdf,
                    estimate = cloneestmerge$estimate,
                    b_clone = cloneestmerge$b_clone,
                    b_prov=provestmerge$b_prov,
                    b_site = siteestmerge$b_site,
                    b_year = yearestmerge$b_year) %>%
    left_join(singledimparest)




# forcing <- seq(from = 175, to=500, length.out=50)
index <- unique(parestdf$index)
# forcingframe <- crossing(index, forcing)





ts <- parestdf %>%
    mutate(betas = b_clone + b_prov + b_site + b_year + beta) %>%
    left_join(clim) %>%
    mutate(t1 = logistic2(sum_forcing, betas, kappa1)) %>%
    mutate(t2 = logistic2(sum_forcing, betas, kappa2)) %>%
    mutate(p2 = pdflogistic(sum_forcing, betas, kappa1)) %>%
    mutate(p3 = pdflogistic(sum_forcing, betas, kappa2)) %>%
    mutate(f501 = kappa1/betas) %>%
    mutate(f502 = kappa2/betas)
    distinct()

# ts_nc <- parestdf %>%
#     mutate(betas = b_prov + b_site + b_year + beta) %>%
#     left_join(clim) %>%
#     mutate(t1 = logistic2(sum_forcing, betas, kappa1)) %>%
#     mutate(t2 = logistic2(sum_forcing, betas, kappa2)) %>%
#     mutate(p2 = pdflogistic(sum_forcing, betas, kappa1)) %>%
#     mutate(p3 = pdflogistic(sum_forcing, betas, kappa2)) %>%
#     distinct()

ts$indexest <- group_indices(ts, groupID, estimate)
#ts_nc$indexest <- group_indices(ts_nc, groupID, estimate)

ggplot(ts, aes(x=SPU_Name, y=f501, fill=estimate)) +
    geom_violin()

ggplot(filter(ts, estimate=="mean"), aes(x=sum_forcing, y=t1, group=indexest) ) +
    geom_line() +
    facet_grid(Site ~ SPU_Name) +
    geom_line(data=filter(ts_nc, estimate=="mean"), aes(x=sum_forcing, y=t1, group=indexest, color="nc"), alpha=0.2)

ggplot(filter(ts, estimate=="mean"), aes(x=sum_forcing, y=t2, group=indexest) ) +
    geom_line() +
    facet_grid(Site ~ SPU_Name) +
    geom_line(data=filter(ts_nc, estimate=="mean"), aes(x=sum_forcing, y=t2, group=indexest, color="nc"), alpha=0.2)

ggplot(filter(ts, estimate=="mean"), aes(x=DoY, y=t1, group=indexest) ) +
    geom_line() +
    facet_grid(Site ~ SPU_Name) +
    geom_line(data=filter(ts, estimate=="mean"), aes(x=DoY, y=t2, group=indexest, color="t2"), alpha=0.2)

ggplot(filter(ts, estimate=="mean"), aes(x=DoY, y=p2, group=indexest, color="p2") ) +
    geom_line() +
    facet_grid(Site ~ SPU_Name) +
    geom_line(data=filter(ts, estimate=="mean"), aes(x=DoY, y=p3, group=indexest, color="p3"), alpha=0.1)

# calculate state probabilities





