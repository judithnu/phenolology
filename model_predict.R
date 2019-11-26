# model predict
# predict fstart and fend from the model.

library(tidyverse)
library(gtools)
# read in climate and model data ##############

source('phenology_functions.R')

phendf <- read_data(slim = FALSE)

fmod_raw <- readRDS("2019-10-28phenologyFEMALE.rds") %>%
    as.data.frame() %>%
    select(-contains("state_rep"))

mmod_raw <- readRDS("2019-10-28phenologyMALE.rds") %>%
    as.data.frame() %>%
    select(-contains("state_rep"))

# GLOBALS ###########
forcingtype="scaled_ristos"

# MODELS #############
fmod <- fmod_raw %>%
    select(contains("kappa"), beta, contains("mean"), contains("sigma")) %>%
    sample_frac(0.5)
fmod$iter <- 1:nrow(fmod)
rm(fmod_raw)

mmod <- mmod_raw %>%
    select(contains("kappa"), beta, contains("mean"), contains("sigma")) %>%
    sample_frac(0.5)
mmod$iter <- 1:nrow(mmod)
rm(mmod_raw)

# Simulate effects from model ##########

# rnorm vectorizes ok, uses recycling so be VERY CAREFUL, e.g.
# > rnorm(4, mean=c(0, 60000), sd=c(1, 1))
# [1] 1.033847e-01 6.000000e+04 1.433197e+00 6.000185e+04

# model is the stan model with effects b_[effect]_mean and sigma_[effect]
# predict the amount of forcing required to start and end the flowering period
# for each iteration from the model, simulate 30 effects each for site, prov, clone, and year
predict_forcing <- function(model, samples, sex) {
    site <- rnorm(samples*nrow(model), model$b_site_mean, model$sigma_site)
    prov <- rnorm(samples*nrow(model), model$b_prov_mean, model$sigma_prov)
    clone <- rnorm(samples*nrow(model), model$b_clone_mean, model$sigma_clone)
    year <- rnorm(samples*nrow(model), model$b_year_mean, model$sigma_year)
    
    # calculate fstart and fend
    
    fbegin <- (logit(0.2) + model$`kappa[1]`)/(model$beta + site + prov + clone + year)
    fend <-  (logit(0.8) + model$`kappa[2]`)/(model$beta + site + prov + clone + year)
    
    beginend <- data.frame(iter=rep(model$iter, samples), draw=sort(rep(1:samples, nrow(model))), fbegin, fend, Sex=sex )
    
    return(beginend)
}

fforce <- predict_forcing(fmod, 30, "FEMALE")
mforce <- predict_forcing(mmod, 30, "MALE")

# # calculate 50% hpdi here
# 
# hpd_lower = function(x, prob) rethinking::HPDI(x, prob)[1]
# hpd_upper = function(x, prob) rethinking::HPDI(x, prob)[2]
# 


forcinglengthdf <- rbind(fforce, mforce) %>%
    mutate(forcinglength = fend-fbegin)

forcingperiod <- forcinglengthdf %>%
    select(-forcinglength) %>%
    pivot_longer(cols=starts_with("f"), names_to = "side", values_to = "sum_forcing")

# forcing_hpdi <- forcingperiod %>%
#     group_by(Sex, side) %>%
#     mutate(hpd_low = hpd_lower(sum_forcing, 0.50), hpd_high = hpd_upper(sum_forcing, 0.50)) %>%
#     filter(sum_forcing > hpd_low & sum_forcing < hpd_high)
# 
clim <- read.csv("data/all_clim_PCIC.csv", header=TRUE, stringsAsFactors=FALSE) %>%
    filter(forcing_type==forcingtype) %>%
    filter(!Site %in% c("Vernon", "Tolko"))
clim$siteyear <- paste(clim$Site, clim$Year, sep='')
# ryears <- sample(unique(clim$siteyear), 20)
# clim <- filter(clim, siteyear %in% ryears)

# CLIMATE PREDICTIONS ###########
splitclim <- split(clim, clim$siteyear)

# find the DoY that a given sum of forcing accumulated on
#x is a dataframe of climate with sum_forcing and the DoY that occured on and
#y is a dataframe with forcing data related to phenology. that has columns for iteration, sex, and draws. forcing cols must be fbegin and fend
ifinder <-  function(x, y) {
    beginindex <- findInterval(y$fbegin, x$sum_forcing)
    endindex <- findInterval(y$fend, x$sum_forcing)
    dbegin <- x$DoY[beginindex]
    dend <- x$DoY[endindex]
    df <- data.frame(iter=y$iter,
                     draw=y$draw,
                     dbegin=dbegin,
                     dend=dend,
                     siteyear=x$siteyear[1],
                     Sex=y$Sex)
    return(df)
}

doyperiod <- purrr::map_df(splitclim, ifinder, forcinglengthdf) %>%
    mutate(doylength = dend-dbegin)

assertthat::assert_that(nrow(doyperiod)==nrow(forcinglengthdf) * length(splitclim))

periodlength <- full_join(doyperiod, forcinglengthdf)

assertthat::assert_that(nrow(doyperiod)==nrow(periodlength))

# mean and standard deviation of period length in days at each forcing length
periodlengthsummary <- periodlength %>%
    group_by(forcinglength, Sex) %>%
    summarise(mean=mean(doylength), sd=sd(doylength))

# x climate timeseries, 7000 iterations, 30 draws, then median of start, end, length per iteration 
# (so median of the x draws per iteration for a given sex and timeseries)
periodlength_drawsummary <- periodlength %>%
    group_by(iter, Sex, siteyear) %>%
    summarise(dbeginm=median(dbegin), dendm=median(dend), doylengthm=median(doylength), 
              fbeginm=median(fbegin), fendm=median(fend), forcinglengthm=median(forcinglength)) # this is extremely slow

#Actually, medians should be calculated separately for forcing and day terms - 
#forcing shouldn't be grouped by siteyear! (thought it shouldn't technically matter here)

# x climate timeseries, 7000 iterations, 30 draws
ggplot(periodlengthsummary, aes(x=forcinglength, y=mean)) +
    geom_point(pch=1, alpha=0.3) +
    facet_grid(Sex ~ .) +
    ylab("Mean DAYS of flowering") +
    theme_bw() +
    ggtitle("Length of flowering period in days vs forcing units")

# x climate timeseries, 7000 iterations, 30 draws, then median of start, end, length per iteration 
# (so median of the x draws per iteration for a given sex and timeseries)
ggplot(periodlength_drawsummary, aes(x=forcinglengthm, y=doylengthm)) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    theme_bw(base_size=18) +
    scale_fill_viridis_c(option="B", begin=0.2) +
    facet_grid(Sex ~ .) +
    ylab("Days spent flowering") +
    xlab("Forcing units") +
    ggtitle("Forcing units required to advance through flowering period", subtitle = "and how many days it takes to do so")

# this is with medians!
ggplot(periodlength_drawsummary, aes(x=forcinglengthm, y=doylengthm)) +
    geom_bin2d(binwidth=c(.1, 1)) +
    theme_bw(base_size=18) +
    scale_fill_viridis_c(option="cividis") +
    facet_grid(Sex ~ .) +
    ylab("Days spent flowering") +
    xlab("Forcing units") +
    ggtitle("Forcing units required to advance through flowering period", subtitle = "and how many days it takes to do so")

# # this is with all the draws!, but shows something slightly different because of grouping in drawsummary calc
ggplot(periodlength, aes(x=forcinglength, y=doylength)) +
    geom_bin2d(binwidth=c(.1, 1)) +
    theme_bw(base_size=18) +
    scale_fill_viridis_c(option="cividis") +
    facet_grid(Sex ~ .) +
    ylab("Days spent flowering") +
    xlab("Forcing units") +
    ggtitle("Forcing units required to advance through flowering period", subtitle = "and how many days it takes to do so")

# maybe plot this as HDPIntervals in 2 dimensions? how would that even work?



fit <- lm(doylength ~ forcinglength, periodlength)
summary(fit)

ggplot(periodlength, aes(x=forcinglength, y=doylength)) +
    geom_point(pch=1, alpha=0.1)

ggplot(periodlength, aes(x=forcinglength, y=doylength)) +
    geom_violin() +
    facet_wrap("Sex")

ggplot(forcingfemale, aes(x=forcing, y= ..scaled.., fill=pp, group=iter)) +
    geom_density(alpha=0.4) +
    stat_ecdf(data=filter(phendf, Phenophase_Derived==2), aes(x=sum_forcing), inherit.aes = FALSE, color="yellow") +
    xlim(c(5, 36)) +
    facet_grid(pp ~ .) +
    scale_fill_viridis_d(option="A", end=0.8) +
    theme_bw(base_size = 18) +
    theme(legend.position="none") +
    ggtitle("Forcing requirements for start and end of phenological period")
