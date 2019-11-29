# model predict
# predict fstart and fend from the model.

library(tidyverse)
library(gtools)
library(ggstance)
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
    sample_frac(0.25)
fmod$iter <- 1:nrow(fmod)
rm(fmod_raw)

mmod <- mmod_raw %>%
    select(contains("kappa"), beta, contains("mean"), contains("sigma")) %>%
    sample_frac(0.25)
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


forcinglengthdf <- rbind(fforce, mforce) %>%
    mutate(forcinglength = fend-fbegin)

# Length of forcing period ###########3
forcingperiod <- forcinglengthdf %>%
    select(-forcinglength) %>%
    pivot_longer(cols=starts_with("f"), names_to = "side", values_to = "sum_forcing")

# Forcing period hpdi and graphed #############
hpd_lower = function(x, prob) rethinking::HPDI(x, prob)[1]
hpd_upper = function(x, prob) rethinking::HPDI(x, prob)[2]

hpd <- function(df, interval) {
    df <- df %>%
        group_by(Sex, side) %>%
        summarise(median = median(sum_forcing), 
                  hpd_low = hpd_lower(sum_forcing, interval), 
                  hpd_high = hpd_upper(sum_forcing, interval)) #%>%
        #pivot_longer(starts_with("hpd"), names_to = "intervalside", values_to = "sum_forcing")
    df$intervalwidth <- interval
    return(df)
}


hpd50 <- hpd(forcingperiod, 0.5)
hpd90 <- hpd(forcingperiod, 0.9)
hpddf <- rbind(hpd50, hpd90) %>%
    mutate(y=case_when(side == "fbegin" ~ 0,
                       side == "fend" ~ 1))
hpd_length # calculate forcing length hpd.


#Forcing accumulation across the flowering period for males and females with model predictions for start and end
ggplot(hpddf, aes(x=median, y=y, color=Sex, group=side)) +
    geom_linerangeh(aes( xmin=hpd_low, xmax=hpd_high, y=y, size=1-intervalwidth), alpha=0.5 )+
    #geom_point(pch=16, size=3) +
    geom_vline(aes(xintercept=median, color=Sex, group=side)) +
    facet_grid(rows = vars(Sex)) +
    #geom_point(aes(x=median, y=y), size=5, alpha=1) +
    theme_bw(base_size = 18) +
    theme(strip.text.y = element_text(angle = 0), legend.position = "none") +
    scale_color_viridis_d(end=0.8) +
    stat_ecdf(data=filter(phendf, Phenophase_Derived==2), aes(x=sum_forcing), inherit.aes = FALSE) +
    xlab("accumulated forcing") +
    ylab("") 

# density plot showing the same thing as the plot above
ggplot(forcingperiod, aes(x=sum_forcing, y=..scaled..,fill=Sex, linetype=side)) +
    geom_density(alpha=0.5) +
    #geom_point(pch=16, size=3) +
    #geom_vline(data=hpddf, aes(xintercept=median, color=Sex, group=side)) +
    #facet_grid(rows = vars(Sex)) +
    #geom_point(aes(x=median, y=y), size=5, alpha=1) +
    theme_bw(base_size = 18) +
    theme(strip.text.y = element_text(angle = 0), legend.position = "none") +
    scale_color_viridis_d(end=0.8) +
    scale_fill_viridis_d(end=0.8) +
    stat_ecdf(data=filter(phendf, Phenophase_Derived==2), aes(x=sum_forcing, color=Sex), inherit.aes = FALSE) +
    xlab("accumulated forcing") +
    ylab("") +
    xlim(c(5,25))


# read in climate data ###############

clim <- read.csv("data/all_clim_PCIC.csv", header=TRUE, stringsAsFactors=FALSE) %>%
    filter(forcing_type==forcingtype) %>%
    filter(Site %in% c("KettleRiver", "Kalamalka", "PGTIS"))
clim$siteyear <- paste(clim$Site, clim$Year, sep='')

climchange <- read.csv("data/all_clim_PCIC_climchange.csv", header=TRUE, stringsAsFactors = FALSE) %>%
    filter(forcing_type==forcingtype) %>%
    filter(Site %in% c("KettleRiver", "Kalamalka", "PGTIS"))
climchange$siteyear <- paste(climchange$Site, climchange$Year, sep = '')

climplus5 <- read.csv("data/all_clim_PCIC_plus5.csv", header=TRUE, stringsAsFactors = FALSE) %>%
    filter(forcing_type==forcingtype) %>%
    filter(Site %in% c("KettleRiver", "Kalamalka", "PGTIS"))
climplus5$siteyear <- paste(climchange$Site, climplus5$Year, sep = '')


# ryears <- sample(unique(clim$siteyear), 20)
# clim <- filter(clim, siteyear %in% ryears)

# CLIMATE PREDICTIONS/MATCHING ###########
splitclim <- split(clim, clim$siteyear)
splitclimcc <- split(climchange, climchange$siteyear)
splitclimplus5 <- split(climplus5, climplus5$siteyear)

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
                     Site = x$Site[1],
                     Year = x$Year[1],
                     Sex=y$Sex)
    return(df)
}

# keep only forcing values in the 90% HPDI
# forcingperiod90 <- forcingperiod %>%
#     left_join(hpd90) %>%
#     filter(sum_forcing > hpd_low && sum_forcing < hpd_high) %>%
#     select(-contains("hpd"), -contains("interval"), -median) %>%
#     pivot_wider(names_from = side, values_from = sum_forcing)

#doyperiod <- purrr::map_df(splitclim, ifinder, forcinglengthdf) %>%
doyperiod <- purrr::map_df(splitclim, ifinder, forcinglengthdf) %>%
    mutate(doylength = dend-dbegin)
doyperiodcc <- purrr::map_df(splitclimcc, ifinder, forcinglengthdf) %>%
    mutate(doylength = dend-dbegin)
doyperiodplus5 <- purrr::map_df(splitclimplus5, ifinder, forcinglengthdf) %>%
    mutate(doylength = dend-dbegin)

assertthat::assert_that(nrow(doyperiod)==nrow(forcinglengthdf) * length(splitclim))


#periodlength <- full_join(doyperiod, forcinglengthdf)
#periodlength <- full_join(doyperiod, forcingperiod90)

#assertthat::assert_that(nrow(doyperiod)==nrow(periodlength))

doyperiod$Site <- factor(doyperiod$Site, ordered=TRUE, levels = c("KettleRiver", "Kalamalka", "PRT", "Sorrento", "PGTIS"))
doyperiodcc$Site <- factor(doyperiodcc$Site, ordered=TRUE, levels = c("KettleRiver", "Kalamalka", "PRT", "Sorrento", "PGTIS"))
doyperiodplus5$Site <- factor(doyperiodplus5$Site, ordered=TRUE, levels = c("KettleRiver", "Kalamalka", "PRT", "Sorrento", "PGTIS"))

doyperiod_hpd_length <- doyperiod %>%
    group_by(Sex, siteyear) %>%
    summarise(median=median(doylength),
              daylength_low_50 = hpd_lower(doylength, 0.5), daylength_high_50=hpd_upper(doylength, 0.5),
              daylength_low_90 = hpd_lower(doylength, 0.9), daylength_high_90=hpd_upper(doylength, 0.9),
              Site=unique(Site), Year=unique(Year))

ggplot(doyperiod_hpd_length, aes(x=Sex, y=median, color=Sex)) +
    geom_point() +
    geom_linerange(aes(ymin=daylength_low_50, ymax=daylength_high_50), size=1) +
    geom_linerange(aes(ymin=daylength_low_90, ymax=daylength_high_90), size=0.5) +
    facet_grid(Site ~ Year) +
    scale_color_viridis_d(end=0.8) +
    theme_bw(base_size = 15) +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylab("Days") 

doyperiod %>% 
    group_by(Site, Sex) %>%
    summarise(medianlength = median(doylength))


ggplot(doyperiod, aes(x=as.factor(Year), y=dbegin)) +
    geom_bin2d(binwidth=c(1,1)) +
    facet_grid(rows=vars(Site), cols= vars(Sex)) +
    scale_fill_continuous(type="viridis", option="B")

ggplot(doyperiod, aes(x=siteyear, y=dbegin, color=Site)) +
    geom_violin() +
    facet_wrap("Sex")
# add hline or grey poly for actual data start. maybe only use the 50% hpdi? or 90? or calculate

# Only HPDIs ################

ifinder_hpd <-  function(x, y) {
    beginindex <- findInterval(y$fbegin, x$sum_forcing)
    endindex <- findInterval(y$fend, x$sum_forcing)
    dbegin <- x$DoY[beginindex]
    dend <- x$DoY[endindex]
    df <- data.frame(dbegin=dbegin,
                     dend=dend,
                     siteyear=x$siteyear[1],
                     Site = x$Site[1],
                     Year = x$Year[1],
                     Sex=y$Sex, 
                     intervalwidth=y$intervalwidth, 
                     hpd_end = y$hpd_end)
    return(df)
}

hpd_long <- hpddf %>%
    pivot_longer(cols=contains("hpd"), names_to = "hpd_end", values_to = "sum_forcing") %>%
    select(-y, -median) %>%
    pivot_wider(names_from = side, values_from = sum_forcing) 

doyperiod_hpd <- purrr::map_df(splitclim, ifinder_hpd, hpd_long) %>%
    mutate(doylength = dend-dbegin) 

doyperiod_hpd_wide <- doyperiod_hpd %>%
    select(-doylength) %>%
    pivot_longer(cols=c("dbegin", "dend"), names_to = "side", values_to = "day") %>%
    pivot_wider(names_from = c(hpd_end, side), values_from = day)

# periodlength_hpd <- full_join(doyperiod_hpd, hpd_long)

doyperiod_hpd$Site <- factor(doyperiod_hpd$Site, ordered=TRUE, levels = c("KettleRiver", "Kalamalka", "PRT", "Sorrento", "PGTIS"))

ggplot(filter(doyperiod_hpd_wide, Year < 2001), aes(xmin=hpd_low, xmax=hpd_high, y=0, color=Sex, size=1-intervalwidth)) +
    geom_linerangeh(alpha=0.5) +
    scale_color_viridis_d(end=0.8) +
    facet_grid(siteyear ~ .)

minmax <- filter(phendf, Phenophase_Derived==2) %>%
    summarise(minday=min(DoY), maxday=max(DoY))

#50 and 90% HPDI for DoY 
ggplot(doyperiod_hpd_wide, aes(ymin=hpd_low_dbegin, ymax=hpd_high_dend, 
                                                    x=Sex, color=Sex, size=1-intervalwidth)) +
    geom_linerange(alpha=0.5) +
    geom_linerange(aes(ymin=hpd_high_dbegin, ymax=hpd_low_dend, 
                       x=Sex, color=Sex, size=1-intervalwidth), alpha=0.5) +
    scale_color_viridis_d(end=0.8) +
    facet_grid(Site ~ Year) +
    geom_hline(yintercept=c(minmax$minday, minmax$maxday)) +
    theme_bw(base_size=18) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none",
         strip.text.x = element_text(size=10, angle=0),
         strip.text.y = element_text(size=8, angle=0),
          panel.border = element_blank()) +
    guides(size=FALSE) 
    
    




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
