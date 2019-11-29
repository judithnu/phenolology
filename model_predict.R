# model predict
# predict fstart and fend from the model.

library(tidyverse)
library(gtools)
library(ggstance)
# read in phenology & model data ##############

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

# Length of forcing period ###########
forcinglengthdf <- rbind(fforce, mforce) %>%
    mutate(forcinglength = fend-fbegin)


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


phendf$Site <- factor(phendf$Site, ordered=TRUE, 
                      levels= rev(c("KettleRiver", "Vernon", "Kalamalka", "PRT", "Tolko", "PGTIS")))
doyperiod_hpd$Site <- factor(doyperiod_hpd$Site, ordered=TRUE, 
                             levels= rev(c("KettleRiver", "Vernon", "Kalamalka", "PRT", "Tolko", "PGTIS")))

doyperiod_hpd_wide <- doyperiod_hpd %>%
    select(-doylength) %>%
    pivot_longer(cols=c("dbegin", "dend"), names_to = "side", values_to = "day") %>%
    pivot_wider(names_from = c(hpd_end, side), values_from = day)

# periodlength_hpd <- full_join(doyperiod_hpd, hpd_long)


minmax <- filter(phendf, Phenophase_Derived==2) %>%
    summarise(minday=min(DoY), maxday=max(DoY))

minmaxbyyrsite <- filter(phendf, Phenophase_Derived==2) %>%
    group_by(Site, Year, Sex) %>%
    summarise(minday=min(DoY), maxday=max(DoY)) %>%
    filter(Site %in% c("KettleRiver", "Kalamalka", "PGTIS"))

#50 and 90% HPDI for DoY 
ggplot(doyperiod_hpd_wide, aes(ymin = hpd_low_dbegin, ymax = hpd_high_dend, 
                               x=Sex, color=Sex, size=1-intervalwidth)) +
    geom_linerange(alpha=0.5) +
    geom_linerange(data=doyperiod_hpd_wide, aes(ymin=hpd_high_dbegin, ymax=hpd_low_dend, 
                                                x=Sex, color=Sex, size=1-intervalwidth), alpha=0.5) +
    scale_color_viridis_d(end=0.8) +
    facet_grid(Site ~ Year) +
    geom_hline(yintercept=c(minmax$minday, minmax$maxday), color="gray") +
    theme_bw(base_size=16) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none",
          strip.text.x = element_text(size=10, angle=0),
          strip.text.y = element_text(size=8, angle=0),
          panel.border = element_blank(), 
          panel.grid=element_blank()) +
    guides(size=FALSE) +
    ylab("Day of Year")+
    geom_point(data=minmaxbyyrsite, aes(x=Sex, y = minday), inherit.aes = FALSE, pch="_", size=5, stroke=2) +
    geom_point(data=minmaxbyyrsite, aes(x=Sex, y = maxday), inherit.aes = FALSE, pch="_", size=5, stroke=2)


head(doyperiod_hpd_wide)


# Overlap analysis ##############

calc_intervals <- function(a, b) {
    assertthat::are_equal(nrow(a), nrow(b))
    ret <- data.frame(overlaps=rep(NA, nrow(a)))
    for (i in 1:nrow(a)) {
        a_begin <- a$dbegin[i]
        b_begin <- b$dbegin[i]
        a_end <- a$dend[i]
        b_end <- b$dend[i]
        if (a_begin > b_end || b_begin > a_end) {
            ret$overlaps[i] <- 0
        } else {
            lo <- max(c(a_begin, b_begin))
            hi <- min(c(a_end, b_end))
            ret$overlaps[i] <- (hi-lo)+1
        }
    }
    return(ret)
}

test_calc_intervals <- function() {
    frame_a <- data.frame(dbegin=c(133, 100, 100, 100, 90, 90),
                          dend=c(147, 110, 100, 105, 100, 100))
    
    frame_b <- data.frame(dbegin=c(135, 100, 100, 105, 110, 101),
                          dend=c(144, 110, 100, 110, 120, 120))
    
    expected <- calc_intervals(frame_a, frame_b)
    
    assertthat::assert_that(
        assertthat::are_equal(expected,
                              data.frame(overlaps=c(10, 11, 1, 1, 0, 0))))
}

test_calc_intervals()

fdoyperiod <- filter(doyperiod_hpd, Sex=="FEMALE")
mdoyperiod <- filter(doyperiod_hpd, Sex=="MALE")

overlap <- calc_intervals(fdoyperiod, mdoyperiod) #REAL slow
overlap2 <- calc_intervals(mdoyperiod, fdoyperiod)

fdoyperiod$overlaps <- overlap$overlaps
mdoyperiod$overlaps <- overlap$overlaps

doyperiod_overlap <- rbind(fdoyperiod, mdoyperiod) 

overlap_hpd <- doyperiod_overlap %>% 
    pivot_wider(names_from = hpd_end, values_from = c(dbegin, dend, doylength, overlaps)) %>%
    mutate(meanlength = (doylength_hpd_low + doylength_hpd_high)/2, 
           meanoverlap= (overlaps_hpd_low + overlaps_hpd_high)/2,
           meanbegin = (dbegin_hpd_low + dbegin_hpd_high)/2,
           meanend = (dend_hpd_low + dend_hpd_high)/2)

hist(doyperiod_overlap$overlaps)

ggplot(doyperiod_overlap, aes(x=doylength, y=overlaps, color=Site)) +
    geom_jitter() +
    facet_wrap(Sex ~ .)

ggplot(doyperiod_overlap, aes(x=dbegin, y=overlaps, color=Site, alpha=1-intervalwidth)) +
    geom_jitter() +
    facet_wrap(Sex ~ .)

ggplot(doyperiod_overlap, aes(x=dend, y=overlaps, color=Site)) +
    geom_jitter() +
    facet_wrap(Sex ~ .)

ggplot(overlap_hpd, aes(xmin=doylength_hpd_low, xmax=doylength_hpd_high, y=meanoverlap, 
                        color = Site, alpha=(1/intervalwidth)/4)) +
    geom_linerangeh() +
    geom_linerange(data=overlap_hpd, aes(ymin=overlaps_hpd_low, ymax=overlaps_hpd_high, x=meanlength, 
                                         color=Site, alpha=(1-intervalwidth)/4), inherit.aes = FALSE) +
    geom_point(aes(meanlength, meanoverlap, color=Site), inherit.aes = FALSE, pch=1) +
    facet_wrap("Sex") +
    scale_color_viridis_d(option="B", end=0.8) +
    theme_bw(base_size=16) +
    guides(alpha=FALSE) +
    ylab("Days of overlap") +
    xlab("Length of flowering period (days)")

ggplot(overlap_hpd, aes(xmin=dbegin_hpd_low, xmax=dbegin_hpd_high, y=meanoverlap, 
                        color = Site, alpha=(1/intervalwidth)/4)) +
    geom_linerangeh() +
    geom_linerange(data=overlap_hpd, aes(ymin=overlaps_hpd_low, ymax=overlaps_hpd_high, x=meanbegin, 
                                         color=Site, alpha=(1/intervalwidth)/4), inherit.aes = FALSE) +
    geom_point(aes(meanbegin, meanoverlap, color=Site), inherit.aes = FALSE, pch=1) +
    facet_wrap("Sex") +
    scale_color_viridis_d(option="B", end=0.8) +
    theme_bw(base_size=16) +
    guides(alpha=FALSE) +
    ylab("days of overlap") +
    xlab("Start of flowering")

ggplot(overlap_hpd, aes(x=meanbegin, y=meanoverlap, color=Site)) +
    geom_point(pch=1) +
    facet_wrap("Sex") +
    scale_color_viridis_d(option="B", end=0.8) +
    theme_bw(base_size=16) +
    theme(legend.position = "top") +
    ylab("Days of overlap") +
    xlab("Start of flowering")

lm(meanoverlap ~ meanbegin, overlap_hpd)
lm(meanoverlap ~ meanend, overlap_hpd)

ggplot(overlap_hpd, aes(x=meanend, y=meanoverlap, color=Site)) +
    geom_point(pch=1) +
    facet_wrap("Sex") +
    scale_color_viridis_d(option="B", end=0.8) +
    theme_bw(base_size=16) +
    ylab("Days of overlap") +
    xlab("End of flowering")

ggplot(overlap_hpd, aes(xmin=dend_hpd_low, xmax=dend_hpd_high, y=meanoverlap, 
                        color = Site, alpha=(1/intervalwidth)/4)) +
    geom_linerangeh() +
    geom_linerange(data=overlap_hpd, aes(ymin=overlaps_hpd_low, ymax=overlaps_hpd_high, x=meanend, 
                                         color=Site, alpha=(1/intervalwidth)/4), inherit.aes = FALSE) +
    geom_point(aes(meanbegin, meanoverlap, color=Site), inherit.aes = FALSE, pch=1) +
    facet_wrap("Sex") +
    scale_color_viridis_d(option="B", end=0.8) +
    theme_bw(base_size=16) +
    guides(alpha=FALSE) +
    ylab("days of overlap") +
    xlab("End of flowering")

#overlap between sites #########

calc_intervals <- function(a, b) {
    assertthat::are_equal(nrow(a), nrow(b))
    ret <- data.frame(overlaps=rep(NA, nrow(a)))
    for (i in 1:nrow(a)) {
        a_begin <- a$dbegin[i]
        b_begin <- b$dbegin[i]
        a_end <- a$dend[i]
        b_end <- b$dend[i]
        if (a_begin > b_end || b_begin > a_end) {
            ret$overlaps[i] <- 0
        } else {
            lo <- max(c(a_begin, b_begin))
            hi <- min(c(a_end, b_end))
            ret$overlaps[i] <- (hi-lo)+1
        }
    }
    return(ret)
}

overlapframer <- function(df) {
    doyperiod_meds <- df %>%
        group_by(Sex, siteyear, Site, Year) %>%
        summarise(dbegin = median(dbegin),
                  dend = median(dend),
                  doylength = median(doylength))
    
    kalamalkaf <- filter(doyperiod_meds, Site=="Kalamalka" & Sex =="FEMALE")
    kalamalkam <- filter(doyperiod_meds, Site=="Kalamalka" & Sex == "FEMALE")
    pgtisf <- filter(doyperiod_meds, Site=="PGTIS" & Sex =="FEMALE")
    pgtism <- filter(doyperiod_meds, Site =="PGTIS" & Sex =="MALE")
    krf <- filter(doyperiod_meds, Site=="KettleRiver" & Sex=="FEMALE")
    krm <- filter(doyperiod_meds, Site=="KettleRiver" & Sex=="MALE")
    
    kalamalkaf$kal_males <- calc_intervals(kalamalkaf, kalamalkam)$overlaps
    kalamalkaf$pg_males <- calc_intervals(kalamalkaf, pgtism)$overlaps
    kalamalkaf$kr_males <- calc_intervals(kalamalkaf, krm)$overlaps
    
    pgtisf$pg_males <- calc_intervals(pgtisf, pgtism)$overlaps
    pgtisf$kal_males <- calc_intervals(pgtisf, kalamalkam)$overlaps
    pgtisf$kr_males <- calc_intervals(pgtisf, krm)$overlaps
    
    krf$kr_males <- calc_intervals(krf, krm)$overlaps
    krf$pg_males <- calc_intervals(krf, pgtism)$overlaps
    krf$kal_males <- calc_intervals(krf, kalamalkam)$overlaps
    
    receptiveoverlap <- rbind(kalamalkaf, pgtisf) 
    receptiveoverlap <- rbind(receptiveoverlap, krf)
    
    roverlap <- receptiveoverlap %>% 
        pivot_longer(cols=contains("males"), names_to = "pollensource", values_to = "overlapdays")
    roverlap$Site <- factor(roverlap$Site, ordered=TRUE, 
                            levels= rev(c("KettleRiver", "Vernon", "Kalamalka", "PRT", "Tolko", "PGTIS")))
    roverlap$pollensource <- factor(roverlap$pollensource, ordered=TRUE, 
                                    levels= rev(c("kr_males", "kal_males", "pg_males")))
    return(roverlap)
}

roverlapreal <- overlapframer(doyperiod) %>%
    mutate(temp="real")
roverlapplus2 <- overlapframer(doyperiodcc) %>%
    mutate(temp="plus2")
roverlapplus5 <- overlapframer(doyperiodplus5) %>%
    mutate(temp="plus5")

roverlap <- rbind(roverlapreal, roverlapplus2)
roverlap <- rbind(roverlap, roverlapplus5)

roverlap_med <- roverlap %>%
    group_by(Site, pollensource, temp) %>%
    summarise(median(overlapdays))

ggplot(roverlap, aes(x=Year, y=overlapdays, color=temp)) +
    geom_point(pch=1) +
    facet_grid(Site ~ pollensource) +
    ggtitle("Pollen sources for 3 sites", subtitle = "") +
    theme_bw(base_size=16) +
    scale_color_viridis_d(option="C", end=0.8)

ggplot(roverlap, aes(x=0, y=overlapdays, fill=temp)) +
    geom_violin(draw_quantiles = 0.5) +
    facet_grid(Site ~ pollensource) +
    scale_fill_viridis_d(option="C", end=0.8) +
    theme_bw() +
    ggtitle("Receptivity overlap with pollen shed at 3 sites") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) #MEDIAN ONLY


ggplot(kalamalkaf, aes(x=doylength, y=pg_males, color="pg_males")) +
    geom_point(pch=1) +
    geom_point(aes(y=kr_males, color="kr_males"), pch=1) +
    ylab("Days of overlap") +
    xlab("Length of period") +
    scale_color_viridis_d(option="B", end=0.8) +
    ggtitle("Kalamalka female overlap with males at other sites")

ggplot(kalamalkaf, aes(x=dbegin, y=pg_males, color="pg_males")) +
    geom_point(pch=1) +
    geom_point(aes(y=kr_males, color="kr_males"), pch=1) +
    ylab("Days of overlap") +
    xlab("Beginning of flowering period") +
    scale_color_viridis_d(option="B", end=0.8) +
    ggtitle("Kalamalka female overlap with males at other sites")

ggplot(kalamalkaf, aes(x=dend, y=pg_males, color="pg_males")) +
    geom_point(pch=1) +
    geom_point(aes(y=kr_males, color="kr_males"), pch=1) +
    ylab("Days of overlap") +
    xlab("Beginning of flowering period") +
    scale_color_viridis_d(option="B", end=0.8) +
    ggtitle("Kalamalka female overlap with males at other sites")

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
