# match DoY to forcing units

#needs tpars from ppc.R
forcingtype = "scaled_ristos"

library(tidyverse)

gclim <- read.csv("data/all_clim_PCIC.csv", header=TRUE, stringsAsFactors=FALSE) %>%
    filter(forcing_type==forcingtype)

#tpars <- read.csv("transformed_parameters.csv", header=TRUE, stringsAsFactors=FALSE)

tpars_clim <- arrange(tpars, sum_forcing)

clim <- gclim
clim$siteyear <- paste(clim$Site, clim$Year, sep='')

# choose cold, hot, and medium years.

climsort <- filter(clim, DoY == 180) %>%
    arrange(sum_forcing)



# KettleRiver 2011 is cold, Kalamalka 1998 is hot, and PRT 2008 is middle of the road.

# index clim and tpars by Site-Year

KR2011 <- filter(clim, siteyear == "KettleRiver2011") %>%
    arrange(sum_forcing)
KAL1998 <- filter(clim, siteyear == "Kalamalka1998") %>%
    arrange(sum_forcing)

#THIS WORKS. This compares forcing units in tpars to sum_forcing in the climate dataset and finds the closest sum_forcing that's greater than a given tpars. Then it extracts they Day of Year that sum_forcing occured on. The additional "+ 1" that follows the findInterval call adds one to the index returned: findInterval by default returns the index of the left-hand (low) side of an interval, and we want the right-hand (high) side of the interval, which is the next index in the interval-vector from the indexes returned.
tpars_clim$DoY_KR2011 <- KR2011$DoY[findInterval(tpars_clim$sum_forcing,
                                                 KR2011$sum_forcing) + 1]
tpars_clim$DoY_KAL1998 <- KAL1998$DoY[findInterval(tpars_clim$sum_forcing,
                                            KAL1998$sum_forcing) + 1]

#
# ggplot(filter(tpars_clim, effect == "all"), aes(x=DoY_KR2011, color=Sex, linetype=param)) +
#     geom_density(alpha = 0.5) +
#     facet_grid(SPU_Name ~ .)
#
# ggplot(filter(tpars_clim, param %in% c("fstart", "fend")), aes(x=DoY_KAL1998, color=Sex, linetype=param)) +
#     geom_density(alpha = 0.5) +
#     facet_grid(SPU_Name ~ .)

tparsgraph <- gather(tpars_clim, key="WeatherRegime", value="DoY", DoY_KAL1998, DoY_KR2011)

# ggplot(filter(tparsgraph, effect=="all"), aes(x=DoY, fill=Sex)) +
#     geom_histogram(alpha=0.7, binwidth=1) +
#     scale_fill_viridis_d() +
#     facet_grid(WeatherRegime ~ .) +
#     theme_bw()

## Calculate phenological period length ##############################

tpars_wide <- select(tparsgraph, -param, -sum_forcing) %>%
    spread(key=side, value=DoY)

tpars_wide$length <- tpars_wide$end-tpars_wide$begin

ggplot(filter(tpars_wide, effect=="all"), aes(x=Sex, y=length, fill=WeatherRegime)) +
    geom_boxplot() +
 #   facet_grid(WeatherRegime ~ Sex) +
    scale_fill_viridis_d(option="C", end=0, begin=0.7, labels=c("hot year", "cold year")) +
    theme_bw(base_size=20) +
    ggtitle("Phenological period length in cold and hot years") +
    facet_wrap(Site ~ SPU_Name)




## prov diff DoY

# tparsgraph <- filter(tpars, str_detect(param, "half")) %>%
#     gather(key="WeatherRegime", value="DoY", DoY_KAL1998, DoY_KR2011)
#
#
# # subtract Provenance from all
#
# halfparams <- filter(tparsgraph, str_detect(param, "_")))
#
# ggplot(halfparams, aes(x=DoY, fill=Sex)) +
#     geom_histogram(position="dodge") +
#     scale_fill_viridis_d() +
#     facet_grid(rows=vars(Sex,param), cols=vars(WeatherRegime), scales="free_y")

# compare with and without provenance effect #############

phd1 <- filter(tparsgraph, effect == "all") %>%
    select(-sum_forcing, -effect)
phd2 <- filter(tparsgraph, effect == "no_prov") %>%
    rename(no_prov=param, DoY_noprov=DoY) %>%
    select(-sum_forcing, -effect)
provcompd <- full_join(phd1, phd2) %>%
    mutate(daydiff = DoY_noprov-DoY)

# ggplot(provcompd, aes(x=DoY, fill="prov", group=side, linetype=WeatherRegime)) +
#     geom_density() +
#     geom_density(aes(x=DoY_noprov, fill="no_prov", group=side, linetype=WeatherRegime), alpha=0.5) +
#     scale_fill_viridis_d(option="B", begin = .1, end=.9) +
#     facet_grid(SPU_Name ~ Sex) +
#     ggtitle("Provenance effect DoY") +
#     xlab("Day of Year") +
#     theme_bw(base_size=18) +
#     theme(strip.text.y = element_text(angle = 0)) +
#     theme(legend.position= "top")


hpd_lower = function(x, prob) rethinking::HPDI(x, prob)[1]
hpd_upper = function(x, prob) rethinking::HPDI(x, prob)[2]

# calculate HPDIs for DoY diffs PROV #############

full <- group_by(provcompd, side, Sex, SPU_Name, WeatherRegime) %>%
    summarize(x=hpd_lower(daydiff, prob=.99), xend=hpd_upper(daydiff, prob=.99)) %>%
    mutate(interval="0.99")

fifty <- group_by(provcompd, side, Sex, SPU_Name, WeatherRegime) %>%
    summarize(x=hpd_lower(daydiff, prob=.75), xend=hpd_upper(daydiff, prob=.75)) %>%
    mutate(interval="0.75")

hpd_df <- full_join(full, fifty) %>%
    gather(key=end, value=value, x, xend) %>%
    arrange(interval)

hpd_df <- mutate(hpd_df, y = case_when(WeatherRegime=="DoY_KAL1998" ~ 0,
                               WeatherRegime=="DoY_KR2011" ~ 1))

# How to do an interval plot

trans1 <- filter(hpd_df, side=="begin")
ggplot(filter(trans1, interval=="0.99"), aes(x=value, y=y, color=WeatherRegime, size=1)) +
    geom_line() +
    geom_line(data = filter(trans1, interval=="0.75"), aes(x=value, y=y, size=1.5)) +
    facet_grid(SPU_Name ~ Sex) +
    theme_bw(base_size=16) +
    scale_color_viridis_d(option="C", end=0, begin=0.7, labels=c("hot year", "cold year")) +
    ylim(c(-1, 2)) +
    geom_vline(xintercept=0, alpha=0.7) +
    ggtitle("Provenance Effects in Days in Cold and Hot Years", subtitle = "Start Day. 75 and 95% HPDI. Negative is days earlier, positive is days later.") +
    theme(strip.text.y = element_text(angle = 0), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position="top") +
    xlab("Day Differences") +
    ylab("") +
    guides(size=FALSE, colour=guide_legend(override.aes = list(size=3)))+
    xlim(c(-20,20))

trans2 <- filter(hpd_df, side=="end")
ggplot(filter(trans2, interval=="0.99"), aes(x=value, y=y, color=WeatherRegime, size=1)) +
    geom_line() +
    geom_line(data = filter(trans1, interval=="0.75"), aes(x=value, y=y, size=1.5)) +
    facet_grid(SPU_Name ~ Sex) +
    theme_bw(base_size=16) +
    scale_color_viridis_d(option="C", end=0, begin=0.7, labels=c("hot year", "cold year")) +
    ylim(c(-1, 2)) +
    geom_vline(xintercept=0, alpha=0.7) +
    ggtitle("Provenance Effects in Days in Cold and Hot Years", subtitle = "End Day. 75 and 95% HPDI. Negative is days earlier, positive is days later.") +
    theme(strip.text.y = element_text(angle = 0), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position="top") +
    xlab("Day Differences") +
    ylab("") +
    guides(size=FALSE, colour=guide_legend(override.aes = list(size=3)))+
    xlim(c(-20,20))

# calculate HPDIs for DoY diffs SITE #############

shd1 <- filter(tparsgraph, effect == "all") %>%
    select(-sum_forcing, -effect)
shd2 <- filter(tparsgraph, effect == "no_site") %>%
    rename(no_site=param, DoY_nosite=DoY) %>%
    select(-sum_forcing, -effect)
sitecompd <- full_join(shd1, shd2) %>%
    mutate(daydiff = DoY_nosite-DoY)

full <- group_by(sitecompd, side, Sex, Site, WeatherRegime) %>%
    summarize(x=hpd_lower(daydiff, prob=.99), xend=hpd_upper(daydiff, prob=.99)) %>%
    mutate(interval="0.99")

fifty <- group_by(sitecompd, side, Sex, Site, WeatherRegime) %>%
    summarize(x=hpd_lower(daydiff, prob=.75), xend=hpd_upper(daydiff, prob=.75)) %>%
    mutate(interval="0.75")

hpd_df <- full_join(full, fifty) %>%
    gather(key=end, value=value, x, xend) %>%
    arrange(interval)


hpd_df <- mutate(hpd_df, y = case_when(WeatherRegime=="DoY_KAL1998" ~ 0,
                                       WeatherRegime=="DoY_KR2011" ~ 1))

trans1 <- filter(hpd_df, side=="begin")
ggplot(filter(trans1, interval=="0.99"), aes(x=value, y=y, color=WeatherRegime, size=1)) +
    geom_line() +
    geom_line(data = filter(trans1, interval=="0.75"), aes(x=value, y=y, size=1.5)) +
    facet_grid(Site ~ Sex) +
    theme_bw(base_size=16) +
    scale_color_viridis_d(option="C", end=0, begin=0.7, labels=c("hot year", "cold year")) +
    ylim(c(-1, 2)) +
    geom_vline(xintercept=0, alpha=0.7) +
    ggtitle("Site Effects in Days in Cold and Hot Years", subtitle = "Start Day. 75 and 95% HPDI. Negative is days earlier, positive is days later.") +
    theme(strip.text.y = element_text(angle = 0), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position="top") +
    xlab("Day Differences") +
    ylab("") +
    guides(size=FALSE, colour=guide_legend(override.aes = list(size=3)))+
    xlim(c(-20,20))

trans2 <- filter(hpd_df, side=="end")
ggplot(filter(trans2, interval=="0.99"), aes(x=value, y=y, color=WeatherRegime, size=1)) +
    geom_line() +
    geom_line(data = filter(trans1, interval=="0.75"), aes(x=value, y=y, size=1.5)) +
    facet_grid(SPU_Name ~ Sex) +
    theme_bw(base_size=16) +
    scale_color_viridis_d(option="C", end=0, begin=0.7, labels=c("hot year", "cold year")) +
    ylim(c(-1, 2)) +
    geom_vline(xintercept=0, alpha=0.7) +
    ggtitle("Provenance Effects in Days in Cold and Hot Years", subtitle = "End Day. 75 and 95% HPDI. Negative is days earlier, positive is days later.") +
    theme(strip.text.y = element_text(angle = 0), axis.text.y = element_blank(), axis.ticks = element_blank(), legend.position="top") +
    xlab("Day Differences") +
    ylab("") +
    guides(size=FALSE, colour=guide_legend(override.aes = list(size=3)))+
    xlim(c(-20,20))

# Graph DoY start and finish ##################

doyprov <- filter(tparsgraph, effect %in% c("all", "no_prov"), Sex == "FEMALE")
ggplot(doyprov, aes(x=DoY, color=effect, linetype=side)) +
    geom_freqpoly(size=1.1, alpha=0.8) +
    stat_ecdf(data=bdatd, aes(x=DoY_obs), inherit.aes=FALSE) +
    scale_color_viridis_d(option="B", end=0.9) +
    facet_grid(SPU_Name ~ WeatherRegime, scales = "free_y") +
    ggtitle("FEMALE Start and end by provenance", subtitle = "in a hot and cold year") +
    xlab("Forcing units") +
    theme_bw(base_size=18) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(legend.position= "top") +
    guides(size=FALSE, colour=guide_legend(override.aes = list(size=3)))

doysite <- filter(tparsgraph, effect %in% c("all", "no_site"), Sex == "FEMALE")
ggplot(doysite, aes(x=DoY, color=effect, linetype=side)) +
    geom_freqpoly(size=1.05, alpha=0.8) +
    #stat_ecdf(data=bdatd, aes(x=DoY_obs), inherit.aes=FALSE) +
    scale_color_viridis_d(option="B", end=.9) +
    facet_grid(Site ~ WeatherRegime, scales = "free_y") +
    ggtitle("FEMALE Start and end by site", subtitle = "in a hot and cold year") +
    xlab("Forcing units") +
    theme_bw(base_size=18) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(legend.position= "top") +
    guides(size=FALSE, colour=guide_legend(override.aes = list(size=3)))

doyprovcold <- filter(tparsgraph, effect %in% c("all", "no_prov"), WeatherRegime== "DoY_KR2011")
ggplot(doyprovcold, aes(x=DoY, color=effect, linetype=side)) +
    geom_freqpoly(alpha=0.8) +
   # stat_ecdf(data=bdatd, aes(x=DoY_obs), inherit.aes=FALSE) +
    scale_color_viridis_d(option="B", end=0.5) +
    facet_grid(SPU_Name ~ Sex, scales = "free_y") +
    ggtitle("Start and end by provenance in a cold year") +
    xlab("Day of Year") +
    theme_bw(base_size=18) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(legend.position= "top")

doysitecold <- filter(tparsgraph, effect %in% c("all", "no_site"), WeatherRegime== "DoY_KR2011")
ggplot(doysitecold, aes(x=DoY, color=effect, linetype=side)) +
    geom_freqpoly(alpha=0.8) +
    # stat_ecdf(data=bdatd, aes(x=DoY_obs), inherit.aes=FALSE) +
    scale_color_viridis_d(option="B", end=0.5) +
    facet_grid(SPU_Name ~ Sex, scales = "free_y") +
    ggtitle("Start and end by site in a cold year") +
    xlab("Day of Year") +
    theme_bw(base_size=18) +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(legend.position= "top")

timing <- filter(tparsgraph, effect == "all")
ggplot(timing, aes(x=DoY, color=WeatherRegime, linetype=side)) +
    geom_freqpoly(bins=40, size=1.5, alpha=0.8) +
    scale_color_viridis_d(option="C", end=0, begin=0.7, labels=c("hot year", "cold year")) +
    facet_grid(Sex ~ .) +
    theme_bw(base_size=18) +
    theme(legend.position= "top") +
    ggtitle("Start and end of flowering in a cold and hot year")

