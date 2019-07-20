# match DoY to forcing units

#needs tpars from ppc.R

library(tidyverse)

gclim <- read.csv("data/all_clim_PCIC.csv", header=TRUE, stringsAsFactors=FALSE) %>%
    filter(forcing_type=="scaled_ristos")

#tpars <- read.csv("transformed_parameters.csv", header=TRUE, stringsAsFactors=FALSE)

tpars <- arrange(tpars, sum_forcing)

clim <- gclim
clim$siteyear <- paste(clim$Site, clim$Year, sep='')

# choose cold, hot, and medium years.

climsort <- filter(clim, DoY == 180) %>%
    arrange(sum_forcing)

ggplot(clim, aes(x=DoY, y=sum_forcing, color=Site, group=siteyear)) +
    geom_line(alpha=0.5, size=1.5) +
    theme(legend.position = "bottom") +
    scale_color_viridis_d() +
    theme_bw() +
    ggtitle("Forcing accumulation at all sites between 1997 and 2012")+
    theme(text=element_text(size=20))


# KettleRiver 2011 is cold, Kalamalka 1998 is hot, and PRT 2008 is middle of the road.

# index clim and tpars by Site-Year

KR2011 <- filter(clim, siteyear == "KettleRiver2011") %>%
    arrange(sum_forcing)
KAL1998 <- filter(clim, siteyear == "Kalamalka1998") %>%
    arrange(sum_forcing)

#THIS WORKS. This compares forcing units in tpars to sum_forcing in the climate dataset and finds the closest sum_forcing that's greater than a given tpars. Then it extracts they Day of Year that sum_forcing occured on. The additional "+ 1" that follows the findInterval call adds one to the index returned: findInterval by default returns the index of the left-hand (low) side of an interval, and we want the right-hand (high) side of the interval, which is the next index in the interval-vector from the indexes returned.
tpars$DoY_KR2011 <- KR2011$DoY[findInterval(tpars$sum_forcing,
                                                 KR2011$sum_forcing) + 1]
tpars$DoY_KAL1998 <- KAL1998$DoY[findInterval(tpars$sum_forcing,
                                            KAL1998$sum_forcing) + 1]


ggplot(filter(tpars, param %in% c("fstart", "fend")), aes(x=DoY_KR2011, color=Sex, linetype=param)) +
    geom_density(alpha = 0.5) +
    facet_grid(SPU_Name ~ .)

ggplot(filter(tpars, param %in% c("fstart", "fend")), aes(x=DoY_KAL1998, color=Sex, linetype=param)) +
    geom_density(alpha = 0.5) +
    facet_grid(SPU_Name ~ .)

#tparsgraph <- filter(tpars, param %in% c("fstart", "fend")) %>%
    gather(key="WeatherRegime", value="DoY", DoY_KAL1998, DoY_KR2011)

ggplot(tparsgraph, aes(x=DoY, fill=Sex, linetype=param)) +
    geom_histogram(alpha=0.5, binwidth=1) +
    scale_fill_viridis_d(end = 0.9) +
    facet_grid(WeatherRegime ~ .) +
    theme_bw()

ggplot(tparsgraph, aes(x=DoY, y=IndSexGroup, color=Sex)) +
    geom_line(alpha=0.3) +
    facet_grid(WeatherRegime ~ .) +
    theme_bw()

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

# compare with and without provenance effect

provnoprov <- filter(tpars, param %in% c("fstart", "fstart_noprov")) %>%
    arrange(IndSexGroup, param, iter)

alleffects <- filter(provnoprov, param=="fstart")
noprov <- filter(provnoprov, param=="fstart_noprov")

provdiffs_hot <- noprov$DoY_KAL1998 - alleffects$DoY_KAL1998
provdiffs_cold <- noprov$DoY_KR2011 - alleffects$DoY_KR2011

provdiffs <- data.frame(alleffects[,c(1:12)], provdiffs_cold, provdiffs_hot)

HPDI(provdiffs_hot)
HPDI(provdiffs_cold)

ggplot(provdiffs, aes(x=provdiffs_cold, fill="cold")) +
    geom_histogram(alpha=.5) +
    geom_histogram(aes(x=provdiffs_hot, fill="hot"), alpha=0.5) +
    scale_fill_viridis_d(option="C") +
    facet_grid(SPU_Name ~ Sex, scales="free_y")
   # geom_line(data = data.frame(x=c(1,5), y=c(1,1)), aes(x=x, y=y), size=2) +
    #xlim(c(-7,15))

# what about larger site effects
sitenosite <- filter(tpars, param %in% c("fhalf1", "fhalf1_nosite")) %>%
    arrange(IndSexGroup, param, iter)

alleffects <- filter(sitenosite, param=="fhalf1")
nosite <- filter(sitenosite, param=="fhalf1_nosite")

sitediffs_hot <- nosite$DoY_KAL1998 - alleffects$DoY_KAL1998
sitediffs_cold <- nosite$DoY_KR2011 - alleffects$DoY_KR2011

sitediffs <- data.frame(alleffects[,c(1:12)], sitediffs_cold, sitediffs_hot)

HPDI(sitediffs_hot)
HPDI(sitediffs_cold)

ggplot(sitediffs, aes(x=sitediffs_cold, fill="cold")) +
    geom_histogram(alpha=.5) +
    geom_histogram(aes(x=sitediffs_hot, fill="hot"), alpha=0.5) +
    scale_fill_viridis_d(option="C") +
    facet_grid(Site ~ Sex, scales="free_y")
# geom_line(data = data.frame(x=c(1,5), y=c(1,1)), aes(x=x, y=y), size=2) +
#xlim(c(-7,15))
