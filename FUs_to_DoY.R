# match DoY to forcing units

#needs tpars from ppc.R

library(tidyverse)

gclim <- read.csv("data/all_clim_PCIC.csv", header=TRUE, stringsAsFactors=FALSE) %>%
    filter(forcing_type=="scaled_ristos")

tpars <- read.csv("transformed_parameters.csv", header=TRUE, stringsAsFactors=FALSE)

clim <- gclim
clim$siteyear <- paste(clim$Site, clim$Year, sep='')

# choose cold, hot, and medium years.

climsort <- filter(clim, DoY == 180) %>%
    arrange(sum_forcing)

ggplot(clim, aes(x=DoY, y=sum_forcing, color=Site, group=siteyear)) +
    geom_line() +
    theme(legend.position = "bottom")

# KettleRiver 2011 is cold, Kalamalka 1998 is hot, and PRT 2008 is middle of the road.

# index clim and tpars by Site-Year
tparstoy <-tpars[1:3000,]
KR2011 <- filter(clim, siteyear == "KettleRiver2011")

tpars <- tpars %>%
    group_by(sum_forcing) %>%
    mutate(closest_sum_forcing = min(KR2011$sum_forcing[KR2011$sum_forcing > sum_forcing])) %>%
    mutate(DoY_KR2011 = KR2011$DoY[which(KR2011$sum_forcing == closest_sum_forcing)]) %>%
    select(-closest_sum_forcing)

climsplit <- clim %>%
    split(.$siteyear)



closestDoY <- data_frame
for (i in 1:length(unique(tparstoy$siteyear))) {
    climtemp <- filter(clim, siteyear==unique(tparstoy$siteyear)[i])
    tparstemp <- filter(tparstoy, siteyear==unique(tparstoy$siteyear)[i])
    for (j in 1:nrow(tparstemp)) {
        closest_sum_forcing <- which.min(climtemp$sum_forcing_real[climtemp$sum_forcing_real > tparstemp$sum_forcing[j]])
        climtemp$DoY[closest_sum_forcing]
    }


}

find_closest_forcing <- function(transformed_pars, climsplit) {
    # choose year site subset from model output
    temptpars <- transformed_pars
    # get correct year site subset from climsplit

    loc <- which(names(climsplit) == names(temptpars))
    tempclim <- climsplit[[loc]]
    temptpars <- tparstoy[[1]]

    DoYhold <- c()
    for (i in 1:nrow(temptpars)) {
        closest_sum_forcing <- min(tempclim$sum_forcing_real[tempclim$sum_forcing_real > temptpars$sum_forcing[i]])
        DoYhold[i] <- tempclim$DoY[which(tempclim$sum_forcing_real == closest_sum_forcing)]
    }

    temptpars$DoY_model <- DoYhold

    return(temptpars)
}
#

map(tparstoy, find_closest_forcing, climsplit=climsplit)
