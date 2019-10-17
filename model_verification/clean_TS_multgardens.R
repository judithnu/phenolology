### 5 March 2019
## Sort of clean NPN data from botanic gardens

# Based off clean_TS.R by Cat
## Updated on 4 July 2019 based on new clean_TS.R code from C Chamberlain ##
## Updated on 25 July 2019 based on email from Cat where she told me what I missed from new clean_TS.R code ##

## How to download data:
# 1) Go to the NPN Data Downloader Tool: https://data.usanpn.org/observations
    # And go to the Phenology Observation Portal
# 2) Select 'Individual Phenometrics' and press NEXT 
# 3) Set Date range applicable to your question and press 'Set Date' and NEXT (I picked 1 Jan 2010 to 31 Dec 2018)
# 4) Select 'Partner Groups' tab on left: press the + next to 
    # 'Botanic Gardens and Arboretums and select: 
       # 'Arnold Arboretum - Tree Spotters' (can also incl. Lilac which is 'Native and Indicator Observation Program')
       # Denver Botanic Gardens
       # New York Botanical Garden Forest Phenology
       # Santa Barbara Botanic Garden
       # Santa Fe Botanical Garden
       # The Morton Arboretum
       # Tucson Botanical Garden 
    # Press 'Set Groups' and NEXT
# 5) Select 'Output fields' tab on left: and select 'ObservedBy Person ID' and 'Multiple Observers' and 'NumYs in Series' and 'Multiple FirstY'
# 6) You then wait a while for the data to download .... 

## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## Load Libraries
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
# library(lubridate)


# Set Working Directory
setwd("~/Documents/git/projects/treegarden/treespotters")

# 1. Get the data - and clean!
d <- read.csv("input/individual_phenometrics_data_multgardens.csv", header=TRUE)

### First let's do the obligatory cleaning checks with citizen scienece data
  d <- d[(d$Multiple_FirstY>=1 | d$Multiple_Observers>0),] ## This selects data where multiple people observed the same phenophase
  d <- d[(d$NumYs_in_Series>=3),] ## This selects data again where the same phenophase was seen 3 times in a row
  d <- d[(d$NumDays_Since_Prior_No>=0 & d$NumDays_Since_Prior_No<=14),] ## And this limits to data where a no is followed by a yes, so that it is a new observation/new phenophase but has been detected within a reasonable timeframe

# change from NPN output to more digestible column names
  bb<-d%>%
    rename(lat=Latitude)%>%
    rename(long=Longitude)%>%
    rename(elev=Elevation_in_Meters)%>%
    rename(state=State)%>%
    rename(year=First_Yes_Year)%>%
    rename(month=First_Yes_Month)%>%
    rename(day=First_Yes_Day)%>%
    rename(doy=First_Yes_DOY)%>%
    rename(numYs=Multiple_Observers)%>%
    rename(phase=Phenophase_Description)%>%
    rename(id=Individual_ID)%>%
    rename(genus=Genus)%>%
    rename(species=Species)
  bb.pheno<-dplyr::select(bb, genus, species, Common_Name, phase, lat, long, elev, state, year, doy, numYs, id)
bb.pheno$phase<-ifelse(bb.pheno$phase=="Breaking leaf buds", "budburst", bb.pheno$phase)
bb.pheno$phase<-ifelse(bb.pheno$phase=="Leaves", "leafout", bb.pheno$phase)
bb.pheno$phase<-ifelse(bb.pheno$phase=="Flowers or flower buds", "flowers", bb.pheno$phase)
bb.pheno$phase<-ifelse(bb.pheno$phase=="Falling leaves", "leaf drop", bb.pheno$phase)

### Now work on finding day of budburst, etc.,
## Edited by Lizzie to include site and toss out some species
bb.pheno.allsp <- filter(bb.pheno, numYs>0) # get only data with multiple observers
bb.pheno.allsp$latbi <- paste(bb.pheno.allsp$genus, bb.pheno.allsp$species)
sort(unique(bb.pheno.allsp$latbi))
# Now I will make a list of species I know off-hand that we should toss. ...
tossme <- c("Achillea millefolium", "Agraulis vanillae", "Alliaria petiolata" , "Aquilegia caerulea", "Artemisia tridentata",
            "Asclepias syriaca" , "Auriparus flaviceps", "Bombus spp.", "Bouteloua dactyloides", "Bouteloua gracilis",
            "Buteo jamaicensis", "Calypte anna", "Chamerion angustifolium", "Chilopsis linearis", "Cyanocitta cristata",
            "Danaus gilippus", "Danaus plexippus", "Ericameria nauseosa", "Erythronium americanum", "Eschscholzia californica",
            "Eurybia divaricata", "Fouquieria splendens", "Impatiens capensis", "Ipomopsis aggregata", "Larrea tridentata",
            "Monarda fistulosa", "Ranunculus ficaria", "Ratibida columnifera", "Scilla siberica", "Simmondsia chinensis",
            "Sphaeralcea coccinea", "Turdus migratorius", "Xylocopa varipuncta",  "Yucca elata")
bb.pheno <- bb.pheno.allsp[which(!bb.pheno.allsp$latbi %in% tossme),]
sort(unique(bb.pheno$latbi))

# Choose the rows with the earliest event date per location and ID ....
# Below, group each individual by phenophase and year to find the first observation (using the slice function), 
# so first day of budburst for that individual for that year
doy_pheno <- bb.pheno %>% 
  group_by(id, latbi, phase, year, lat, long, elev, state) %>% 
  slice(which.min(doy))
doy_pheno<-doy_pheno[!duplicated(doy_pheno),]

###################################
# plots to check for outliers ... #
###################################
# these are okay
doy_pheno_starttoripe <- subset(doy_pheno, phase=="budburst" | phase=="Ripe fruits")
ggplot(subset(doy_pheno_starttoripe, state=="MA"), aes(phase, doy, color=latbi)) + 
  geom_point() +
  facet_grid(latbi ~ year) +
  labs(x = "phases", y = "doy")
ggplot(subset(doy_pheno_starttoripe, state=="NY"), aes(phase, doy, color=latbi)) + 
  geom_point() +
  facet_grid(latbi ~ year) +
  labs(x = "phases", y = "doy")
# these are a little better
ggplot(subset(doy_pheno_starttoripe, state=="MA" & phase=="budburst"), aes(doy)) + 
  geom_histogram() +
  facet_grid(latbi ~ year)
ggplot(subset(doy_pheno_starttoripe, state=="MA" & phase=="Ripe fruits"), aes(doy)) + 
  geom_histogram() +
  facet_grid(latbi ~ year)
ggplot(subset(doy_pheno_starttoripe, state=="NY" & phase=="budburst"), aes(doy)) + 
  geom_histogram() +
  facet_grid(latbi ~ year)
ggplot(subset(doy_pheno_starttoripe, state=="NY" & phase=="Ripe fruits"), aes(doy)) + 
  geom_histogram() +
  facet_grid(year ~ latbi) # ripe fruit observations at start of year? This seems odd ...

checkme <- subset(doy_pheno, latbi=="Liriodendron tulipifera" & state=="NY" & phase=="Ripe fruits")
hist(checkme$doy)

# Sometimes, base R graphics are the best option
nydatbb <- subset(doy_pheno_starttoripe, state=="NY" & phase=="budburst")
madatbb <- subset(doy_pheno_starttoripe, state=="MA" & phase=="budburst")
nydatrf <- subset(doy_pheno_starttoripe, state=="NY" & phase=="Ripe fruits")
madatrf <- subset(doy_pheno_starttoripe, state=="MA" & phase=="Ripe fruits")
xlimnyrf <- c(min(nydatrf$doy, na.rm=TRUE), max(nydatrf$doy, na.rm=TRUE))

pdf(file.path("figures/cleaningplots_NY_bb.pdf"), width = 6, height = 4)
for (i in 1:length(sort(unique(nydatbb$latbi)))){
    subby <- subset(nydatbb, latbi==sort(unique(nydatbb$latbi))[i])
    hist(subby$doy, main=sort(unique(nydatbb$latbi))[i], xlab="doy of budburst")
    }
dev.off()
pdf(file.path("figures/cleaningplots_NY_rf.pdf"), width = 6, height = 4)
for (i in 1:length(sort(unique(nydatrf$latbi)))){
    subby <- subset(nydatrf, latbi==sort(unique(nydatrf$latbi))[i])
    hist(subby$doy, main=sort(unique(nydatrf$latbi))[i], xlab="doy of ripe fruit")
    }
dev.off()
pdf(file.path("figures/cleaningplots_MA_bb.pdf"), width = 6, height = 4)
for (i in 1:length(sort(unique(madatbb$latbi)))){
    subby <- subset(madatbb, latbi==sort(unique(madatbb$latbi))[i])
    hist(subby$doy, main=sort(unique(madatbb$latbi))[i], xlab="doy of budburst")
    }
dev.off()
pdf(file.path("figures/cleaningplots_MA_rf.pdf"), width = 6, height = 4)
for (i in 1:length(sort(unique(madatrf$latbi)))){
    subby <- subset(madatrf, latbi==sort(unique(madatrf$latbi))[i])
    hist(subby$doy, main=sort(unique(madatrf$latbi))[i], xlab="doy of ripe fruit")
    }
dev.off()

pdf(file.path("figures/cleaningplots_MA_bb_byyr.pdf"), width = 6, height = 4)
for (i in 1:length(sort(unique(madatbb$latbi)))){
    subby <- subset(madatbb, latbi==sort(unique(madatbb$latbi))[i])
    par(mfrow=c(2,3))
    for(yr in 1:length(unique(subby$year))){
        subbyyr <- subset(subby, year==unique(subby$year)[yr])
    hist(subbyyr$doy, main=paste(sort(unique(madatbb$latbi))[i], unique(subby$year)[yr]), xlab="doy of budburst")
    }
}
dev.off()
pdf(file.path("figures/cleaningplots_NY_bb_byyr.pdf"), width = 6, height = 4)
for (i in 1:length(sort(unique(nydatbb$latbi)))){
    subby <- subset(nydatbb, latbi==sort(unique(nydatbb$latbi))[i])
    par(mfrow=c(2,3))
    for(yr in 1:length(unique(subby$year))){
        subbyyr <- subset(subby, year==unique(subby$year)[yr])
    hist(subbyyr$doy, main=paste(sort(unique(nydatbb$latbi))[i], unique(subby$year)[yr]), xlab="doy of budburst")
    }
}
dev.off()
pdf(file.path("figures/cleaningplots_MA_rf_byyr.pdf"), width = 6, height = 4)
for (i in 1:length(sort(unique(madatrf$latbi)))){
    subby <- subset(madatrf, latbi==sort(unique(madatrf$latbi))[i])
    par(mfrow=c(2,3))
    for(yr in 1:length(unique(subby$year))){
        subbyyr <- subset(subby, year==unique(subby$year)[yr])
    hist(subbyyr$doy, main=paste(sort(unique(madatrf$latbi))[i], unique(subby$year)[yr]), xlab="doy of ripe fruit")
    }
}
dev.off()
pdf(file.path("figures/cleaningplots_NY_rf_byyr.pdf"), width = 6, height = 4)
for (i in 1:length(sort(unique(nydatrf$latbi)))){
    subby <- subset(nydatrf, latbi==sort(unique(nydatrf$latbi))[i])
    par(mfrow=c(2,3))
    for(yr in 1:length(unique(subby$year))){
        subbyyr <- subset(subby, year==unique(subby$year)[yr])
    hist(subbyyr$doy, main=paste(sort(unique(nydatrf$latbi))[i], unique(subby$year)[yr]), xlab="doy of ripe fruit")
    }
}
dev.off()
pdf(file.path("figures/cleaningplots_NY_rf_byyr_fixedxlim.pdf"), width = 6, height = 4)
for (i in 1:length(sort(unique(nydatrf$latbi)))){
    subby <- subset(nydatrf, latbi==sort(unique(nydatrf$latbi))[i])
    par(mfrow=c(2,3))
    for(yr in 1:length(unique(subby$year))){
        subbyyr <- subset(subby, year==unique(subby$year)[yr])
    hist(subbyyr$doy, main=paste(sort(unique(nydatrf$latbi))[i], unique(subby$year)[yr]), xlim=xlimnyrf,
        xlab="doy of ripe fruit")
    }
}
dev.off()


checkme <- subset(doy_pheno, latbi=="Nyssa sylvatica" & state=="NY" & phase=="Ripe fruits")
checkme

######################
## end of plotting ###
######################

# some random cleaning (found through plotting!) 
doy_pheno$doy<-ifelse(doy_pheno$species=="alleghaniensis" & doy_pheno$year==2016 & doy_pheno$doy==59,
    NA, doy_pheno$doy)


# Calculate mean dates
meandates <-
      ddply(doy_pheno, c("phase", "genus", "species", "latbi", "year", "lat", "long", "elev", "state"), summarise,
      mean = mean(doy))

meandatesplus <-
      ddply(doy_pheno, c("phase", "genus", "species", "latbi", "year", "lat", "long", "elev", "state"), summarise,
      mean = mean(doy),
      sd = sd(doy),
      sem = sd(doy)/sqrt(length(doy)))


longdf <- spread(meandates, phase, mean)
names(longdf)[names(longdf)=="Ripe fruits"] <- "ripefruits"

# quick look at ripe fruit
howmanyripe <- subset(longdf, is.na(ripefruits)==FALSE)
ny <- subset(howmanyripe, state=="NY")
table(ny$year)
ny2018 <- subset(ny, year=="2018")
unique(ny2018$latbi)

write.csv(longdf, file="output/quickdirty_botanicgardens.csv", row.names=FALSE)



