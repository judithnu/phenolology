# try logit model with subset of real data
library(dplyr)

rawdat <- read.csv('/home/sus/Documents/research_phd/data/PhenologyAndPollenCounts/from Rita Wagner/data_cleaned/PGTIS_pheno_1997_2012_cleaned.csv', stringsAsFactors = FALSE)
rclim <- read.csv('/home/sus/Documents/research_phd/data/Climate/formatted/PrinceGeorgeSTP.csv', header = TRUE)

# phenology data: one year, male

pdat <- subset(rawdat, Year == 1997 & Sex == "MALE")
pdat$DayofYear <- lubridate::yday(pdat$Date)

rdat <- subset(rdat, Year == 1997 & Sex == "FEMALE")
rdat$DayofYear <- lubridate::yday(rdat$Date)

# climate data
clim <- subset(rclim, Year == 1997)

#calculate amount of heat per day assume no heating below 5 degrees and linear heating starting at 5
no_heat <- clim %>%
    filter(MeanTempC < 5) %>%
    mutate(Heat = 0)

heat <- clim %>%
    filter(MeanTempC >= 5) %>%
    mutate(Heat = MeanTempC - 5)

clim <- rbind(no_heat, heat) %>%
    arrange(DayofYear) %>%
    mutate(Heatsum = cumsum(Heat)) %>% # add heatsum
    select(DayofYear, Heat, Heatsum)

# combine climate and phenology data

df <- merge(pdat, clim) %>%
    select(DayofYear, Clone, Tree, Phenophase, Heat, Heatsum) %>%
    arrange(Tree, DayofYear) %>%
    filter(!Phenophase==0)# drop unexplained 0s

rdf <- merge(rdat, clim) %>%
    select(DayofYear, Clone, Tree, Phenophase, Heat, Heatsum) %>%
    arrange(Tree, DayofYear)

# transform phenology data into 1 (not started), 2 (active), (finished)


by_tree <- group_by(df, Tree)
rby_tree <- group_by(rdf, Tree)

fo <- by_tree %>%
    filter(Phenophase == '4') %>%
    summarise(First_Occurence = min(DayofYear))

rfo <- rby_tree %>%
    filter(Phenophase == '4') %>%
    summarise(First_Occurence = min(DayofYear))

intermed <- merge(df, fo) %>%# first occurance recorded
    arrange(Tree, DayofYear)
rintermed <- merge(rdf, rfo) %>%# first occurance recorded
    arrange(Tree, DayofYear)



intermed_ind <- which(intermed$DayofYear < intermed$First_Occurence) # not started
intermed$Phenophase_Simp <- NA
intermed$Phenophase_Simp[intermed_ind] <- 0
intermed$Phenophase_Simp[is.na(intermed$Phenophase_Simp) == TRUE & intermed$Phenophase == '4'] <- 1
intermed$Phenophase_Simp[is.na(intermed$Phenophase_Simp) == TRUE] <- 2
intermed

rintermed_ind <- which(rintermed$DayofYear < rintermed$First_Occurence) # not started
rintermed$Phenophase_Simp <- NA
rintermed$Phenophase_Simp[rintermed_ind] <- 0
rintermed$Phenophase_Simp[is.na(rintermed$Phenophase_Simp) == TRUE & rintermed$Phenophase == '4'] <- 1
rintermed$Phenophase_Simp[is.na(rintermed$Phenophase_Simp) == TRUE] <- 2
rintermed

df <- intermed %>%
    filter(Phenophase_Simp < 2)
rdf <- rintermed %>%
    filter(Phenophase_Simp < 2)

#males
flist <- alist(
    Phenophase_Simp ~ dbinom(1,  prob = p),
    logit(p) <- k * (Heatsum - (h + h_ind[Tree])),
    h_ind[Tree] ~ dnorm(0, sigma_ind),
    k ~ dnorm(mean = 0, sd = 0.25),
    h ~ dnorm(mean = 150, sd = 25),
    sigma_ind ~ dnorm(0,10)
)

m_bin <- map2stan(flist,
                  data = df,
                  iter = 1000,
                  chains = 1,
                  start = list(h = 150, k = -0.1, sigma_ind = 10)
)


post <- extract.samples(m_bin)
total_h_ind <- sapply(1:100, function(ind) post$h + post$h_ind[ , ind])
round(apply(total_h_ind, 2, mean), 2)

dens(post$k)
dens(post$h)
dens(unlist(total_h_ind), show.HPDI = .80)

for (i in 1:length(post)) {
    dens(post[[i]])
}

ggplot(df, aes(x = Heatsum, y = as.factor(Phenophase_Simp))) +
    geom_count() +
    xlab("Heatsum (Celsius)") +
    ylab("Shedding Pollen?") +
    ggtitle("Pollen shed vs heatsum \n 1997 at Prince George")

#females

ggplot(rdf, aes(x = Heatsum, y = as.factor(Phenophase_Simp))) +
    geom_count() +
    xlab("Heatsum (Celsius)") +
    ylab("Receptive?") +
    ggtitle("Receptivity vs heatsum \n 1997 at Prince George")

flist <- alist(
    Phenophase_Simp ~ dbinom(1,  prob = p),
    logit(p) <- k * (Heatsum - (h + h_ind[Tree])),
    h_ind[Tree] ~ dnorm(0, sigma_ind),
    k ~ dnorm(mean = 0.5, sd = 0.25),
    h ~ dnorm(mean = 150, sd = 25),
    sigma_ind ~ dnorm(0,10)
)

m_bin <- map2stan(flist,
                  data = rdf,
                  iter = 100,
                  chains = 1,
                  start = list(h = 150, k = .4, sigma_ind = 10)
)

