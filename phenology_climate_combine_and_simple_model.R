# try logit model with subset of real data
library(dplyr)
library(lubridate)
library(tidyr)
library(viridis)
library(ggplot2)

# Data --------------------------------------------------------------------


#rawdat <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/from Rita Wagner/data_cleaned/PGTIS_pheno_1997_2012_cleaned.csv', stringsAsFactors = FALSE)
phendat <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/data_formatted_and_derived//derived_phenophase.csv', stringsAsFactors = FALSE)
climdat <- read.csv('~/Documents/research_phd/data/Climate/formatted/PrinceGeorgeSTP.csv', header = TRUE)

# Functions ----------------------------------

#function to calculate heat sum using a linear threshold model of heat sum accumulation. If mean temp for the day is below a given threshold, no heat is accumulated. If mean temp for the day is above a given threshold, the heat accumulated = mean temp - threshold temp. This requires a dataframe with MeanTempC columns (Mean temperature, numeric) and Year and DoY (Day of Year) columns. threshold temp is a number. Returns a df with year, day of year, daily heat and heatsum columns

#eventually will need to modify so site id in here somewhere

calculate_heat_sum <- function(climate_df, threshold_temp) {
    # calculate heat added on a given day
    # days where no heat is added bc threshold temp not reached
    no_heat <- climate_df %>%
        filter(MeanTempC < threshold_temp) %>%
        mutate(Heat = 0)
    # days where heat is added bc threshold temp is reached
    heat <- climate_df %>%
        filter(MeanTempC >= threshold_temp) %>%
        mutate(Heat = MeanTempC - threshold_temp)

    clim <- rbind(no_heat, heat) %>%
        arrange(Year, DoY) %>%
        group_by(Year) %>%
        mutate(Heatsum = cumsum(Heat)) %>% # add heatsum
        select(Year, DoY, Heat, Heatsum)

    return(clim)
}

# Calculate heatsum -----------------

# climate data
clim <- subset(climdat, Year %in% unique(phendat$Year))
colnames(clim)[2] <- "DoY" #rename DayofYear to DoY

#calculate amount of heat per day assume no heating below 5 degrees and linear heating starting at 5

clim <- calculate_heat_sum(clim, 5)

ggplot(clim, aes(x=DoY, y=Heatsum, color = Year)) +
    geom_point() #check nothing weird happened

# RUN FOR ONLY WAGNER DATA -------------------------------
# Combine climate and phenology data

# phenology data

phen <- subset(phendat, Source == "Rita Wagner") #only Wagner
phen_pgtis <- subset(phendat, Site = "PGTIS") #all PGTIS

# mdat <- subset(rawdat, Sex == "MALE" & Source == "Rita Wagner") #male
#
# fdat <- subset(rawdat, Sex == "FEMALE" & Source == "Rita Wagner") #female


phendf <- merge(phen, clim) %>%
    select(Sex, Year, DoY, Index, Clone, Phenophase_Derived, Heat, Heatsum, Orchard) %>%
    arrange(Year, Sex, Index, Clone, DoY) %>%
    filter(!Phenophase_Derived==0) # drop trees that didn't flower
phendf$Phenophase_Derived <- as.factor(phendf$Phenophase_Derived)

mdf <- subset(phendf, Sex == "MALE")

fdf <- subset(phendf, Sex == "FEMALE")

#write.csv(phendf, "~/Documents/research_phenolology/data/phenology_heatsum.csv")

#Test that no data dropped unintentionally
nrow(mdf) + nrow(fdf) == nrow(phendf)

# Visualize the data

ggplot(phendf, aes(x = Phenophase_Derived, y=Heatsum, color = Sex)) +
    geom_violin() +
    ggtitle("Distribution of Heatsums at each phenophase") +
    facet_wrap(~ Year)

ggplot(phendf, aes(x = Heatsum, y = Phenophase_Derived, color = Sex)) +
    geom_jitter(shape = 1, alpha = .5)

ggplot(phendf, aes(x = Heatsum, y = Phenophase_Derived, color = Sex)) +
    geom_jitter(shape = 1, alpha = .5) +
    facet_wrap(~ Year)

ggplot(phendf, aes(x = Heatsum, y = Phenophase_Derived, color = Sex)) +
    geom_jitter(shape = 1, alpha = .5) +
    facet_wrap(Orchard ~ .)

# Ordered Logit Model (rethinking)

# read in data for priors
heatsum_priors_dat <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/orchard_heatsums_WalshWebber2008.csv')
# calculate priors

simple_beta_prior = 0.5
pre_pollination_summary <- heatsum_priors_dat %>%
    filter(period == 'pre-pollination') %>%
    summarise(mean=mean(Tsum_day, na.rm=TRUE), sd=sd(Tsum_day, na.rm=TRUE))
pre_pollination_summary <- pre_pollination_summary*simple_beta_prior

# all_poll_summary <- heatsum_priors_dat %>%
#     group_by(year) %>%
#     summarise(Tsum_day = sum(Tsum_day)) %>%
#     summarise(mean=mean(Tsum_day, na.rm=TRUE), sd=sd(Tsum_day, na.rm=TRUE))
# all_poll_summary <- all_poll_summary*simple_beta_prior

#fix Phenophase


#fit pollen model
m1 <- map(
    alist(
        Phenophase_Derived ~ dordlogit(phi, c(a1, a2)),
        phi <- b*Heatsum,
        b ~ dbeta(2,5),
        c(a1,a2) ~ dnorm(pre_pollination_summary$mean, pre_pollination_summary$sd)
    ),
    data = mdf,
    start= list(a1=100*.2, a2=250*.2)
)

#view model
post_m1 <- extract.samples(m1)
precis(m1)

post_m1 <-post_m1 %>%
    mutate(k = b, h1 = a1/b, h2 = a2/b)

dens(post_m1, show.HPDI=.5)

# now predict data from the model (there are only 6 sets of parameters chosen - not sure where that gets picked)
shortheat <- seq(from=0, to = max(phendf$Heatsum), length.out=100)
s <- numeric(0)
for (i in 1:length(post_m1)) {
    si <- apply(data.frame(shortheat = shortheat), 1,
                function(y) rordlogit(100, phi = post_m1[i,3]*y, a = post_m1[i,1:2]))# simulate data at all input values (shortheat) at a given parameter set
    si <- cbind(i, si)
    s <- rbind(s, si)
}

s <- data.frame(s)
colnames(s) <- c("param_set", shortheat)
post_pred_m1 <- gather(s, key = heatsum, value = state, -param_set)
post_pred_m1$heatsum <- as.numeric(post_pred_m1$heatsum)

ggplot(post_pred_m1, aes(y = heatsum, x = as.factor(state), color = as.factor(param_set))) +
    geom_violin() +
    scale_color_viridis(discrete = TRUE) +
    ggtitle('heatsum distribution at each phenophase\n for 6 sets of parameters')
   # facet_grid(state ~ param_set)

ggplot(mdf, aes(x = Heatsum, y = Phenophase_Derived)) +
    geom_point() +
    ggtitle('heatsum vs. state, real data from Prince George')

ggplot(post_pred_m1, aes(factor(state), y = heatsum, fill = as.factor(state))) +
    geom_violin(trim=FALSE) +
    scale_fill_viridis(discrete=TRUE) +
    ggtitle("Data simulated from parameters from fitted model")

ggplot(mdf, aes(factor(Phenophase_Derived), y = Heatsum, fill = as.factor(Phenophase_Derived))) +
    geom_violin(trim=FALSE) +
    scale_fill_viridis(discrete=TRUE) +
    ggtitle("Real data from Prince George")

#plot real data and predictions on same graph

data_for_merging <- mdf %>% select(Heatsum, Phenophase_Derived) %>%
    mutate(source = "data") %>%
    dplyr::rename(heatsum = Heatsum, state = Phenophase_Derived)

data_for_merging$state <- as.numeric(data_for_merging$state)

post_pred_for_merging <- post_pred_m1 %>% mutate(source = "model") %>%
    select(-param_set)

combined_preds_data <- rbind(data_for_merging, post_pred_for_merging)

ggplot(combined_preds_data, aes(factor(state), y = heatsum, fill = source)) +
    geom_violin(trim=FALSE) +
    scale_fill_viridis(discrete=TRUE)


# RUN FOR ALL PGTIS DATA -----------------------
# Combine climate and phenology data

# phenology data

#phen <- subset(phendat, Source == "Rita Wagner") #only Wagner
phen_pgtis <- subset(phendat, Site == "PGTIS") #all PGTIS

# mdat <- subset(rawdat, Sex == "MALE" & Source == "Rita Wagner") #male
#
# fdat <- subset(rawdat, Sex == "FEMALE" & Source == "Rita Wagner") #female


phendf <- merge(phen_pgtis, clim) %>%
    select(Sex, Year, DoY, Index, Clone, Phenophase_Derived, Heat, Heatsum, Orchard) %>%
    arrange(Year, Sex, Index, Clone, DoY) %>%
    filter(!Phenophase_Derived==0) # drop trees that didn't flower
phendf$Phenophase_Derived <- as.factor(phendf$Phenophase_Derived)

mdf <- subset(phendf, Sex == "MALE")

fdf <- subset(phendf, Sex == "FEMALE")

#write.csv(phendf, "~/Documents/research_phenolology/data/phenology_heatsum.csv")

#Test that no data dropped unintentionally
nrow(mdf) + nrow(fdf) == nrow(phendf)

# Visualize the data ------------------

ggplot(phendf, aes(x = Phenophase_Derived, y=Heatsum, color = Sex)) +
    geom_violin() +
    ggtitle("Distribution of Heatsums at each phenophase") +
    facet_wrap(~ Year)

ggplot(phendf, aes(x = Heatsum, y = Phenophase_Derived, color = Sex)) +
    geom_jitter(shape = 1, alpha = .5)

ggplot(phendf, aes(x = Heatsum, y = Phenophase_Derived, color = Sex)) +
    geom_jitter(shape = 1, alpha = .5) +
    facet_wrap(~ Year)

ggplot(phendf, aes(x = Heatsum, y = Phenophase_Derived, color = Sex)) +
    geom_jitter(shape = 1, alpha = .5) +
    facet_wrap(Orchard ~ .)

# Ordered Logit Model (rethinking)-----------------------------------------------------

# read in data for priors
heatsum_priors_dat <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/orchard_heatsums_WalshWebber2008.csv') #heat sum data from p. 49 of 2008 Tree Improvement Program Project Report. Victoria, B.C. The data is contained in Tables 12 and 13 of Section 5.2.9, Enhancing Lodgepole Pine Seed Set submitted by Chris Walsh and prepared by Joe Webber..

# calculate priors

simple_beta_prior = 0.5
pre_pollination_summary <- heatsum_priors_dat %>%
    filter(period == 'pre-pollination') %>%
    summarise(mean=mean(Tsum_day, na.rm=TRUE), sd=sd(Tsum_day, na.rm=TRUE))
pre_pollination_summary <- pre_pollination_summary*simple_beta_prior

# all_poll_summary <- heatsum_priors_dat %>%
#     group_by(year) %>%
#     summarise(Tsum_day = sum(Tsum_day)) %>%
#     summarise(mean=mean(Tsum_day, na.rm=TRUE), sd=sd(Tsum_day, na.rm=TRUE))
# all_poll_summary <- all_poll_summary*simple_beta_prior

#fix Phenophase


#fit pollen model
m1 <- map(
    alist(
        Phenophase_Derived ~ dordlogit(phi, c(a1, a2)),
        phi <- b*Heatsum,
        b ~ dbeta(2,5),
        c(a1,a2) ~ dnorm(pre_pollination_summary$mean, pre_pollination_summary$sd)
    ),
    data = mdf,
    start= list(a1=100*.2, a2=250*.2)
)

#view model
post_m1 <- extract.samples(m1)
precis(m1)

post_m1 <-post_m1 %>%
    mutate(k = b, h1 = a1/b, h2 = a2/b)

dens(post_m1, show.HPDI=.5)

# now predict data from the model (there are only 6 sets of parameters chosen - not sure where that gets picked)
shortheat <- seq(from=0, to = max(phendf$Heatsum), length.out=100)
s <- numeric(0)
for (i in 1:length(post_m1)) {
    si <- apply(data.frame(shortheat = shortheat), 1,
                function(y) rordlogit(100, phi = post_m1[i,3]*y, a = post_m1[i,1:2]))# simulate data at all input values (shortheat) at a given parameter set
    si <- cbind(i, si)
    s <- rbind(s, si)
}

s <- data.frame(s)
colnames(s) <- c("param_set", shortheat)
post_pred_m1 <- gather(s, key = heatsum, value = state, -param_set)
post_pred_m1$heatsum <- as.numeric(post_pred_m1$heatsum)

ggplot(post_pred_m1, aes(y = heatsum, x = as.factor(state), color = as.factor(param_set))) +
    geom_violin() +
    scale_color_viridis(discrete = TRUE) +
    ggtitle('heatsum distribution at each phenophase\n for 6 sets of parameters')
# facet_grid(state ~ param_set)

ggplot(mdf, aes(x = Heatsum, y = Phenophase_Derived)) +
    geom_point() +
    ggtitle('heatsum vs. state, real data from Prince George')

ggplot(post_pred_m1, aes(factor(state), y = heatsum, fill = as.factor(state))) +
    geom_violin(trim=FALSE) +
    scale_fill_viridis(discrete=TRUE) +
    ggtitle("Data simulated from parameters from fitted model")

ggplot(mdf, aes(factor(Phenophase_Derived), y = Heatsum, fill = as.factor(Phenophase_Derived))) +
    geom_violin(trim=FALSE) +
    scale_fill_viridis(discrete=TRUE) +
    ggtitle("Real data from Prince George")

#plot real data and predictions on same graph

data_for_merging <- mdf %>% select(Heatsum, Phenophase_Derived) %>%
    mutate(source = "data") %>%
    dplyr::rename(heatsum = Heatsum, state = Phenophase_Derived)

data_for_merging$state <- as.numeric(data_for_merging$state)

post_pred_for_merging <- post_pred_m1 %>% mutate(source = "model") %>%
    select(-param_set)

combined_preds_data <- rbind(data_for_merging, post_pred_for_merging)

ggplot(combined_preds_data, aes(factor(state), y = heatsum, fill = source)) +
    geom_violin(trim=FALSE) +
    scale_fill_viridis(discrete=TRUE)

# Ordered Logit Model with levels-----------------------------------------------------
#This can't work because of limitations of the rethinking package
m2 <- map2stan(
    alist(
        #likelihood
        Phenophase_Simp ~ dordlogit(phi, cutpoints),
        #model
        phi <- b*Heatsum + a_clone[Clone],
        #priors
        b ~ dbeta(2,5),
        cutpoints ~ dnorm(pre_pollination_summary$mean, pre_pollination_summary$sd),
        a_clone[Clone] ~ dnorm(0, sigma_clone),
        #hyperpriors
        sigma_clone ~ dcauchy(0,1)
    ),
    data = df,
    start= list(a1=100*.2, a2=250*.2),
    warmup = 10, inter=500, chains = 1, cores = 1
)
