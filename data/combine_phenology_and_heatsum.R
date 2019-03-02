# Combine phenology and climate data

# combine phenology and weather station data
library(dplyr)
library(lubridate)
library(viridis)
library(ggplot2)

# Data --------------------------------------------------------------------


#rawdat <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/from Rita Wagner/data_cleaned/PGTIS_pheno_1997_2012_cleaned.csv', stringsAsFactors = FALSE)
phendat <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/data_formatted_and_derived//derived_phenophase.csv', stringsAsFactors = FALSE)
climdat <- read.csv('~/Documents/research_phd/data/Climate/formatted/PrinceGeorgeSTP_VernonNorth.csv', header = TRUE)

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

    clim <- rbind(no_heat, heat)
    clim <- clim %>%
        arrange(Year, DoY) %>%
        group_by(Site, Year) %>%
        mutate(Heatsum = cumsum(Heat)) %>% # add heatsum
        select(Site, Year, DoY, Heat, Heatsum)

    return(clim)
}

# Calculate heatsum -----------------

# climate data
clim <- subset(climdat, Year %in% unique(phendat$Year))


#calculate amount of heat per day assume no heating below 5 degrees and linear heating starting at 5

clim <- calculate_heat_sum(clim, 5)

ggplot(clim, aes(x=DoY, y=Heatsum, color = Year)) +
    geom_point() +
    facet_wrap(Site ~ .) #check nothing weird happened

# RUN FOR ONLY WAGNER DATA -------------------------------
# Combine climate and phenology data

# phenology data

#phen <- subset(phendat, Source == "Rita Wagner") #only Wagner
phen <- subset(phendat, Site %in% c("Kalamalka","PGTIS")) #Kalamalka and PGTIS

# mdat <- subset(rawdat, Sex == "MALE" & Source == "Rita Wagner") #male
#
# fdat <- subset(rawdat, Sex == "FEMALE" & Source == "Rita Wagner") #female


phendf <- merge(phen, clim) %>%
    select(Site, Sex, Year, DoY, Index, Clone, Phenophase_Derived, Heat, Heatsum, Orchard) %>%
    arrange(Site, Year, Sex, Index, Clone, DoY) %>%
    filter(!Phenophase_Derived==0) # drop trees that didn't flower
phendf$Phenophase_Derived <- as.factor(phendf$Phenophase_Derived)

mdf <- subset(phendf, Sex == "MALE")

fdf <- subset(phendf, Sex == "FEMALE")

write.csv(phendf, "~/Documents/research_phenolology/data/stan_input/phenology_heatsum.csv", row.names = FALSE)

#Test that no data dropped unintentionally
nrow(mdf) + nrow(fdf) == nrow(phendf)

# Visualize the data

ggplot(phendf, aes(x = Phenophase_Derived, y=Heatsum, color = Sex, fill = Site)) +
    geom_boxplot() +
    ggtitle("Distribution of Heatsums at each phenophase") +
    facet_wrap(. ~ Year) +
    scale_fill_viridis_d() +
    scale_color_manual(values = c("black", "white"))

ggplot(phendf, aes(x = Heatsum, y = Phenophase_Derived, color = Sex)) +
    geom_jitter(shape = 1, alpha = .5) +
    facet_grid(Sex ~ Site)

ggplot(phendf, aes(x = Heatsum, y = Phenophase_Derived, color = Sex)) +
    geom_jitter(shape = 1, alpha = .5) +
    facet_wrap(Site ~ Year)

ggplot(phendf, aes(x = Heatsum, y = Phenophase_Derived, color = Sex)) +
    geom_jitter(shape = 1, alpha = .5) +
    facet_wrap(Orchard ~ .)
