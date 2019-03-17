# Combine phenology and climate data

# combine phenology and weather station data
library(dplyr)
library(lubridate)
library(viridis)
library(ggplot2)
library(tidyr)

# Data --------------------------------------------------------------------


#rawdat <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/from Rita Wagner/data_cleaned/PGTIS_pheno_1997_2012_cleaned.csv', stringsAsFactors = FALSE)
phendat <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/data_formatted_and_derived/derived_phenophase.csv', stringsAsFactors = FALSE)
phendat$TreeID <- group_indices(phendat, Site, Orchard, Clone, Tree, X, Y)
phendat <- filter(phendat, Site !="Quesnel") #drop quesnel until climate data obtained

climdat <- read.csv('~/Documents/research_phd/data/Climate/formatted/PCIC_all_seed_orchard_sites.csv', header = TRUE) %>%
    mutate(DoY = yday(Date))

# Functions ----------------------------------

#function to calculate heat sum using a linear threshold model of heat sum accumulation. If mean temp for the day is below a given threshold, no heat is accumulated. If mean temp for the day is above a given threshold, the heat accumulated = mean temp - threshold temp. This requires a dataframe with MeanTempC columns (Mean temperature, numeric) and Year and DoY (Day of Year) columns. threshold temp is a number. Returns a df with year, day of year, daily heat and heatsum columns

#eventually will need to modify so site id in here somewhere

calculate_heat_sum <- function(climate_df, threshold_temp) {
    # calculate heat added on a given day
    # days where no heat is added bc threshold temp not reached
    no_heat <- climate_df %>%
        filter(mean_temp < threshold_temp) %>%
        dplyr::mutate(Heat = 0)
    # days where heat is added bc threshold temp is reached
    heat <- climate_df %>%
        filter(mean_temp >= threshold_temp) %>%
        dplyr::mutate(Heat = mean_temp - threshold_temp)

    clim <- rbind(no_heat, heat)
    clim <- clim %>%
        dplyr::arrange(Site, Year, DoY) %>%
        dplyr::group_by(Site, Year) %>%
        dplyr::mutate(Heatsum = cumsum(Heat)) %>% # add heatsum
        select(Site, Year, DoY, Heatsum)

    return(clim)
}

# Calculate heatsum -----------------

# climate data
clim <- mutate(climdat, Year = year(Date))
clim <- subset(clim, Year %in% unique(phendat$Year))





#calculate amount of heat per day assume no heating below 5 degrees and linear heating starting at 5

clim0 <- calculate_heat_sum(clim, 0)
clim1 <- calculate_heat_sum(clim, 1)
clim2 <- calculate_heat_sum(clim, 2)
clim3 <- calculate_heat_sum(clim, 3)
clim4 <- calculate_heat_sum(clim, 4)
clim5 <- calculate_heat_sum(clim, 5)
clim6 <- calculate_heat_sum(clim, 6)

clim <- clim5

#check nothing weird happened
ggplot(heatsum_vary_thresholds, aes(x=DoY, y=Heatsum, color = Year)) +
    geom_point() +
    facet_grid(Site ~ threshold) +
    theme_bw(base_size=18)

climdat$Year <- lubridate::year(climdat$Date)
climdat$Month <- lubridate::month(climdat$Date)
ggplot(climdat, aes(x=as.factor(Month), y=mean_temp, fill=Site)) +
    geom_violin(draw_quantiles = c(0.5)) +
    #theme_bw() +
    xlab("Month") +
    #theme(axis.text.x=element_blank()) +
    theme_bw(base_size=20) +
    #facet_wrap("Site") +
    theme(legend.position = "top") +
    scale_fill_viridis_d()+
    ggtitle("Mean Temperature 1997-2011")


# Combine climate and phenology data ---------------

# phenology data


phen <- phendat

# mdat <- subset(rawdat, Sex == "MALE" & Source == "Rita Wagner") #male
#
# fdat <- subset(rawdat, Sex == "FEMALE" & Source == "Rita Wagner") #female


phendf <- merge(phen, clim5) %>%
    rename(Heatsum=Heatsum5) %>%
    select(-Tree, -X, -Y, -Source, -First_RF, -Last_RF) %>%
    #select(Site, Sex, Year, DoY, Index, Clone, TreeID, Phenophase_Derived, Heat, Heatsum, Orchard) %>%
    arrange(Heatsum, Site, Year, Sex, Index, Clone, DoY) %>%
    filter(!Phenophase_Derived==0) # drop trees that didn't flower
phendf$Phenophase_Derived <- as.factor(phendf$Phenophase_Derived)

write.csv(phendf, "~/Documents/research_phenolology/data/stan_input/phenology_heatsum.csv", row.names = FALSE)

## Visualizations, checks, and threshold test

mdf <- subset(phendf, Sex == "MALE")

fdf <- subset(phendf, Sex == "FEMALE")

# Determine best threshold
ggplot(mdf, aes(x=Heatsum, color=Site))+
    stat_ecdf() +
    facet_grid(threshold ~ Phenophase_Derived) +
    ggtitle("Male phenophase transitions with varying heating thresholds") +
    theme_bw()

ggplot(fdf, aes(x=Heatsum, color=Site))+
    stat_ecdf() +
    facet_grid(threshold ~ Phenophase_Derived) +
    ggtitle("Female phenophase transitions with varying heating thresholds") +
    theme_bw()



#Test that no data dropped unintentionally
nrow(mdf) + nrow(fdf) == nrow(phendf)

# Visualize the data
ggplot(mdf, aes(x = Heatsum, color=Site)) +
    stat_ecdf() +
    facet_wrap("Phenophase_Derived")

phendf <- phendf %>%
    filter(Phenophase_Derived %in% c(2,3))

pre2006 <- phendf %>% #MAYBE PRESENTATION?
    filter(Year < 2006)
ggplot(phendf, aes(x=Heatsum, linetype=Sex, color=Site)) +
    stat_ecdf() +
    facet_grid(Year ~ Phenophase_Derived) +
    theme_bw() +
    theme(legend.position="bottom") +
    scale_color_viridis_d()

ggplot(pre2006, aes(x=Heatsum, linetype=Sex, color=Site)) +
    stat_ecdf() +
    facet_grid(Year ~ Phenophase_Derived) +
    theme_bw() +
    theme(legend.position="bottom") +
    scale_color_viridis_d() +
    ggtitle("PGTIS only")

ggplot(mdf, aes(x = Phenophase_Derived, y=Heatsum, color = Sex, fill = Site)) +
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
