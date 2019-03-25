# Combine phenology and climate data

# combine phenology and weather station data
library(dplyr)
library(lubridate)
library(viridis)
library(ggplot2)
library(tidyr)

# Data --------------------------------------------------------------------


#rawdat <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/from Rita Wagner/data_cleaned/PGTIS_pheno_1997_2012_cleaned.csv', stringsAsFactors = FALSE)
phendat <- read.csv('~/Documents/research_phd/data/PhenologyAndPollenCounts/data_formatted_and_derived/inferred_derived_phenology.csv', stringsAsFactors = FALSE)
phendat$TreeID <- group_indices(phendat, Site, Orchard, Clone, Tree, X, Y)

climdat <- read.csv('~/Documents/research_phd/data/Climate/formatted/PCIC_all_seed_orchard_sites.csv', header = TRUE) %>%
    mutate(DoY = yday(Date))

# Functions ----------------------------------

#function to calculate heat sum using a linear threshold model of heat sum accumulation. If mean temp for the day is below a given threshold, no heat is accumulated. If mean temp for the day is above a given threshold, the heat accumulated = mean temp - threshold temp. This requires a dataframe with MeanTempC columns (Mean temperature, numeric) and Year and DoY (Day of Year) columns. threshold temp is a number. Returns a df with year, day of year, daily heat and heatsum columns

#eventually will need to modify so site id in here somewhere

calculate_GDD <- function(climate_df, threshold_temp) {
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
        select(Site, Year, DoY, Heat, Heatsum)

    return(clim)
}

#calculate forcing units according to Sarvas 1972/4/Hanninen 1990/Chuine ForcSar 1999
calculate_forcing_units <- function(climate_df) { #climatedf is a dataframe with a column named "mean_temp" containing mean temperatures in Celcius. columns named Site, Year, and DoY must also exist
    #calculate the amount of forcing any mean daily temperature provices
    climate_df$ristos <- 28.4/(1+exp(-0.185*(climate_df$mean_temp-18.4)))
    #sum the forcing units
    clim <- climate_df
    clim <- clim %>%
        dplyr::arrange(Site, Year, DoY) %>%
        dplyr::group_by(Site, Year) %>%
        dplyr::mutate(forcing_units = cumsum(ristos)) %>%
        select(Site, Year, DoY, mean_temp, ristos, forcing_units)
    return(clim)
}

# Calculate heatsum -----------------


# Climate data processing
clim <- mutate(climdat, Year = year(Date))
clim <- subset(clim, Year %in% unique(phendat$Year)) # drop climate data not in phenology dataset




# Calculate Forcing Units and Forcing sums

risto <- calculate_forcing_units(climate_df = clim)
clim5 <- calculate_GDD(clim, threshold_temp = 5)

clim <- dplyr::full_join(risto, clim5)

#check nothing weird happened
ggplot(clim, aes(x=DoY, y=Heatsum, color = Year)) +
    geom_point() +
    facet_wrap("Site") +
    theme_bw(base_size=18)

ggplot(clim, aes(x=DoY, y=forcing_units, color=Year)) +
    geom_point() +
    facet_wrap("Site") +
    theme_bw(base_size=18)

# Visualizations that should be moved somewhere else ##################
activeperiod <- filter(clim, DoY < 181) #Jan thru June only
ggplot(activeperiod, aes(x=Heatsum, y=forcing_units, color=Year)) +
    geom_point() +
    facet_wrap("Site") +
    theme_bw(base_size=18)

# Forcing units begin to accumulate before heatsum does

activeperiod <- filter(clim, DoY < 181) #Jan thru June only
ggplot(activeperiod, aes(x=Heat, y=ristos)) +
    geom_point() +
    #facet_wrap("Site") +
    theme_bw(base_size=18) +
    geom_abline()


ggplot(activeperiod_g, aes(x=mean_temp, y=daily_forcing, color=daily_forcing_type)) +
    geom_line(lwd=1.2) +
    facet_wrap("Site") +
    theme_bw(base_size=14) +
    geom_abline()


# Tidy forcing unit df ##########
gclim <- clim %>%
    gather(key=daily_forcing_type, value=daily_forcing, ristos, Heat) %>%
    tidyr::gather(key=cummulation_type, value=forcing_accum, forcing_units, Heatsum) %>%
    dplyr::filter(!(daily_forcing_type == "ristos" & cummulation_type == "Heatsum")) %>%
    dplyr::filter(!(daily_forcing_type == "Heat" & cummulation_type == "forcing_units"))


# More stuff that should be moved out of this file ##########
earliest_accummulation <- gclim %>% group_by(Site, Year, daily_forcing_type) %>%
    filter(forcing_accum > 5) %>%
    dplyr::summarize(earliest_accummulation = min(DoY))

gclim <- full_join(gclim, earliest_accummulation)

ggplot(gclim, aes(x=DoY, y=forcing_accum, color=daily_forcing_type, line_type=as.factor(Year))) +
    geom_line(alpha=0.7) +
    facet_wrap("Site") +
    geom_vline(aes(xintercept=earliest_accummulation, color=daily_forcing_type, alpha=0.5)) +
    theme_bw(base_size = 14) +
    scale_color_viridis_d(end=.85, option="A") +
    xlim(0,180) +
    ylim(0,1000)

# Plot climate data ####################
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


phendf <- merge(phen, gclim) %>%
    select(-Tree, -X, -Y, -Source, -First_RF, -Last_RF, -cummulation_type) %>%
    #select(Site, Sex, Year, DoY, Index, Clone, TreeID, Phenophase_Derived, Heat, Heatsum, Orchard) %>%
    arrange(daily_forcing_type, Site, Year, Sex, Index, Clone, DoY) %>%
    filter(!Phenophase_Derived==0) # drop trees that didn't flower

phendf$Phenophase_Derived <- as.factor(phendf$Phenophase_Derived)

write.csv(phendf, "~/Documents/research_phenolology/data/stan_input/phenology_heatsum.csv", row.names = FALSE)

## Visualizations, checks, and threshold test

mdf <- subset(phendf, Sex == "MALE")

fdf <- subset(phendf, Sex == "FEMALE")


#Test that no data dropped unintentionally
nrow(mdf) + nrow(fdf) == nrow(phendf)

# Visualize the data # MOVE
ggplot(mdf, aes(x = forcing_accum, color=Site)) +
    stat_ecdf() +
    facet_grid(daily_forcing_type ~ Phenophase_Derived)

ggplot(mdf, aes(x=forcing_accum, fill=Site)) +
    geom_density() +
    facet_grid(daily_forcing_type ~ Phenophase_Derived) +
    ggtitle("Male")

ggplot(fdf, aes(x=forcing_accum, fill=Site)) +
    geom_density() +
    facet_grid(daily_forcing_type ~ Phenophase_Derived) +
    ggtitle("Female")

ggplot(mdf, aes(x=forcing_units, color=Site)) +
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
