# Compare ristos to heatsum

# Owens reports that heatsums of 500 are enough to get phenology going

library(tidyr)
library(ggplot2)
library(dplyr)

forcingraw <- read.csv("data/all_clim_PCIC.csv")
phenology <- read.csv("data/phenology_heatsum_all.csv") %>%
    filter(Phenophase_Derived==2)



minmaxdoy <- c(min(phenology$DoY), max(phenology$DoY))

forcing_gdd_ristos <- filter(forcingraw, forcing_type %in% c("gdd", "ristos"))

forcingwide <- pivot_wider(forcingraw, names_from = forcing_type, values_from = c(forcing, sum_forcing))
forcingwide$group <- group_indices(forcingwide, Site, Year)

forcingraw$siteyear <- group_indices(forcingraw, Site, Year, forcing_type)
ggplot(forcingwide, aes(x=sum_forcing_gdd, y=sum_forcing_scaled_ristos, group=Year)) +
    geom_line() +
    facet_wrap("Site") +
    xlim(c(0,600)) +
    ylim(c(0, 30)) +
    geom_vline(xintercept = 500)

ggplot(forcing_gdd_ristos, aes(x=DoY, y=sum_forcing, color=forcing_type, group=siteyear)) +
    geom_line(alpha=0.5) +
    geom_hline(yintercept = 500) +
    scale_color_viridis_d(option="cividis") +
    ggtitle("Ristos and GDD5 sums throughout the year")

ggplot(filter(forcing_gdd_ristos, DoY < 200), aes(x=DoY, y=sum_forcing, color=forcing_type, group=siteyear)) +
    geom_line(alpha=0.5) +
    geom_hline(yintercept = 500) +
    scale_color_viridis_d(option="cividis") +
    ggtitle("Ristos and GDD5 sums before about July 18") +
    geom_vline(xintercept = minmaxdoy)

ggplot(ggplot(filter(forcing_gdd_ristos, DoY < 200), aes(x=DoY, y=sum_forcing, color=forcing_type, group=siteyear)))


