# Drop redundant data
# Observations at Prince George that aren't the last recorded "not started" day, the first or last recorded active day, or the first recorded "finished" day are redundant.
# Keep for each unique Tree, Clone, Year, Site, Prov combination at PGTIS
## Last recorded 1
## First recorded 2
## Last recorded 2
## First recorded 1

library(tidyverse)


data <- read.csv("~/Documents/research_phenolology/data/phenology_heatsum.csv", header=TRUE, stringsAsFactors = FALSE)
data$Phenophase_Derived <- as.factor(data$Phenophase_Derived)

# Drop redundant data from PGTIS
slim_pgtis <- filter(data, Site == "PGTIS") %>%
    group_by(Index, Phenophase_Derived) %>%
    mutate(FR = min(DoY)) %>%
    mutate(LR = max(DoY)) %>%
    mutate(keep = case_when(Phenophase_Derived == 1 & DoY==LR ~ 1, #last recorded not started
                            Phenophase_Derived == 2 & DoY==FR ~ 1, #first recorded flowering
                            Phenophase_Derived == 2 & DoY==LR ~ 1, #last recorded flowering
                            Phenophase_Derived == 3 & DoY==FR ~ 1)) %>% #first recorded finished
    filter(keep == 1) %>%
    select(-FR, -LR, -keep)

# Combine PGTIS and rest of data

no_pgtis <- filter(data, Site !="PGTIS")

data_slimmed <- full_join(slim_pgtis, no_pgtis)

write.csv(data_slimmed, "~/Documents/research_phenolology/data/phenology_heatsum.csv", row.names = FALSE)

# don't drop data you don't mean to
nrow(data)
nrow(data_slimmed)

nrow(no_pgtis)+nrow(slim_pgtis) == nrow(data_slimmed)
