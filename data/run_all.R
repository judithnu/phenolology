# Script to run all data processing scripts

# pick a computer
machine <- "72"
#machine <- "mycon"

# Combine phenology and heatsum data

## Reads in phenology data from
### ~/Documents/research_phd/data/PhenologyAndPollenCounts/data_formatted_and_derived/inferred_derived_phenology.csv' and
## weather data from
### '~/Documents/research_phd/data/Climate/formatted/PCIC_all_seed_orchard_sites_corrected.csv'

source('data/combine_phenology_and_heatsum_PCIC.R')

# Writes out a combined heatsum and phenology file to /data/phenology_heatsum_all.csv" and a climate only file to
# ~/Documents/research_phenolology/data/all_clim_PCIC.csv", row.names=FALSE)

# Remove redundant data

# Reads in ~/Documents/research_phenolology/data/phenology_heatsum_all.csv"

source('data/slim_data.R')

# Writes out ~/Documents/research_phenolology/data/phenology_heatsum.csv"