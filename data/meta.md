# data

Data used for phenology modeling

Some scripts in this file do data processing. They should be run in the following order
-1) `combine_phenology_and_heatsum_PCIC.R`
-2) `slim_data.R`

## stan_input

Data used by stan directly

## all_clim_PCIC.csv

Forcing units at all sites for 1997-2011

## combine_phenology_and_heatsum_PCIC.R

Calculate forcing unit data and combine with phenology data. Reads in `research_phd/data/PhenologyAndPollenCounts/data_formatted_and_derived/inferred_derived_phenology.csv` and `PCIC_all_seed_orchard_sites_adjusted.csv` Writes combined phenology and heatsum to `phenology_heatsum_all.csv` and only the weather and forcing unit data to `all_clim_PCIC.csv`

## combine_phenology_and_heatsum_weather_station.R

Combine phenology data with forcing unit data calculated from weather station data.

## phenology_heatsum.csv

Combined phenology and forcing unit data. Used by `run_stan.R` to calculate `.rdump` files in `stan_input` folder.

## slim_data.R

Prince George data is highly redundant due to method of recording data for all trees on all observation days. This script drops redundant data - keeping only records for the last recorded "not started" day, the first recorded "flowering" day, the last recorded "flowering" day, and the first recorded "done flowering" day. Reads in from `phenology_heatsum_all.csv` and outputs to `phenology_heatsum.csv`
