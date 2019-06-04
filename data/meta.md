# data

Data used for phenology modeling

## stan_input

Data used by stan directly

## all_clim_PCIC.csv

Forcing units at all sites for 1997-2011

## combine_phenology_and_heatsum_PCIC.R

Combine phenology data and forcing unit data calculated from PCIC data.

## combine_phenology_and_heatsum_weather_station.R

Combine phenology data with forcing unit data calculated from weather station data.

## phenology_heatsum.csv

Combined phenology and forcing unit data. Used by `run_stan.R` to calculate `.rdump` files in `stan_input` folder.
