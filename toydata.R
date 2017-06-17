# Make simple dataset for building & testing models

## CLIMATE DATA
clim <- data.frame(
    date = seq.Date(as.Date("2017-01-01"), by = "day", length.out = 365),
    temp = 5
)

## EVENT DATA
## Data for 1 tree at one site 20 and 80% phenology data

phene <- data.frame(
    tree = 1,
    pct_flowering = c(0.2, 0.8),
    date = as.Date(c("2017-06-01", "2017-06-02"))
)

## SURVEY DATA
## Data for 1 tree at one site. below+ "below 20% flowering" between = "at least 20% flowering and below 80% flowering" and above = "at or above 80% flowering"

phens <- data.frame(
    tree = 1,
    date = seq.Date(as.Date("2017-06-01"), as.Date("2017-06-07"), length.out = 4),
    phase = c("below", "between", "below", "above")
)

