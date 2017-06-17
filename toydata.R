# Make simple dataset for building & testing models

## EVENT DATA
## Data for 1 tree at one site with constant temperatures and yes/no phenology data

phen <- data.frame(
    tree = 1,
    pct_flowering = c(0.2, 0.8),
    date = as.Date(c("2017-06-01", "2017-06-02"))
)

clim <- data.frame(
    date = seq.Date(as.Date("2017-01-01"), by = "day", length.out = 365),
    temp = 5
)

