# Make simple dataset for building & testing models

## EVENT DATA
## Data for 1 tree at one site with constant temperatures and yes/no phenology data
n = 5 # How many days

dat <- data.frame(
    day = seq(1:n),
    temp = rep(5, n),
    tree = 1,
    pct_budburst = c(0,0,0,1,1),
    accumulated_temp = NA)

