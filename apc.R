# Average Predictive comparisons

fmod <- readRDS("slopes_nc_scaled_ristos_FEMALE2019-10-04climatena.rds") %>%
    as.data.frame()

mmod <- readRDS("slopes_nc_scaled_ristos_MALE2019-10-04climatena.rds") %>%
    as.data.frame()

state_repf <- t(as.matrix(dplyr::select(fmod, contains("state_rep"))))
#state_repm <- t(as.matrix(dplyr::select(mmod, contains("state_rep"))))

reps <- c(1:1800)
#colnames(state_repm) <- paste("staterep", reps, sep="_")
colnames(state_repf) <- paste("staterep", reps, sep="_")

# original data (this code should match relevant bits in run_stan)
phendf <- read_data()

# identify combinations of effects that actually occur

fdf <- splitsex(phendf, "FEMALE")
#mdf <- splitsex(phendf, "MALE")


#udf <- unique_grouper(fdf, mdf)
