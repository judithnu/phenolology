# Average Predictive comparisons

library(tidyverse)
library(gtools)

fmod <- readRDS("slopes_nc_scaled_ristos_FEMALE2019-10-04climatena.rds") %>%
    as.data.frame() %>%
  select(-contains("state_rep")) %>%
  sample_frac(0.2)

mmod <- readRDS("slopes_nc_scaled_ristos_MALE2019-10-04climatena.rds") %>%
    as.data.frame() %>%
  select(-contains("state_rep")) %>%
  sample_frac(0.2)


# original data (this code should match relevant bits in run_stan)
phendf <- read_data(slim = FALSE)

# identify combinations of effects that actually occur

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")


udf <- unique_grouper(fdf, mdf)

fpardf <- build_par_df(mcmcdf = fmod, datdf = fdf, sex = "FEMALE")
mpardf <- build_par_df(mcmcdf = mmod, datdf = mdf, sex = "MALE")

pardf <- rbind(fpardf, mpardf)


# calculate the start of forcing for all model configurations for the first "observation"
udf[1,]


site1 <- select(fmod, `b_site[1]`)
site2 <- select(fmod, `b_site[2]`)
site3 <- select(fmod, `b_site[3]`)
prov <- select(fmod, matches(paste("b_prov", fdf$ProvenanceID[1], sep="\\[")))
clone <- select(fmod, matches(paste("b_clone", fdf$CloneID[1], sep="\\[")))
year <- select(fmod, matches(paste("b_year", fdf$YearID[1], sep="\\[")))
beta <- select(fmod, beta)

betas1 <- beta + site1 + prov + clone + year
betas2 <- beta + site2 + prov + clone + year
betas3 <- beta + site3 + prov + clone + year
kappa1 <- select(fmod, `kappa[1]`)
kappa2 <- select(fmod, `kappa[2]`)

calcstageforcing <- function(p=0.2, beta=betas, kappa) {
  prob <- (logit(p) + kappa)/beta
  return(prob)
}

start1 <- calcstageforcing(p=0.2, beta=betas1, kappa=kappa1)
start2 <- calcstageforcing(p=0.2, beta=betas2, kappa=kappa1)
start3 <- calcstageforcing(p=0.2, beta=betas3, kappa=kappa1)
end <- calcstageforcing(p=0.8, kappa=kappa2)

diff1 <- (start1-start2)^2
diff2 <- (start1-start3)^2
diff3 <- (start2-start3)^2
diff <- cbind(diff1, diff2, diff3)
diffsum <- colSums(diff) # sum over model configs
usum <- sum(diffsum)# sum over param combinations
denom <- nrow(udf) * nrow(fmod) * ncol(diff)

siteapc <- usum/denom

##############################
x <- 0:50
mu1 <- .5*x 
mu2 <- 2*x 


site.1 <- rnorm(length(x), mu1, 0.5)
site.2 <- rnorm(length(x), mu2, 0.5)


dat <- data.frame(cbind(x, site.1, site.2) ) %>%
  pivot_longer(contains("site")) %>%
  separate(name, into = c("site", "siteID")) %>%
  select(-site)

library(rethinking)
tfit <- ulam(
  alist(
    value ~ normal(mu, sigma),
    mu <- beta[siteID]*x,
    beta[siteID] ~ normal(0,1),
    sigma ~ exponential(2)
    ),
  data=dat, declare_all_data = FALSE)

summary(tfit)

samples <- as.data.frame(extract.samples(tfit))

lin <- function(b,x) {
  y <- b*x
  return(y)
}

e1 <- sapply(dat$x, lin, b=samples$beta.1)
e2 <- sapply(dat$x, lin, b=samples$beta.2)

dat %>% group_by(siteID) %>%
  summarise(m = mean(value))

diff <- ((e2-e1)^2)
modelconfigsums <- rowSums(diff) #maybe should be rowsums!
datasums <- sum(modelconfigsums)
denom <- 500*102
sqrt(datasums/denom)



# On average, site 2 differs from site 1 by 43 birds/frog seen.

#OR (not RMS)

diff <- e2-e1
modelconfigsums <- rowSums(diff) #maybe should be rowsums!
datasums <- sum(modelconfigsums)
denom <- 500*102
datasums/denom

# For every frog you see, you'll see 37 more birds at site 2.