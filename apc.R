# Average Predictive comparisons
# references Gelman & Pardoe 2007 + dchudz predcomps

source('phenology_functions.R')

library(tidyverse)
library(gtools)

fmod <- readRDS("2019-10-28phenologyFEMALE.rds") %>%
    as.data.frame() %>%
  select(-contains("state_rep")) %>%
  sample_frac(0.05)

mmod <- readRDS("2019-10-28phenologyMALE.rds") %>%
    as.data.frame() %>%
  select(-contains("state_rep")) %>%
  sample_frac(0.05)


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


#############################

# WEIGHTS

fmod <- readRDS("2019-10-28phenologyFEMALE.rds") %>%
  as.data.frame() %>%
  select(-contains("state_rep")) %>%
  sample_frac(0.05)

mmod <- readRDS("2019-10-28phenologyMALE.rds") %>%
  as.data.frame() %>%
  select(-contains("state_rep")) %>%
  sample_frac(0.05)

# identify combinations of effects that actually occur

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")

# Site effects
udf <- unique_grouper(fdf, mdf) %>%
  filter(Sex=="FEMALE") %>%
  distinct() %>%
  select("SPU_Name", "Clone", "Year") # v

weights <- matrix(nrow=nrow(udf), ncol = nrow(udf))

for (i in 1:nrow(udf)) { # This calculation is slow af
  for (j in i:nrow(udf)) {
    weights[i,j] <- sum(udf[i,] %in% udf[j,])
  }
}

#make weights matrix symmetric
weights[is.na(weights)] <- 0 # NA bottom half should be 0s so that addition works later
tweights <- t(weights)
diag(tweights) <- 0 # replace diagonal so self doesn't get counted twice

# calculate weights
weights <- (weights + tweights)/ncol(udf)

# renormalize weights
wdf <- as.data.frame(weights)
wdf$totweights <- rowSums(weights)
wdf <- wdf/wdf$totweights # 



############################


library(predcomps)

priceCoef <- -.12
qualityCoef <- .1
qualityNoiseStdDev <- 5
nWines=50000
nRowsForPlottingSample <- 1000

numForTransitionStart <- 500
numForTransitionEnd <- 10000
onlyIncludeNearestN = 100

priceQualitySlope <- .4

df1 <- local({
  price <- sample(20:120, nWines, replace=TRUE)
  quality <- price * priceQualitySlope + 22 + rnorm(nWines, sd=qualityNoiseStdDev)
  purchaseProbability <- logistic(priceCoef*(price - 70) + qualityCoef*(quality - 50)  )
  purchased  <- rbinom(n = length(purchaseProbability), size=1, prob=purchaseProbability)
  data.frame(Quality = quality, 
             Price = price, 
             PurchaseProbability = purchaseProbability, 
             Purchased = purchased)
})


df1Sample <- df1[sample.int(nWines, size=nRowsForPlottingSample), ]
qplot(Price, Quality, alpha=I(.5), data = df1Sample) + 
  expand_limits(y=c(0,100))

logitFit1 <- glm(Purchased ~ Price + Quality, data = df1, family = "binomial")
logitFit1

predfun <- function(df) {
  with(df, Price + Quality)
}
df1$y <- predfun(df1)

apc1 <- GetPredCompsDF(predfun, df1, inputVars = c("Price", "Quality"),
                       numForTransitionStart = numForTransitionStart,
                       numForTransitionEnd = numForTransitionEnd,
                       onlyIncludeNearestN = onlyIncludeNearestN)

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