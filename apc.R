# Average Predictive comparisons
# references Gelman & Pardoe 2007 + dchudz predcomps

source('phenology_functions.R')

library(tidyverse)
library(gtools)
library(assertthat)

########## FUNCTIONS ####################

calcstageforcing <- function(p=0.2, beta=betas, kappa) {
  prob <- (logit(p) + kappa)/beta
  return(prob)
}

# this function is based on a similar function in predcomps ###########
GetPairs <- function(X, u, v,
                     numForTransitionStart = NULL,
                     numForTransitionEnd = NULL,
                     onlyIncludeNearestN = NULL) {
  
  assertthat::assert_that(length(u) == 1) # make sure we have exactly 1 input var of interest
  for (columnName in c(u,v)) {
    assert_that(columnName %in% names(X))
    columnClass <- class(X[[columnName]])
    # if (!(columnClass) %in% c("integer", "numeric")) {
    #   stop(sprintf("Sorry, column %s is of class %s. I can only deal with integer and numeric types for now.", columnName, columnClass))
    # }
  }
  
  if (!is.null(numForTransitionStart)) {
    X1 <- X[sample.int(nrow(X), size=numForTransitionStart), c(v,u)] 
  } else {
    X1 <- X[c(v,u)] # all v with u
  }
  
  if (!is.null(numForTransitionEnd)) {
    X2 <- X[sample.int(nrow(X), size=numForTransitionEnd), c(v,u)]
  } else {
    X2 <- X[c(v,u)] # all v with u
  }
  
  X1$OriginalRowNumber <- 1:nrow(X1)
  X2$OriginalRowNumber.B <- 1:nrow(X2)
  
  vMatrix1 <- as.matrix(X1[,v])
  vMatrix2 <- as.matrix(X2[,v])
  
  
  # covV=cov(vMatrix2)
  # 
  #distMatrix <- apply(vMatrix1, 1, function(row) mahalanobis(vMatrix2, row, covV))
  
  # calculate weights for unordered categorical variables
  weights <- matrix(nrow=nrow(X), ncol = nrow(X)) #empty matrix for weights with a row and column 
  # for every row in udf (for every v)
  
  # calculate a matrix where the upper half holds every row-col comparison
  for (i in 1:nrow(vMatrix1)) { # This calculation is slow af
    for (j in i:nrow(vMatrix1)) {
      weights[i,j] <- sum(vMatrix1[i,] %in% vMatrix1[j,])
    }
  }
  
  # now make the matrix look like a times table where every row-row comparison is available 
  # starting from the row or the column
  #make weights matrix symmetric. fill out the weights matrix by flipping the upper triangle matrix across
  # the diagona
  weights[is.na(weights)] <- 0 # NA bottom half should be 0s so that addition works later
  tweights <- t(weights)
  diag(tweights) <- 0 # replace diagonal so self doesn't get counted twice
  # calculate weights for bottom triangle
  weights <- (weights + tweights)/length(v)
  
  distMatrix <- weights
  dim(distMatrix)
  
  colnames(distMatrix) <- 1:ncol(distMatrix)
  rownames(distMatrix) <- 1:nrow(distMatrix)
  distDF <- as.data.frame(as.table(distMatrix))
  names(distDF) <- c("OriginalRowNumber.B", "OriginalRowNumber", "Weight")
  
  
  if (!is.null(onlyIncludeNearestN)) {
    distDF <- distDF %>%
      group_by(OriginalRowNumber) %>%
      filter(rank(Weight, ties.method="random") < onlyIncludeNearestN)
  }
  
  
  pairs <- merge(X1, distDF, by = "OriginalRowNumber")
  pairs <- merge(X2, pairs, by = "OriginalRowNumber.B", suffixes = c(".B", ""))
  #pairs$Weight <- 1/(mahalanobisConstantTerm + pairs$MahalanobisDistance)
  
  # If we haven't sampled, then OriginalRowNumber == OriginalRowNumber.B means that 
  # the transition start and end are the same, so we should remove those rows.
  if (is.null(numForTransitionStart) && is.null(numForTransitionEnd)) {
    pairs <- subset(pairs, OriginalRowNumber != OriginalRowNumber.B)
  }
  
  # Renormalize weights:
  pairs <- pairs %>% group_by(OriginalRowNumber) %>% mutate(Weight = Weight/sum(Weight))
  pairs$Weight[is.nan(pairs$Weight)] <- 0
  
  return(data.frame(pairs))
}

#################

start_time <- Sys.time()
fmod <- readRDS("2019-10-28phenologyFEMALE.rds") %>%
  as.data.frame() %>%
  select(-contains("state_rep")) %>%
  #sample_frac(0.05)
  sample_n(10) 
fmod$iter <- 1:nrow(fmod)

mmod <- readRDS("2019-10-28phenologyMALE.rds") %>%
  as.data.frame() %>%
  select(-contains("state_rep")) %>%
  # sample_frac(0.05)
  sample_n(10)
mmod$iter <- 1:nrow(mmod)


# original data 
phendf <- read_data(slim = FALSE) 

# identify combinations of effects that actually occur

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")


udf <- unique_grouper(fdf, mdf) %>%
  filter(Sex=="FEMALE") %>%
  sample_n(10)

# pairs of rows and their weights
pairs <- GetPairs(X=udf, u="SiteID", v=c("ProvenanceID", "CloneID", "YearID"))

nrow(pairs) == (nrow(udf) * nrow(udf)) - nrow(udf)

# I need to rethink these functions to be smaller because (now that I'm not doing ppc) I don't need to keep every observation
# just udf

# pardf is a dataframe with nrow(fmod) x nrow(udf) rows that contains identifying information for data 
# and the parameter values associated with the data
pardf <- build_par_df(mcmcdf = fmod, datdf = fdf, sex = "FEMALE") %>%
  dplyr::select(-Index, -DoY, -contains("Phenophase"), -Date, -TreeUnique, -mean_temp, -OrchardID, -contains("forcing"), -TreeID, -contains("_mean"), -contains("sigma_")) %>%
  distinct()
#mpardf <- build_par_df(mcmcdf = mmod, datdf = mdf, sex = "MALE")


# calculate apc for site 1
fstart <- mutate(pardf, fstart_real = calcstageforcing(p=0.2, 
                                                       beta= beta + b_site + b_prov + b_clone + b_year,
                                                       kappa = kappa1)) 
#add a column to pardf with the estimated fstart from the model


# calculate fstart as if all trees had grown at site 1
site1 <- select(fmod, "b_site[1]", iter) # THIS IS WRONG. See what happens in GetComparisonDFFromPairs.
# I wonder if I could use GetComparisonDFFfrompairs if I created 2 pardfs and then joined them, one with original and one with "B" from the pairs. 
# Then I could use the predictfunction. I really don't understand how that function works though.

# Try merging pardf with pairs, changing the .B, and then merging again with pairs.
# Then I should be able to use the predict function for yhat1 and yhat2 (GetComparisonDFFFromPairs)

pairsWithParams <- left_join(pairs,pardf)  # adds 1 row to pairs for every sample, so dim should be nrow(fmod) * nrow(pairs)
nrow(pairsWithParams) == nrow(fmod)*nrow(pairs)



# merge with pardf (should not add any rows)
# readd B suffix to anything that doesn't end in .A
# remove .A suffix

# add a temporary suffix .A to nonB columns EXCEPT for iter
foo <- pairsWithParams %>%
  select(-ends_with(".B"), -iter)
colnames(foo) <- paste(colnames(foo), ".A", sep="")

# remove B suffix
foob <- pairsWithParams %>%
  select(ends_with(".B"), -iter)
colnames(foob) <- str_extract(colnames(foob), "\\w+")

#recombine and add params
foo <- cbind(foo, foob)
foo <- left_join(foo, pardf) # join after combining to maintain order and weights

#correct colnames (re-add B suffix)
iter <- foo$iter
Bs <- select(foo, -iter, -ends_with(".A"))
colnames(Bs) <- paste(colnames(Bs), ".B", sep="")
  # delete A suffix
As <- select(foo, ends_with(".A"))
colnames(As) <- str_extract(colnames(As), "\\w+")
foo <- data.frame(Bs, iter, As) %>%
  select(contains("Original"), contains("ID"), contains("b_"), contains("beta"), contains("kappa1"), Weight, iter) %>% #can drop some of those other cols earlier i think
  select(OriginalRowNumber.B, beta.B, SiteID.B, b_site.B, ProvenanceID.B, b_prov.B, YearID.B, b_year.B, CloneID.B, b_clone.B, kappa1.B,
         OriginalRowNumber, beta, SiteID, b_site, ProvenanceID, b_prov, YearID, b_year, CloneID, b_clone, kappa1, iter, Weight) #arrange cols

calc_fstart <- function(df) {
  forcing <- with(df, {
    (logit(0.2) + kappa1)/(beta + b_site + b_prov + b_clone + b_year)
  })
  return(forcing)
}

hist(calc_fstart(foo))

GetComparisonDFFromPairs.function <- function(predictionFunction, pairs, u, v) {
  uNew <- paste(u,".B",sep="")
  pairs$yHat1 <- predictionFunction(pairs)
  pairsNew <- structure(pairs[,c(v,uNew)], names=c(v,u)) #renaming u in pairsNew so we can call predictionFunction
  pairs$yHat2 <- predictionFunction(pairsNew)  
  return(pairs)
}


compdf <- GetComparisonDFFromPairs.function(calc_fstart, foo, u=c("SiteID", "b_site"), v=c("beta", "ProvenanceID", "b_prov", "CloneID", "b_clone", "YearID", "b_year", "kappa1"))

ggplot(compdf, aes(x=yHat1)) +
  geom_histogram() +
  geom_histogram(aes(x=yHat2), color="blue", alpha=0.5)

ComputeApcFromPairs <- function(predictionFunction, pairs, u, v, absolute=FALSE, impact=FALSE) {
  uNew <- paste(u,".B",sep="")
  ComparisonDF <- GetComparisonDFFromPairs(predictionFunction, pairs, u, v)
  absoluteOrIdentity <- if (absolute) abs else identity
  uDiff <- ComparisonDF[[uNew]] - ComparisonDF[[u]]
  denom <- if (impact) sum(ComparisonDF$Weight) else sum(ComparisonDF$Weight * uDiff * sign(uDiff))
  Apc <- sum(absoluteOrIdentity(ComparisonDF$Weight * (ComparisonDF$yHat2 - ComparisonDF$yHat1) * sign(uDiff))) / denom
  return(Apc)
}

# This needs to account for samples and it doesn't right now
# Summation order don't matter https://www.math.ubc.ca/~feldman/m321/twosum.pdf

num <- (sum(compdf$Weight * (compdf$yHat2 - compdf$yHat1)))^2
denom <- sum(compdf$Weight) * nrow(fmod)
apc <- (num/denom)^1/2

# calculate standard error
seframe <- compdf %>%
  group_by(iter) %>%
  mutate(numnk = Weight*(yHat2 - yHat1)^2) %>%
  summarise(apciter = sum(numnk)/sum(Weight))

t1 <- 1/(2*apc)
t2 <- 1/(nrow(seframe)-1)
t3 <- sum(seframe$apciter - apc^2)^2
seApc <- t1*sqrt(t2*t3)


############# garbage below


# compdf <- select(fstart, -contains("Site")) %>%
#   left_join(site1, by="iter") %>%
#   rename(b_site = 'b_site[1]') %>%
#   mutate(fstart_site1 = calcstageforcing(p=0.2, 
#                                          beta= beta + b_site + b_prov + b_clone + b_year,
#                                          kappa = kappa1))
# 

fstart$fstart_site1 <- compdf$fstart_site1 # add expected value for u1 to fstart df


# IF I'M GOING WRONG, THIS IS LIKELY WHERE. I should be merging in only comparisons between site 1 data and all other data
# So maybe that means doing a different kind of join? But I think this is right. 
compdfwithweights <- full_join(fstart, pairs) # requires A LOT of RAM # get the weight for each comparison

#signed (identity)
difference <- compdfwithweights$fstart_site1 - compdfwithweights$fstart_real
diffweighted <- difference * compdfwithweights$Weight
num <- sum(diffweighted)
denom <- sum(compdfwithweights$Weight) * nrow(fmod)
apc <- num/denom
print(apc)

#absolute
difference <- compdfwithweights$fstart_site1 - compdfwithweights$fstart_real
diffabs <- abs(difference)
diffweighted <- diffabs * compdfwithweights$Weight
num <- sum(diffweighted)
denom <- sum(compdfwithweights$Weight) * nrow(fmod)
apc <- num/denom
print(apc)

#gelman
difference <- compdfwithweights$fstart_site1 - compdfwithweights$fstart_real
diffsq <- difference^2
diffweighted <- diffsq * compdfwithweights$Weight
num <- sum(diffweighted)
denom <- sum(compdfwithweights$Weight) * nrow(fmod)
apc <- sqrt(num/denom)

end_time <- Sys.time()

end_time-start_time

# next up is calculating APC for all sites - turning the above into a function

#################
# ( (data * data) - data) * samples

##################3

allcomb <- expand.grid(udf$IndSexGroup, udf$IndSexGroup) %>%
  filter(Var1 != Var2) %>% # get a list of all row combinations
  arrange(Var1) %>%
  rename(IndSexGroup = Var2) %>%
  left_join(udf) %>%
  rename(IndSexGroup.B = IndSexGroup) %>%
  rename(IndSexGroup = Var1) %>%
  left_join(udf, suffix=c("", ".B"), by="IndSexGroup") %>%
  arrange(Sex, Sex.B, Site, Site.B, SPU_Name, SPU_Name.B, Clone, Year)






# Now to calculate the predictive comparison

# Expected value real

pardf <- pardf %>% mutate(fstart_real = calcstageforcing(p=0.2, beta = beta+b_clone+b_site+b_prov+b_year, kappa=kappa1))

ComputeApcFromPairs <- function(predictionFunction, pairs, u, v, absolute=FALSE, impact=FALSE) {
  uNew <- paste(u,".B",sep="")
  ComparisonDF <- GetComparisonDFFromPairs(predictionFunction, pairs, u, v)
  absoluteOrIdentity <- if (absolute) abs else identity
  uDiff <- ComparisonDF[[uNew]] - ComparisonDF[[u]]
  denom <- if (impact) sum(ComparisonDF$Weight) else sum(ComparisonDF$Weight * uDiff * sign(uDiff))
  Apc <- sum(absoluteOrIdentity(ComparisonDF$Weight * (ComparisonDF$yHat2 - ComparisonDF$yHat1) * sign(uDiff))) / denom
  return(Apc)
}

calc_fstart <- function(p=0.2, beta=betas, kappa, model, pairs) {
  # calc expected value
  b_prov <- select(matches(paste("b_prov", pairs$ProvenanceID[1], sep="\\[")))[1] #choose the first sample and the first row
  prob <- (logit(p) + kappa)/beta
  return(prob)
}

GetComparisonDFFromPairs.function <- function(predictionFunction, pairs, u, v) {
  uNew <- paste(u,".B",sep="")
  pairs$yHat1 <- predictionFunction(pairs)
  pairsNew <- structure(pairs[,c(v,uNew)], names=c(v,u)) #renaming u in pairsNew so we can call predictionFunction
  pairs$yHat2 <- predictionFunction(pairsNew)  
  return(pairs)
}



# now add weights

# calculate the start of forcing for all model configurations for the first "observation"
udf[1,]

# select all site betas (model configs) from model
site1 <- select(fmod, `b_site[1]`)
site2 <- select(fmod, `b_site[2]`)
site3 <- select(fmod, `b_site[3]`)

#select all prov, clone, and year betas (model configs) that for the provenance, clone, and year 
# in the first row of the dataset (at one v)
prov <- select(fmod, matches(paste("b_prov", fdf$ProvenanceID[1], sep="\\[")))
clone <- select(fmod, matches(paste("b_clone", fdf$CloneID[1], sep="\\[")))
year <- select(fmod, matches(paste("b_year", fdf$YearID[1], sep="\\[")))
#select the overall beta
beta <- select(fmod, beta)

# calculate the total beta for each site and 1 v
betas1 <- beta + site1 + prov + clone + year
betas2 <- beta + site2 + prov + clone + year
betas3 <- beta + site3 + prov + clone + year

# select the kappas from the model
kappa1 <- select(fmod, `kappa[1]`)
kappa2 <- select(fmod, `kappa[2]`)



# calculate the expected value for starting for each site at v1
start1 <- calcstageforcing(p=0.2, beta=betas1, kappa=kappa1)
start2 <- calcstageforcing(p=0.2, beta=betas2, kappa=kappa1)
start3 <- calcstageforcing(p=0.2, beta=betas3, kappa=kappa1)
end <- calcstageforcing(p=0.8, kappa=kappa2)

# calculate the squared difference in expected value for starting for each site at v1
diff1 <- (start1-start2)^2
diff2 <- (start1-start3)^2
diff3 <- (start2-start3)^2
diff <- cbind(diff1, diff2, diff3) #combine into a dataframe

diffsum <- colSums(diff) # sum over model configs
usum <- sum(diffsum)# sum over param combinations
denom <- nrow(udf) * nrow(fmod) * ncol(diff)

siteapc <- usum/denom

# WEIGHTS ###################


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

# Site effects - create the "v" dataframe. Select all data that isn't the site.
udf <- unique_grouper(fdf, mdf) %>%
  #filter(Sex=="FEMALE") %>%
  distinct() %>%
  select("SPU_Name", "Clone", "Year", "Sex") # v

# create weights dataframe
weights <- matrix(nrow=nrow(udf), ncol = nrow(udf)) #empty matrix for weights with a row and column 
# for every row in udf (for every v)

# calculate a matrix where the upper half holds every row-col comparison
for (i in 1:nrow(udf)) { # This calculation is slow af
  for (j in i:nrow(udf)) {
    weights[i,j] <- sum(udf[i,] %in% udf[j,])
  }
}

# now make the matrix look like a times table where every row-row comparison is available 
# starting from the row or the column
#make weights matrix symmetric. fill out the weights matrix by flipping the upper triangle matrix across
# the diagona
weights[is.na(weights)] <- 0 # NA bottom half should be 0s so that addition works later
tweights <- t(weights)
diag(tweights) <- 0 # replace diagonal so self doesn't get counted twice
# calculate weights for bottom triangle
weights <- (weights + tweights)/ncol(udf)
#STOP

# normalize weights
wdf <- as.data.frame(weights)
wdf$totweights <- rowSums(weights) # total weight for each v
#wdf <- wdf/max(wdf$totweights) # dataframe with weights (normalized) for each comparison of v's that occur in the dataset

wdf <- data.frame(OriginalRowNumber = rownames(wdf), Weight = wdf$totweights)

# test : are all totweights == 1?

# Predcomps ###############

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