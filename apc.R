# Average Predictive comparisons
# references Gelman & Pardoe 2007 + dchudz predcomps

source('phenology_functions.R')

library(tidyverse)
library(gtools)
library(assertthat)

s = 1000
u = "b_clone"
uid = "CloneID"
v = c("b_site", "b_prov", "b_year")
vid = c("SiteID", "ProvenanceID", "YearID")

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

calc_fstart <- function(df) {
  forcing <- with(df, {
    (logit(0.2) + kappa1)/(beta + b_site + b_prov + b_clone + b_year)
  })
  return(forcing)
}

GetComparisonDFFromPairs.function <- function(predictionFunction, pairs, u, v) {
  uNew <- paste(u,".B",sep="")
  pairs$yHat1 <- predictionFunction(pairs)
  pairsNew <- structure(pairs[,c(v,uNew)], names=c(v,u)) #renaming u in pairsNew so we can call predictionFunction
  pairs$yHat2 <- predictionFunction(pairsNew)  
  return(pairs)
}


#################

start_time <- Sys.time()

fmod <- readRDS("2019-10-28phenologyFEMALE.rds") %>%
  as.data.frame() %>%
  select(-contains("state_rep")) %>%
  #sample_frac(0.05)
  sample_n(s) 
fmod$iter <- 1:nrow(fmod)

mmod <- readRDS("2019-10-28phenologyMALE.rds") %>%
  as.data.frame() %>%
  select(-contains("state_rep")) %>%
  # sample_frac(0.05)
  sample_n(s)
mmod$iter <- 1:nrow(mmod)


# original data 
phendf <- read_data(slim = FALSE) 

# identify combinations of effects that actually occur

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")


udf <- unique_grouper(fdf, mdf) %>%
  filter(Sex=="FEMALE") 

# pairs of rows and their weights
#opairs <- GetPairs(X=udf, u="SiteID", v=c("ProvenanceID", "CloneID", "YearID"))
opairs <- GetPairs(X=udf, u=uid, v=vid) #General

assert_that(nrow(opairs) == (nrow(udf) * nrow(udf)) - nrow(udf), msg="Pairs dataframe is the wrong size")

# I need to rethink these functions to be smaller because (now that I'm not doing ppc) I don't need to keep every observation
# pardf is a dataframe with nrow(fmod) x nrow(udf) rows that contains identifying information for data 
# and the parameter values associated with the data
pardf <- build_par_df(mcmcdf = fmod, datdf = fdf, sex = "FEMALE") %>% # since I'm removing so much after, I think I can make a better build function here
  dplyr::select(-Index, -DoY, -contains("Phenophase"), -Date, -TreeUnique, -mean_temp, -contains("Orchard"), 
                -contains("forcing"), -TreeID, -contains("_mean"), -contains("sigma_"), -contains("kappa2"), -Sex,
                -Site, -SPU_Name, -Clone, -Year) %>%
  distinct()
#BUILD_PAR_DF IS SLOW SLOW OMG FIX IT
assert_that(nrow(pardf) == nrow(fmod)*nrow(udf))

# it might be faster to fold this into for loop, or at least allow more samples nrow(fmod)*nrow(opairs) will get big fast
#pairsWithParams <- left_join(opairs,pardf)  # adds 1 row to pairs for every sample, so dim should be nrow(fmod) * nrow(pairs)
#assert_that(nrow(pairsWithParams) == nrow(fmod)*nrow(opairs), msg = "pairsWithParams is the wrong size")

# calculate average predictive comparison for each sample
apc_sample <- data.frame(iter=rep(NA, nrow(fmod)), pc=NA)

# this loop could be parallelized if you were clever
for (i in 1:nrow(fmod)) {
  print(paste("Working on sample", i, "of", nrow(fmod)))
  apc_sample$iter[i] <- fmod$iter[i]
  
  # Join the pairs with parameter values from the model so I can calculate expected values (yHats)
    # Get parameter values for rows without B suffix (u_i)
  pardf_iteri <- filter(pardf, iter==i) # grab parameter values from the first model iteration/sample
  pairsWithParams <- left_join(opairs, pardf_iteri, by = c(uid, vid))
  assert_that(nrow(pairsWithParams) == nrow(opairs), msg = "pairsWithParams is the wrong size")
  
  tpairs <- filter(pairsWithParams, iter==i)
  assert_that(nrow(tpairs) == (nrow(udf)*nrow(udf)) - nrow(udf), msg = "You've selected too many rows from pairsWithParams")
  #tpardf <- filter(pardf, iter==i) %>%
   # select( contains("ID"), iter, contains("site") ) # NEEDS GENERALIZATION
  tpardf <- pardf_iteri %>% #GENERAL
    select( contains("ID"), iter, u)
  
    # Format the pairs df and add u parameter values for comparison rows B (u^(k)) 
    # add a temporary suffix .A to nonB columns EXCEPT for iter
  As <- tpairs %>%
    select(-ends_with(".B"), -iter)
  colnames(As) <- paste(colnames(As), ".A", sep="")
  
    # remove B suffix
  Bs <- tpairs %>%
    select(ends_with(".B"), -iter)
  colnames(Bs) <- str_extract(colnames(Bs), "\\w+")
  
    #recombine and add params
  tpairs <- cbind(As, Bs) %>%
    left_join(tpardf, by=c(uid,vid)) # join after combining to maintain order and weights
  
  #correct colnames (re-add B suffix)
  iter <- tpairs$iter
  assert_that(length(unique(tpairs$iter))==1, msg="you should only have 1 sample from the model, but you have too many iterations")
  Bs <- select(tpairs, -iter, -ends_with(".A"))
  colnames(Bs) <- paste(colnames(Bs), ".B", sep="")
  # delete A suffix
  As <- select(tpairs, ends_with(".A"))
  colnames(As) <- str_extract(colnames(As), "\\w+")
  pairs <- data.frame(Bs, iter, As) %>%
    select(contains("Original"), contains("ID"), contains("b_"), contains("beta"), contains("kappa1"), Weight, iter) %>% #can drop some of those other cols earlier i think
   #arrange cols
    select(OriginalRowNumber.B, paste(uid, ".B", sep=""), paste(u, ".B", sep=""),
           OriginalRowNumber, beta, kappa1, u, uid, v, vid, iter, Weight) #arrange cols
           
  assert_that(nrow(pairs)==nrow(tpairs), msg = "why isn't pairsWithParams the same size as pairs?")
  
  
  compdf <- GetComparisonDFFromPairs.function(calc_fstart, pairs, u=c(u, uid), v=c("beta", "kappa1", v, vid))
  
  # Summation order don't matter https://www.math.ubc.ca/~feldman/m321/twosum.pdf
  
  num <- sum(compdf$Weight * (compdf$yHat2 - compdf$yHat1)^2)
  denom <- sum(compdf$Weight) #* nrow(fmod)
  
  apc_sample$pc[i] <- num/denom
}

# calculate apc by summing over all samples and averaging
apc <- sqrt(sum(apc_sample$pc)/nrow(apc_sample))
print(apc)
# calculate standard error

t1 <- 1/(2*apc)
t2 <- 1/(s-1)
diff <- apc_sample$pc^2 - apc^2
diffsq <- diff^2
t3 <- sum(diffsq)
seApc <- t1*sqrt(t2*t3)
print(seApc)



end_time <- Sys.time()
end_time-start_time
# print(apc)
# print(seApc)
# 

