# Average Predictive comparisons
# references Gelman & Pardoe 2007 + dchudz predcomps

# tidybayes might help here with some of the weird pardf stuff
library(tidyverse)
library(gtools)
library(assertthat)

source('phenology_functions.R')

fmod_raw <- readRDS("2019-10-28phenologyFEMALE.rds") %>%
  as.data.frame() %>%
  select(-contains("state_rep"))

mmod_raw <- readRDS("2019-10-28phenologyMALE.rds") %>%
  as.data.frame() %>%
  select(-contains("state_rep"))

# GLOBALS ############

forcingtype = "scaled_ristos"
# 7000
samplefrac=0.1
r=20


########## FUNCTIONS ####################

# this function is based on a similar function in predcomps. Considering a dataset X with variable
# of interest uid and other variables vid, this function creates a dataframe where each row is a pairwise
# comparison to another row with a weight.
GetPairs <- function(X, uid, vid,
                     numForTransitionStart = NULL,
                     numForTransitionEnd = NULL,
                     onlyIncludeNearestN = NULL) {

  assertthat::assert_that(length(uid) == 1) # make sure we have exactly 1 input var of interest
  for (columnName in c(uid,vid)) {
    assertthat::assert_that(columnName %in% names(X))
    columnClass <- class(X[[columnName]])
    # if (!(columnClass) %in% c("integer", "numeric")) {
    #   stop(sprintf("Sorry, column %s is of class %s. I can only deal with integer and numeric types for now.", columnName, columnClass))
    # }
  }

  if (!is.null(numForTransitionStart)) {
    X1 <- X[sample.int(nrow(X), size=numForTransitionStart), c(vid,uid)]
  } else {
    X1 <- X[c(vid,uid)] # all v with u
  }

  if (!is.null(numForTransitionEnd)) {
    X2 <- X[sample.int(nrow(X), size=numForTransitionEnd), c(vid,uid)]
  } else {
    X2 <- X[c(vid,uid)] # all v with u
  }

  X1$OriginalRowNumber <- 1:nrow(X1)
  X2$OriginalRowNumber.B <- 1:nrow(X2)

  vMatrix1 <- as.matrix(X1[,vid])
  vMatrix2 <- as.matrix(X2[,vid])

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
  weights <- (weights + tweights)/length(vid)

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


# calculate standard error where apc is the apc (1 number) and apc_sample is the per sample apc
calc_seApc <- function(apc, apc_sample, n) {
  t1 <- 1/(2*apc)
  t2 <- 1/(n-1)
  #diff <- apply(apc_sample, MARGIN=1, FUN=function(x) (x^2 - apc^2)) # RMS
  diff <- apply(apc_sample, MARGIN=1, FUN=function(x) (x - apc))
  diffsq <- diff^2
  t3 <- rowSums(diffsq)
  seApc <- t1*sqrt(t2*t3)
  return(seApc)
}


# REWRITE WITH MU AND ALPHAS - OR ACTUALLY THIS IS JUST YPPC?
# calculate expected value of fstart
calc_fstart <- function(df) {
  forcing <- with(df, {
    (logit(0.2) + kappa1)/(beta + b_site + b_prov + b_clone + b_year)
  })
  return(forcing)
}

# make a dataframe that adds expected values to the pairs dataframe
GetComparisonDFFromPairs.function <- function(predictionFunction, pairs, u, v) {
  uNew <- paste(u,".B",sep="")
  pairs$yHat1 <- predictionFunction(pairs)
  pairsNew <- structure(pairs[,c(v,uNew)], names=c(v,u)) #renaming u in pairsNew so we can call predictionFunction
  pairs$yHat2 <- predictionFunction(pairsNew)
  return(pairs)
}
#
# calc_sample_apc <- function(df, yhat1, yhat2) {
#   num <- sum(df$Weight * (df[yhat2] - df[yhat1])^2)
#   denom <- sum(df$Weight) #* nrow(fmod)
#   sapc <- num/denom
#   return(sapc)
# }

match_fu_to_DoY <- function(climdf, yhat1 = compdf$yHat1, yhat2=compdf$yHat2) {
  yHat1 <- climdf$DoY[findInterval(yhat1,
                                   climdf$sum_forcing) + 1]
  yHat2 <- climdf$DoY[findInterval(yhat2,
                                   climdf$sum_forcing) + 1]
  df <- data.frame(yHat1, yHat2)
  #colnames(df) <- paste(colnames(df), nameappend, sep="_")
  return(df)
}

# calculate the per-sample apc. Can be calculated either with absolute differences ("type = "absolute")
# or signed differences (type = "signed")
calc_sample_apc <- function(weight, yhat1, yhat2, type = "absolute") {
  #num <- sum(weight * (yhat2 - yhat1)^2) # squared difference for RMS calculation
  if (type=="absolute") {
    num <- sum(weight * abs(yhat2 - yhat1)) # absolute difference
  } else {
    num <- sum(weight * (yhat2-yhat1)) # signed difference
    }
  denom <- sum(weight) #* nrow(fmod)
  sapc <- num/denom
  return(sapc)
}

calc_doy_apc <- function(colname, apc_sample=apc_sample, df=compdf, type="absolute") {
  yhats <- select(df, contains(colname), Weight)
  sapc <- calc_sample_apc(yhats$Weight, yhats[,1], yhats[,2], type=type)
}


# calculate apc from sample apcs

calc_apc <- function(vec) {
 # apc <- sqrt(sum(vec)/length(vec)) #RMS
  apc <- sum(vec)/length(vec) # absolute or signed comparison
  return(apc)
}

calc_apc_of_variable <- function(u, uid, v, vid, phendata = udf,
                                 model_params = pardf, n=nrow(fmod), climlist=climlist) {

  # pairs of rows and their weights
  opairs <- GetPairs(X=phendata, uid=uid, vid=vid) #slow step if dataset is large
  assert_that(nrow(opairs) == (nrow(phendata) * nrow(phendata)) - nrow(phendata), msg="Pairs dataframe is the wrong size")

  # calculate average predictive comparison for each sample
  names_absolute <- names(climlist)
  names_signed <- paste(names(climlist), "_signed", sep="")
  names <- c(names_absolute, names_signed)
  apc_sample <- data.frame(iter=NA, forcing=NA, forcing_signed=NA, setNames(rep(list(NA), length(names)), names))

  # this loop could be parallelized if you were clever
  for (i in 1:n) {
    print(paste("Working on sample", i, "of", n))

    # Join the pairs with parameter values from the model so I can calculate expected values (yHats)
    # Get parameter values for rows without B suffix (u_i)
    pardf_iteri <- filter(model_params, iter==i) # grab parameter values from the first model iteration/sample
    pairsWithParams <- left_join(opairs, pardf_iteri, by = c(uid, vid))

    assert_that(nrow(pairsWithParams) == nrow(opairs), msg = "pairsWithParams is the wrong size")
    assert_that(nrow(pairsWithParams) == (nrow(phendata)*nrow(phendata)) - nrow(phendata), msg = "You've selected too many rows from pairsWithParams")

    tpardf <- pardf_iteri %>% #GENERAL
      select( contains("ID"), iter, u)

    # Format the pairs df and add u parameter values for comparison rows B (u^(k))

    # add a temporary suffix .A to nonB columns EXCEPT for iter
    As <- pairsWithParams %>%
      select(-ends_with(".B"), -iter)
    colnames(As) <- paste(colnames(As), ".A", sep="")

    # remove B suffix
    Bs <- pairsWithParams %>%
      select(ends_with(".B"), -iter)
    colnames(Bs) <- str_extract(colnames(Bs), "\\w+")

    #recombine and add params
    pairsWithParams <- cbind(As, Bs) %>%
      left_join(tpardf, by=c(uid,vid)) # join after combining to maintain order and weights

    #correct colnames (re-add B suffix)
    iter <- pairsWithParams$iter
    assert_that(length(unique(pairsWithParams$iter))==1, msg="you should only have 1 sample from the model, but you have too many iterations")
    Bs <- select(pairsWithParams, -iter, -ends_with(".A"))
    colnames(Bs) <- paste(colnames(Bs), ".B", sep="")
    # delete A suffix
    As <- select(pairsWithParams, ends_with(".A"))
    colnames(As) <- str_extract(colnames(As), "\\w+")
    pairs <- data.frame(Bs, iter, As) %>%
      select(contains("Original"), contains("ID"), contains("b_"), contains("beta"), contains("kappa1"), Weight, iter) %>% #can drop some of those other cols earlier i think
      #arrange cols
      select(OriginalRowNumber.B, paste(uid, ".B", sep=""), paste(u, ".B", sep=""),
             OriginalRowNumber, beta, kappa1, u, uid, v, vid, iter, Weight) #arrange cols

    assert_that(nrow(pairs)==nrow(pairsWithParams), msg = "why isn't pairsWithParams the same size as pairs?")


    compdf <- GetComparisonDFFromPairs.function(calc_fstart, pairs, u=c(u, uid), v=c("beta", "kappa1", v, vid))

    # translate fstart into DoY
    climlist <- climlist[order(names(climlist))]

    doy_yhats <- sapply(climlist, FUN = match_fu_to_DoY, USE.NAMES = TRUE, yhat1=compdf$yHat1, yhat2=compdf$yHat2)
    names(doy_yhats) <- paste(rep(attributes(doy_yhats)$dimnames[[1]],length(climlist)),
                                 sort(rep(attributes(doy_yhats)$dimnames[[2]],2)), sep="_") #name with the correct yhat # and descriptor
    doy_yhats <- do.call(cbind.data.frame, doy_yhats)

    # add doy yhats to comparison df
    compdf <- cbind(compdf, doy_yhats)

    # calculate sample apc for fstart
    sapc <- calc_sample_apc(compdf$Weight, compdf$yHat1, compdf$yHat2, type = "absolute")
    sapc_signed <- calc_sample_apc(compdf$Weight, compdf$yHat1, compdf$yHat2, type = "signed")

    # calculate sample apc for doy in various years
    doy_apc <- sapply(names(climlist), FUN = calc_doy_apc, USE.NAMES=TRUE, df=compdf, type="absolute")
    doy_apc_signed <- sapply(names(climlist), FUN = calc_doy_apc, USE.NAMES=TRUE, df=compdf, type="signed")

    apc_sample[i,] <- c(iter=i, sapc=sapc, sapc_signed, doy_apc, doy_apc_signed)


    # apc_sample$sapc_cold[i] <- calc_sample_apc(compdf$Weight, compdf$yHat1_cold, compdf$yHat2_cold)
    # apc_sample$sapc_hot[i] <- calc_sample_apc(compdf$Weight, compdf$yHat1_hot, compdf$yHat2_hot)
    # apc_sample$sapc_med[i] <- calc_sample_apc(compdf$Weight, compdf$yHat1_med, compdf$yHat2_med)
    # apc_sample$sapc_median[i] <- calc_sample_apc(compdf$Weight, compdf$yHat1_median, compdf$yHat2_median)
  }

  apcs <- apply(apc_sample[,-1], MARGIN=2, FUN=calc_apc)
  apcs <- data.frame(comparison = names(apcs), apcs )


  # apc <- calc_apc(apc_sample, "sapc")
  # coldapc <- calc_apc(apc_sample, "sapc_cold")
  # medapc <- calc_apc(apc_sample, "sapc_med")
  # medianapc <- calc_apc(apc_sample, "sapc_median")
  # hotapc <- calc_apc(apc_sample, "sapc_hot")
  # apcs <- data.frame(comparison = c("apc", "cold_apc", "med_apc", "median_apc", "hot_apc"),
  #                    apc=c(apc, coldapc, medapc, medianapc, hotapc))
  apcs$se <- calc_seApc(apcs$apc, apc_sample[,-1], n)

  return(apcs)

}

# read in data #############

## Model data #####


fmod <- sample_frac(fmod_raw, samplefrac)
fmod$iter <- 1:nrow(fmod)

mmod <- sample_frac(mmod_raw, samplefrac)
mmod$iter <- 1:nrow(mmod)

## Phenology data ####
phendf <- read_data(slim = FALSE)

# identify combinations of effects that actually occur

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")

fudf <- unique_grouper(fdf, mdf) %>%
  filter(Sex=="FEMALE")

mudf <- unique_grouper(fdf, mdf) %>%
  filter(Sex=="MALE")

## Climate data ####
clim <- read.csv("data/all_clim_PCIC.csv", header=TRUE, stringsAsFactors=FALSE) %>%
  filter(forcing_type==forcingtype)

### choose r climate timeseries from coldest to hottest

clim$siteyear <- paste(clim$Site, clim$Year, sep='') # index clim and tpars by Site-Year

# choose cold, hot, and medium years.

phendf2 <- filter(phendf, Phenophase_Derived==2)

climsort <- filter(clim, DoY == round(mean(phendf2$DoY))) %>%
  arrange(sum_forcing) %>%
  mutate(coldtohot=order(sum_forcing))
#
# coldyear <- filter(clim, siteyear == climsort$siteyear[1]) %>%
#   arrange(DoY)
# hotyear <- filter(clim, siteyear == climsort$siteyear[nrow(climsort)]) %>%
#   arrange(DoY)

#ryears <- sample(unique(climsort$siteyear), r) #only run this line to regenerate years
#
# randomyears <- list()
# for (i in 1:length(ryears)) {
#   randomyears[[i]] <- filter(clim, siteyear==ryears[i]) %>%
#     select(DoY, sum_forcing, siteyear) %>%
#     arrange(DoY)
# }
# names(randomyears) <- ryears
# climlist <- list(hot=hotyear, cold=coldyear) %>%
#   append(randomyears)

spaceseq <- round(seq(from=1, to =length(unique(climsort$siteyear)), length.out = r))
syears <- unique(climsort$siteyear)[spaceseq]

spacedyears <- list()
for (i in 1:length(syears)) {
  spacedyears[[i]] <- filter(clim, siteyear==syears[i]) %>%
    select(DoY, sum_forcing, siteyear) %>%
    arrange(DoY)
}

names(spacedyears) <- syears
climlist <- spacedyears

# graph different accumulated forcing timeseries

apcclimdf <- dplyr::bind_rows(climlist) %>%
  left_join(select(climsort, coldtohot, siteyear))

ggplot(clim, aes(x=DoY, y=sum_forcing, group=siteyear)) +
  geom_line(color="darkgray") +
  #geom_line(data=apcclimdf, aes(x=DoY, y=sum_forcing, color=coldtohot), alpha=0.5) +
  xlim(c(0,180)) +
  ylim(c(0,35)) +
  geom_vline(xintercept = mean(phendf2$DoY), color="dark gray") +
 # geom_jitter(data=filter(phendf, Phenophase_Derived==2), aes(x=DoY, y=sum_forcing),
           #  inherit.aes = FALSE, pch=2, alpha=0.1) +
  ggtitle("Forcing accumulation in 30 of 105 accumulated forcing time series",
          subtitle = "with stage 2 phenology") +
  xlab("Day of Year") +
  ylab("Accumulated forcing") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_viridis_c("cividis") +
  geom_line(data = filter(clim, siteyear %in% weirdyears), aes(x=DoY, y=sum_forcing, group=siteyear), color="red")


weirdyears <- c("Kalamalka2003", "Tolko2009", "Vernon2011")



######################

# pardf is a dataframe with nrow(fmod) x nrow(udf) rows that contains identifying information for data
# and the parameter values associated with the data
fpardf <- build_par_df(mcmcdf = fmod, datdf = fudf, sex = "FEMALE") %>%
  dplyr::select(-contains("Ind"), -contains("_mean"), -contains("sigma_"), -contains("kappa2"), -Sex,
             -Site, -SPU_Name, -Clone, -Year)
assert_that(nrow(fpardf) == nrow(fmod)*nrow(fudf))

mpardf <- build_par_df(mcmcdf = mmod, datdf = mudf, sex = "MALE") %>%
  dplyr::select(-contains("Ind"), -contains("_mean"), -contains("sigma_"), -contains("kappa2"), -Sex,
                -Site, -SPU_Name, -Clone, -Year)
assert_that(nrow(mpardf) == nrow(mmod)*nrow(mudf))

# you can source all the apc calculation files OR run each one as a local job, which is much faster
# could probably parallelize somehow, but I don't know
# move all apc files to apc folder

# store apcs ############3

#collect apcs for all sexes and u's and save to disk because that took a long-ass time to calculate
apcs_male <- list(site=apc_male_site, prov=apc_male_prov, year=apc_male_year, clone=apc_male_clone)
apcs_female <- list(site=apc_female_site, prov=apc_female_prov, year=apc_female_year, clone=apc_female_clone)
mapcsdf <- map2_df(apcs_male, names(apcs_male), ~ mutate(.x, ID = .y)) %>%
  arrange(ID, comparison) %>%
  mutate(sex="MALE")
fapcsdf <- map2_df(apcs_female, names(apcs_female), ~ mutate(.x, ID = .y)) %>%
  arrange(ID, comparison) %>%
  mutate(sex="FEMALE")
apcsdf <- rbind(mapcsdf, fapcsdf)

write.csv(apcsdf, "apc_values.csv", row.names = FALSE)

# visualize apc ###############


site_apct <- site_apct %>%
  separate(comparison, into = c("comparison", "apc_type"), sep="_")
site_apct[which(is.na(site_apct$apc_type)),]$apc_type <- "absolute"

site_apc_err <- site_apct %>%
  mutate(upper = apcs + se) %>%
  mutate(lower= apcs - se) %>%
  pivot_longer(cols = c(upper, lower), names_to = "side", values_to = "errorbound")

apc_err <- apcsdf %>%
  mutate(upper = apcs + se) %>%
  mutate(lower= apcs - se) %>%
  pivot_longer(cols = c(upper, lower), names_to = "side", values_to = "errorbound")

library(ggExtra)
ggplot(apcsdf, aes(x=apcs, y=comparison)) +
  geom_point() +
  geom_line(data=apc_err, aes(x=errorbound, y=comparison)) +
  geom_vline(xintercept=0) +
  theme_bw() +
  removeGrid(y=TRUE, x=FALSE) +
  facet_grid(ID ~ .)

ggplot(site_apct, aes(x=apcs, y=comparison, color=apc_type)) +
  geom_point() +
  geom_line(data=site_apc_err, aes(x=errorbound, y=comparison)) +
  geom_vline(xintercept=0) +
  theme_bw() +
  removeGrid(y=TRUE, x=FALSE)
