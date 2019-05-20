# F50 analysis

# Difference between 50% transition by provenance

# LIBRARIES #################################
library(tidyr)
library(dplyr)
library(stringr)
library(rstan)
library(rethinking)
library(purrr)

# FUNCTIONS ###################################33
stanindexer <- function(df) {
    df$CloneID <- group_indices(df, Clone)
    df$OrchardID <- group_indices(df, Orchard)
    df$ProvenanceID <- group_indices(df, SPU_Name)
    df$SiteID <- group_indices(df, Site)
    df$YearID <- group_indices(df, Year)
    df$Tree <- group_indices(df,TreeID)
    return(df)
}

logistic2 <- function(x,b,c) {
    1/(1 + exp(-(b * (x-(c/b)))))
}

pdflogistic <- function(x,b,c) {
    num <- b * exp(-b*(x-(c/b)))
    denom <- (1 + exp(-(b * (x-(c/b)))))^2
    val <- num/denom
    return(val)
}

# MODEL AND ORIGINAL DATA ############
#model
ffit.stan <- readRDS("female_slopes.rds")
#mfit.stan <- readRDS("male_slopes.rds")

fmod <- as.data.frame(ffit.stan)
#mmod <- as.data.frame(ffit.stan)


#data
phenology_data <- read.csv("data/stan_input/phenology_heatsum.csv",
                           stringsAsFactors = FALSE, header = TRUE) %>%
    filter(forcing_type=="ristos")

clim <- read.csv("data/all_clim_PCIC.csv", stringsAsFactors = FALSE, header=TRUE) %>%
    filter(forcing_type=="ristos") %>%
    filter(DoY %in% c(90:190))

## provenance
SPU_dat <- read.csv("../research_phd/data/OrchardInfo/LodgepoleSPUs.csv",
                    header=TRUE, stringsAsFactors = FALSE) %>%
    dplyr::select(SPU_Name, Orchard)

# Data Processing ##################
# join provenance and phenology data

phendf <- phenology_data %>%
    na.omit()
phendf <- dplyr::left_join(phenology_data, SPU_dat) %>%
    unique()

# separate into male and female dataframes and turn factors into integers
fdf <- filter(phendf, Sex == "FEMALE")
fdf <- stanindexer(fdf)
fdf$groupID <- group_indices(fdf, SiteID, ProvenanceID, YearID, CloneID)
fdf$recordID <- group_indices(fdf, groupID, Date )
# mdf <- filter(phendf, Sex == "MALE")
# mdf <- stanindexer(mdf)

# identify combinations of effects & predictor (sum_forcing) that actually occur
ufdf <- fdf %>%
    select(groupID, recordID, Site, SiteID, SPU_Name, ProvenanceID, Clone, CloneID, Year, YearID, sum_forcing) %>%
    distinct()

# Build a dataframe where each row is a unique combination of parameters for each unique combination of forcing, site, clone, year, and provenance #########
#split by param

#draws <- base::sample(1:nrow(fmod), 300) # select a subset of model draws
draws <- c(1:nrow(fmod))
fmod_sampled <- fmod[draws,]

singledimpars <- fmod_sampled %>%
    mutate(draw=draws) %>%
    select(draw, beta, sigma_site, sigma_prov, sigma_clone, sigma_year, contains("kappa")) %>%
    rename(kappa1 = `kappa[1]`) %>%
    rename(kappa2 = `kappa[2]`)

tidypar <- function(stanframe, param, id) {

    #take a stan model dataframe and create a tidy dataframe for a given parameter. takes a dataframe, a string with the parameter name, and a string describing the param ID (e.g. par = "Site" and id="SiteID")

    #wide to long format for parameters
    par <- stanframe %>% select(contains(param)) %>%
        mutate(draw = draws) %>%
        gather(key = key, value = value, contains("b_")) %>%
        mutate(id = str_extract(key, "[0-9]{1,}"))
    colnames(par) <- c("draw", "name", param, id)
    par[,4] <- as.integer(par[,4])
    return(par)
}

kappa <- select(fmod_sampled, contains("kappa"))
siteb <- tidypar(fmod_sampled, "b_site", "SiteID")
provb <- tidypar(fmod_sampled, "b_prov", "ProvenanceID")
cloneb <- tidypar(fmod_sampled, "b_clone", "CloneID")
yearb <- tidypar(fmod_sampled, "b_year", "YearID")

clonemerge <- left_join(ufdf, cloneb)
provmerge <- left_join(ufdf, provb)
sitemerge <- left_join(ufdf, siteb)
yearmerge <- left_join(ufdf, yearb)

pardf <- data.frame(ufdf,
                    draw = clonemerge$draw,
                    b_clone = clonemerge$b_clone,
                    b_prov = provmerge$b_prov,
                    b_site = sitemerge$b_site,
                    b_year = yearmerge$b_year) %>%
    left_join(singledimpars)

#pardf is a dataframe with N rows of parameters for each sum_forcing, site, provenance, clone, and year combination that appear in the data, where N = number of draws from the posterior distribution

# Now simulate state predictions ##################

# calculate transitions, cdf, and pdf
pred_df <- pardf %>%
    mutate(betas = b_clone + b_prov + b_site + b_year + beta) %>%
    mutate(f501 = kappa1/betas) %>%
    mutate(f502 = kappa2/betas) %>%
    distinct()

# get all provenance combinations
provcombos <- crossing(1:6, 1:6)
provcombos <- unique(t(apply(provcombos, 1, sort)))

pred_df %>% group_by(ProvenanceID)
