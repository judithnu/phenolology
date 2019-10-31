# functions for prepping data for stan and model analysis

# Stan can only take consecutive integers for groups, so turn factors into consecutive integers
stanindexer <- function(df) {
    df$CloneID <- group_indices(df, Clone)
    df$OrchardID <- group_indices(df, Orchard)
    df$ProvenanceID <- group_indices(df, SPU_Name)
    df$SiteID <- group_indices(df, Site)
    df$YearID <- group_indices(df, Year)
    df$TreeID <- group_indices(df, TreeUnique)
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

read_data <- function(slim=TRUE) { #choose true if slimmed data and false for full dataset
    if (isTRUE(slim)) {
    phenology_data <- read.csv("data/phenology_heatsum.csv",
                               stringsAsFactors = FALSE, header = TRUE
    ) %>%
        filter(!(Year==2011 & Site=="KettleRiver"))
    } else {
        phenology_data <- read.csv("data/phenology_heatsum_all.csv",
                                   stringsAsFactors = FALSE, header = TRUE
        )
    }


    ## provenance
    SPU_dat <- read.csv("../phd/data/OrchardInfo/LodgepoleSPUs.csv",
                        header=TRUE, stringsAsFactors = FALSE) %>%
        dplyr::select(SPU_Name, Orchard)

    # join provenance and phenology data

    phendf <- phenology_data %>%
        na.omit() %>%
        left_join(SPU_dat) %>%
        unique()

    return(phendf)
}


# Filter the phenology dataframe by sex. Add grouping vars a la stan and extract unique combinations of grouping variables. df is a dataframe and sex is "MALE" or "FEMALE". Stan requires grouping vars to be consecutive integers, so indexes will differ/represent different underlying groups for any groups that are not identical across sexes
splitsex <- function(df = phendf, sex) {
    df <- filter(phendf, Sex == sex)
    # add grouping columns to match with stan output
    df <- stanindexer(df)
}


# join male and female sex
unique_grouper <- function(df1 = fdf, df2 = mdf) {
    df <- full_join(df1, df2)
    #extract unique groups
    udf <- df %>%
        dplyr::select(Sex, Site, SiteID, SPU_Name, ProvenanceID, Clone, CloneID, Year, YearID) %>%
        distinct()
    # Group by individual clone and group by individuals WITH sex
    udf$IndGroup <- group_indices(udf, SiteID, ProvenanceID, YearID, CloneID)
    udf$IndSexGroup <- group_indices(udf, Sex, IndGroup)
    return(udf)
}

# Build a dataframe where each row is a unique combination of parameters for each unique combination of sex, site, clone, year, and provenance that occurs in the data (IndSexGroup). #pardf is a dataframe with N rows of parameters for each sum_forcing, site, provenance, clone, and year combination that appear in the data, where N = number of draws from the posterior distribution
build_par_df <- function(mcmcdf, datdf = udf, sex) {

    udf <- dplyr::filter(datdf, Sex==sex)
    mcmcdf$iter <- 1:nrow(mcmcdf)

    singledimpars <- mcmcdf %>%
        dplyr::select(iter, beta, sigma_site, sigma_prov, sigma_clone, sigma_year, contains("kappa"), contains("mean")) %>%
        rename(kappa1 = `kappa[1]`) %>%
        rename(kappa2 = `kappa[2]`)

    siteb <- tidypar(mcmcdf, "b_site", "SiteID")
    provb <- tidypar(mcmcdf, "b_prov", "ProvenanceID")
    cloneb <- tidypar(mcmcdf, "b_clone", "CloneID")
    yearb <- tidypar(mcmcdf, "b_year", "YearID")

    clonemerge <- left_join(udf, cloneb) %>%
        select(-name)
    pardf <- left_join(udf, clonemerge)
    rm(clonemerge, cloneb)
    sitemerge <- left_join(udf, siteb) %>%
        select(-name)
    pardf <- left_join(pardf, sitemerge)
    rm(sitemerge, siteb)
    provmerge <- left_join(udf, provb) %>%
        select(-name)
    pardf <- left_join(pardf, provmerge)
    rm(provmerge, provb)
    yearmerge <- left_join(udf, yearb) %>%
        select(-name)
    pardf <- left_join(pardf, yearmerge)
    rm(yearmerge, yearb)
    pardf <- left_join(pardf, singledimpars)

    # pardf <- left_join(udf, clonemerge) %>%
    #     left_join(provmerge) %>%
    #     left_join(sitemerge) %>%
    #     left_join(yearmerge) %>%
    #     left_join(singledimpars)
}

#take a stan model dataframe and create a tidy dataframe for a given parameter. takes a dataframe, a string with the parameter name, and a string describing the param ID (e.g. par = "Site" and id="SiteID")
tidypar <- function(stanframe, param, id) {
    #wide to long format for parameters
    par <- stanframe %>%
        dplyr::select(contains(param), iter) %>%
        tidyr::gather(key = key, value = value, contains("b_")) %>%
        dplyr::mutate(id = str_extract(key, "[0-9]{1,}"))
    colnames(par) <- c("iter", "name", param, id)
    par[,ncol(par)] <- as.integer(par[,ncol(par)])
    return(par)
}


