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

