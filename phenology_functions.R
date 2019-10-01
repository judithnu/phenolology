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

read_data <- function() {
    phenology_data <- read.csv("data/phenology_heatsum.csv",
                               stringsAsFactors = FALSE, header = TRUE
    ) %>%
        filter(!(Year==2011 & Site=="KettleRiver"))

    ## provenance
    SPU_dat <- read.csv("../research_phd/data/OrchardInfo/LodgepoleSPUs.csv",
                        header=TRUE, stringsAsFactors = FALSE) %>%
        dplyr::select(SPU_Name, Orchard)

    # join provenance and phenology data

    phendf <- phenology_data %>%
        na.omit() %>%
        left_join(SPU_dat) %>%
        unique()

    return(phendf)
}
