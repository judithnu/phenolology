library(tidyverse)
library(gtools)
library(assertthat)

apc_male_year <- calc_apc_of_variable(u="b_year", uid = "YearID", v=c("b_site", "b_prov", "b_clone"), 
                                 vid = c("SiteID", "ProvenanceID", "CloneID"), 
                                 phendata = mudf, model_params = mpardf, n=nrow(mmod), climlist=climlist)