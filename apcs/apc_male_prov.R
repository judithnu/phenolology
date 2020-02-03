library(tidyverse)
library(gtools)
library(assertthat)


apc_male_prov <- calc_apc_of_variable(u="b_prov", uid = "ProvenanceID", v=c("b_site", "b_year", "b_clone"), 
                                 vid = c("SiteID", "YearID", "CloneID"),
                                 phendata = mudf, model_params = mpardf, n=nrow(mmod), climlist=climlist)
