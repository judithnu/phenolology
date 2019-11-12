library(tidyverse)
library(gtools)
library(assertthat)

apc_male_site <- calc_apc_of_variable(u="b_site", uid = "SiteID", v=c("b_year", "b_prov", "b_clone"), 
                                 vid = c("YearID", "ProvenanceID", "CloneID"), 
                                 phendata = mudf, model_params = mpardf, n=nrow(mmod), climlist=climlist)
