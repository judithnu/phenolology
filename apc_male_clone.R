library(tidyverse)
library(gtools)
library(assertthat)

apc_male_clone <- calc_apc_of_variable(u="b_clone", uid = "CloneID", v=c("b_site", "b_prov", "b_year"), 
                                  vid = c("SiteID", "ProvenanceID", "YearID"), 
                                  phendata = mudf, model_params = mpardf, n=nrow(mmod), climlist=climlist)