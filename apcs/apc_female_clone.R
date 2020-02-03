library(tidyverse)
library(gtools)
library(assertthat)

apc_female_clone <- calc_apc_of_variable(u="b_clone", uid = "CloneID", v=c("b_site", "b_prov", "b_year"), 
                                  vid = c("SiteID", "ProvenanceID", "YearID"), 
                                  phendata = fudf, model_params = fpardf, n=nrow(fmod), climlist=climlist)