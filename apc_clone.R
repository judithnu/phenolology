library(tidyverse)
library(gtools)
library(assertthat)

clone_apc <- calc_apc_of_variable(u="b_clone", uid = "CloneID", v=c("b_site", "b_prov", "b_year"), 
                                  vid = c("SiteID", "ProvenanceID", "YearID"), 
                                  phendata = udf, model_params = pardf, n=nrow(fmod), climlist=climlist)