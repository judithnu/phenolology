library(tidyverse)
library(gtools)
library(assertthat)

year_apc <- calc_apc_of_variable(u="b_year", uid = "YearID", v=c("b_site", "b_prov", "b_clone"), 
                                 vid = c("SiteID", "ProvenanceID", "CloneID"), 
                                 phendata = udf, model_params = pardf, n=nrow(fmod), climlist=climlist)