library(tidyverse)
library(gtools)
library(assertthat)

apc_female_site <- calc_apc_of_variable(u="b_site", uid = "SiteID", v=c("b_year", "b_prov", "b_clone"), 
                                 vid = c("YearID", "ProvenanceID", "CloneID"), 
                                 phendata = fudf, model_params = fpardf, n=nrow(fmod), climlist=climlist)
