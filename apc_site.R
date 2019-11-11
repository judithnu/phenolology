library(tidyverse)
library(gtools)
library(assertthat)

site_apc <- calc_apc_of_variable(u="b_site", uid = "SiteID", v=c("b_year", "b_prov", "b_clone"), 
                                 vid = c("YearID", "ProvenanceID", "CloneID"), 
                                 phendata = udf, model_params = pardf)