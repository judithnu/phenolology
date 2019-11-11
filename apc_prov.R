library(tidyverse)
library(gtools)
library(assertthat)


prov_apc <- calc_apc_of_variable(u="b_prov", uid = "ProvenanceID", v=c("b_site", "b_year", "b_clone"), 
                                 vid = c("SiteID", "YearID", "CloneID"),
                                 phendata = udf, model_params = pardf)