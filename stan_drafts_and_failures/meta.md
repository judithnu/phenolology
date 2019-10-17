# meta for stan_drafts_and_failures

# 2019-03_heatsum_intercepts

Stan model with 
 * heatsum 
 * intercepts- provenance, clone, site, and tree. 
 * slopes - sex

# 2019-03_ristos_intercepts

Stan models with
 * ristos
 * intercepts - provenance, clone, site, tree
 * slopes - sex

# 2019-04_correlation
Stan models with slopes and intercepts and attempts to account for the correlation between slopes and intercepts

# 2019-04_cutpoints
Stan model without hierarchy and priors on difference between cutpoints

# 2019_04_multilevel
Stan model with hierarchy, priors on cutpoints diff, 
 * intercepts - Sex, Prov, Clone, Site
 * slopes - Sex
Also includes a version with intercepts only

# 2019-04_separatesex
Stan models with hierarchy. Separate models for separate sexes. Priors on cutpoint diffs. Versions with effects on both slope and intercepts and with effects on only intercepts

#2019-08-25_sitexpcic_droppedKR
model run where weather correction was based on reconstructed site temperatures regressed on pcic pnwnamet data and with Kettle River's 2011 data dropped.

# baby ordered logistic

## generative_model.rds

## generative_model.stan
simulates model configurations from the prior and then simulates data from the generating process. Output to `simulated_data_from_gen_model.R`

## heatsum_intercepts
Rather similar to 2019-03_ristos_intercepts, but 

Stan model with
 * heatsum
 * intercepts - clone, provenance, site, tree
 * slopes - sex

## labmeetingfeb2019stanworkflow.R
version of `phenology.R` for lab meeting on March 4, 2019

## logit_and_ordered_logit_with_glm_and_map2stan.R
all with Prince George data. gradually exploring logit and ordered logit with glm and rethinking.

## phenology.R
workflow for phenology model

1. Design a generative model and use it to simulate data
2. Build a fitting model based on the generative model and fit it to the simulated data from the generative mode
3. Check the fit of the fitting model to the generative model
4. Fit the fitting model to real data
5. Evaluate the fit
6. Simulate data from the fitted model and compare to real data

## phenology.rds

## phenology_hier.stan

hierarchical ordered logistic model in stan
## simple ordered logistic with rethinking

An ordered logistic model with no hierarchy done with the rethinking package.

# forcing_to_doy.R

2019-04-27. attempting to model the relationship between forcing and day of year

## simulated_data_from_gen_model.R
output from `generative_model.stan` 

# simulations
Simulating data and fitting stan models to it.

# slopes_tree.stan

2019-07-11. added individual trees to stan and tightened priors on beta effects.

# slopes.stan

2019-07. Model for committee meeting 2019-07-30.

# workspaces

