# research_phenolology
 
## simple ordered logistic with rethinking

An ordered logistic model with no hierarchy done with the rethinking package.

## data

Data used in models

## hmm

Hidden Markov Model attempt

## paper

Paper for phenology project

## workspaces

Workspaces saved for slow models for later eval

## graphsforsally

Graphs for meetings with supervisor

## ClarkCDMcode
Partial code for models in Clark's 2014 papers

1. Clark JS, Melillo J, Mohan J, Salk C. The seasonal timing of warming that controls onset of the growing season. Global Change Biology. 2014 Oct 1;20(4):1136–45. 
2. Clark JS, Salk C, Melillo J, Mohan J. Tree phenology responses to winter chilling, spring warming, at north and south range limits. Anten N, editor. Functional Ecology. 2014 Jul 22;28(6):1344–55. 

## brms_hierarchical.R
March 2019 attempt to create a hierarchical ordinal logistic model using the brms package - extremely slow

## brms_hierarchical.stan
Stan model version of `brms_hierarchical.R`

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

## simulated_data_from_gen_model.R
output from `generative_model.stan` 

## generative_model.rds

## generative_model.stan
simulates model configurations from the prior and then simulates data from the generating process. Output to `simulated_data_from_gen_model.R`

## labmeetingfeb2019stanworkflow.R
version of `phenology.R` for lab meeting on March 4, 2019

## logit_and_ordered_logit_with_glm_and_map2stan.R
all with Prince George data. gradually exploring logit and ordered logit with glm and rethinking.



