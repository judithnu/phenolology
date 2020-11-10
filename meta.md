# research_phenolology
 
## apcs
Average Predictive Comparisons
## data

Data used in models

## dropped_models

abandoned modeling approaches

## paper

Paper for phenology project

## priors

Use [Betancourt's phenology case study](https://betanalpha.github.io/assets/case_studies/ordinal_regression.html) to understand and test a Dirichlet prior

## stan_drafts_and_failures
Learning stan and refining the stan model

## FUs_to_DoY.R

Transform forcing unit predictions from model into day of year predictions. Requires tpars from ppc.R

## model_diagnostics.R

Look at stan model diagnostics and estimates with shinystan and bayesplot

## modeling_clarity*

Simple versions of models and data viz with different predictors to determine how to structure model in stan.

## model_predict.R

Predict fstart and fend in general from the model. Simulate effects based on model parameters. Then simulate fend and fstart based on those effects along with the other parameter values for each iteration.

## model_results.R

Calculate fstart and fend for each iteration and site/prov/year/clone combination in the data.

## phenology.stan


## ppc
NEEDS RENAMING
calculate transformed parameters - calculate params in terms of forcing units.

## ppc.Rmd
Was going to be R markdown version of ppc, but abandoned.

## priorinvestigation.R
What do the priors look like

## README.md

## run_stan.R
R script to process data, run, and save stan model.

## slopes_ristos_scaled*.rds
Stan model output

## transformed_parameters.csv
Stan model output parameters calculated in terms of forcing units. Downsampled.

















