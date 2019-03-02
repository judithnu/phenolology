# replicating my model in brms

library(brms)

phendf <- read.csv("~/Documents/research_phenolology/data/stan_input/phenology_heatsum.csv", header=TRUE)

mdf <- subset(phendf, Sex == "MALE")

fit <- brm(
    formula = Phenophase_Derived ~ Heatsum + (1|Clone) + (Heatsum | Site),
       family = "cumulative",
       prior = c(set_prior("beta(.5,5)", class = "b", lb=0, ub=1)),
       data = mdf,
       chains = 5,
       cores = 5,
    save_model = "brms_hierachical.stan"
)
