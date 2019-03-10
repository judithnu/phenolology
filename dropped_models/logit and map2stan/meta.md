#logit and map2stan

##maptry.R
Simulates fake phenology data with one transition only then fits it with logit, glmer, and map2stan (with individual effects on h and k)

## logit_and_ordered_logit_with_glm_and_map2stan.R
Basically maptry, but with real data. Also added an ordered logit model and started using prepollination priors. Couldn't add individual effects though in map2stan.

## processsimulationR
Generate phenofakes object (simulated phenophases) used in `maptry.R`

## simulateddata.R
Early attempt at simulating data from a logit
