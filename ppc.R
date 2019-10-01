#### posterior predictive checks
## for forcing units across the range observed in my data (length X), simulate 30 phenological states for each model configuration (Y).
## This will result in X x 30 x Y states

source('phenology_functions.R')

forcingtype = 'scaled_ristos'



library(tidyverse)
library(rstan)
library(bayesplot)
library(rethinking)


#model
fmod <- readRDS("slopes_nc_FEMALE2019-09-16gq.rds") %>%
    as.data.frame()

mmod <- readRDS("slopes_nc_MALE2019-09-16gq.rds") %>%
    as.data.frame()

state_rep <- t(as.matrix(dplyr::select(fmod, contains("state_rep"))))

reps <- c(1:2400)
colnames(state_rep) <- paste("staterep", reps, sep="_")


# original data (this code should match relevant bits in run_stan)
phendf <- read_data()

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")

postcheckf <- cbind(fdf, state_rep)
postcheckf$recordID <- group_indices(postcheckf, DoY, Sex, Year, Site, Orchard, Clone, TreeUnique)

postcheckf <- gather(postcheckf, key = "yrep", value = "state", starts_with("staterep")) %>%
    mutate(repid = str_extract(yrep, "\\d+"))
postcheckf$staterep <- group_indices(postcheckf, state, yrep)
postcheckref <- filter(postcheckf, repid %in% sample(1:2400, 200)) # make it faster


fdf2 <- filter(fdf, Phenophase_Derived==2)
frep2 <- filter(postcheckf, state==2)

ggplot(frep2, aes(x=sum_forcing, group=yrep), alpha=0.01, color="gray") +
    stat_ecdf() +
    stat_ecdf(data=fdf2, aes(x=sum_forcing), color="red", inherit.aes=FALSE)

# interesting
ggplot(postcheckf, aes(x=DoY, group=staterep, color=as.factor(state))) +
    stat_ecdf(linetype=3) +
    stat_ecdf(aes(x=DoY, color=as.factor(Phenophase_Derived)), inherit.aes = FALSE) +
    facet_wrap("Year")

ggplot(fdf, aes(x=DoY, y=Phenophase_Derived)) +
    geom_point()

ggplot(data=postcheckf, aes(x=factor(state), y=DoY, group=yrep)) +
    geom_violin(position="identity", fill=NA, alpha=0.1) +
    geom_boxplot(data=fdf, aes(x=as.factor(Phenophase_Derived), y=DoY, fill=as.factor(Phenophase_Derived)), width=0.1, inherit.aes = FALSE) +
    scale_fill_viridis_d() +
    facet_grid(Year ~ Phenophase_Derived)

ggplot(data=frep2, aes(x=DoY, group=yrep)) +
    geom_density(alpha=0.1) +
    facet_wrap("Year")

ppc_bars(state, state_rep)
ppc_bars_grouped(state, state_rep, group = df$SiteID, facet_args = list(scales="free_y"))

ppc_bars_grouped(state, state_rep, group = df$ProvenanceID, facet_args = list(scales="free_y")) +
    ggtitle("by Provenance")





