#### posterior predictive checks
## for forcing units across the range observed in my data (length X), simulate 30 phenological states for each model configuration (Y).
## This will result in X x 30 x Y states

source('phenology_functions.R')

forcingtype = 'scaled_ristos'


library(dplyr)
library(rstan)
library(bayesplot)
library(rethinking)


#model
fmod <- readRDS("slopes_nc_scaled_ristos_FEMALE2019-10-04climatena.rds") %>%
    as.data.frame()

mmod <- readRDS("slopes_nc_scaled_ristos_MALE2019-10-04climatena.rds") %>%
    as.data.frame()

state_repf <- t(as.matrix(dplyr::select(fmod, contains("state_rep"))))
state_repm <- t(as.matrix(dplyr::select(mmod, contains("state_rep"))))

reps <- c(1:1800)
colnames(state_repm) <- paste("staterep", reps, sep="_")
colnames(state_repf) <- paste("staterep", reps, sep="_")

# original data (this code should match relevant bits in run_stan)
phendf <- read_data()

fdf <- splitsex(phendf, "FEMALE")
mdf <- splitsex(phendf, "MALE")

#merge data and predicted phenological states for females
postcheckf <- cbind(fdf, state_repf)
postcheckf$recordID <- group_indices(postcheckf, DoY, Sex, Year, Site, Orchard, Clone, TreeUnique)

postcheckf <- gather(postcheckf, key = "yrep", value = "state", starts_with("staterep")) %>%
    mutate(repid = str_extract(yrep, "\\d+"))
postcheckf$state <- as.factor(postcheckf$state)
postcheckf$Phenophase_Derived <- as.factor(postcheckf$Phenophase_Derived)
postcheckf$staterep <- group_indices(postcheckf, state, yrep)
#postcheckref <- filter(postcheckf, repid %in% sample(1:2400, 200)) # make it faster

#merge data and predicted phenological states for males

postcheckm <- cbind(mdf, state_repm)
postcheckm$recordID <- group_indices(postcheckm, DoY, Sex, Year, Site, Orchard, Clone, TreeUnique)

postcheckm <- gather(postcheckm, key = "yrep", value = "state", starts_with("staterep")) %>%
  mutate(repid = str_extract(yrep, "\\d+"))
postcheckm$state <- as.factor(postcheckm$state)
postcheckm$Phenophase_Derived <- as.factor(postcheckm$Phenophase_Derived)
postcheckm$staterep <- group_indices(postcheckm, state, yrep)


# calculate proportion of correct predictions by observation and phenophase
ppcf <- postcheckf %>%
  mutate(correct = case_when(state==Phenophase_Derived ~ 1,
                             state!=Phenophase_Derived ~ 0)) %>%
  group_by(recordID) %>%
  mutate(prop_correct = sum(correct)/length(correct)) %>%
  select(recordID, DoY, Sex, Year, Site, Orchard, Clone, TreeUnique, sum_forcing,
         SPU_Name, Phenophase_Derived, prop_correct) %>%
  distinct()

ppcm <- postcheckm %>%
  mutate(correct = case_when(state==Phenophase_Derived ~ 1,
                             state!=Phenophase_Derived ~ 0)) %>%
  group_by(recordID) %>%
  mutate(prop_correct = sum(correct)/length(correct)) %>%
  select(recordID, DoY, Sex, Year, Site, Orchard, Clone, TreeUnique, sum_forcing,
         SPU_Name, Phenophase_Derived, prop_correct) %>%
  distinct()

ppc <- rbind(ppcf, ppcm)

#combine both sex's posterior prediction dfs
postrep <- rbind(postcheckf, postcheckm)

# calculate proportion of correct predictions by phenophase
ppc_phenophase <- postrep %>%
  mutate(correct = case_when(state==Phenophase_Derived ~ 1,
                             state!=Phenophase_Derived ~ 0)) %>%
  group_by(Sex, Phenophase_Derived) %>%
  mutate(prop_correct = sum(correct)/length(correct)) %>%
  select(Sex, Phenophase_Derived, prop_correct) %>%
  distinct() %>%
  arrange(Sex, Phenophase_Derived)

ppc_phenophase_site <- postrep %>%
  mutate(correct = case_when(state==Phenophase_Derived ~ 1,
                             state!=Phenophase_Derived ~ 0)) %>%
  group_by(Sex, Phenophase_Derived, Site) %>%
  mutate(prop_correct = sum(correct)/length(correct)) %>%
  select(Sex, Phenophase_Derived, Site, Site, prop_correct) %>%
  distinct() %>%
  arrange(Sex, Phenophase_Derived)


# plot model predictions
ggplot(ppc, aes(x=sum_forcing, y=prop_correct)) +
    geom_point(alpha=0.2, pch=16) +
    facet_grid(Sex ~ Phenophase_Derived) +
  ggtitle("Proportion of predictions that match data across forcing units")

ggplot(postrep, aes(x=Phenophase_Derived, fill=state)) +
  geom_bar() +
  ggtitle("Model predictions at each observed phenophase") +
  scale_fill_viridis_d(option="cividis") +
  facet_wrap(Sex ~ .)

ggplot(ppc_phenophase, aes(x=Sex, y=prop_correct, fill=Phenophase_Derived)) +
  geom_bar(stat="identity", position="dodge") +
  ylim(c(0,1)) +
  scale_fill_viridis_d(option="A") +
  ggtitle("Proportion of states correctly predicted by the model for each observed phenophase")

ggplot(ppc_phenophase_site, aes(x=Phenophase_Derived, y=prop_correct, fill=Site)) +
  geom_bar(stat="identity", position="dodge") +
  ylim(c(0,1)) +
  scale_fill_viridis_d(option="A") +
  facet_wrap("Sex") +
  ggtitle("Proportion of states correctly predicted by the model for each observed phenophase")

  #####################################
obs2 <- filter(postcheckf, Phenophase_Derived==2)
pred2 <- filter(postcheckf, state==2)

# When stage 2 is observed, the model predicts what stage?
ggplot(obs2, aes(x=sum_forcing, fill=state)) +
    geom_histogram(position = "identity", alpha=0.5) +
    ggtitle("Model predictions when observations are stage 2", subtitle = "female")

# When the model predicts stage 2, what is observed?
ggplot(pred2, aes(x=sum_forcing, fill=Phenophase_Derived)) +
    geom_histogram(position = "identity", alpha=0.5) +
    ggtitle("Observations when model predicts stage 2", subtitle = "female")



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

ppc_bars(state, state_repf)
ppc_bars_grouped(state, state_repf, group = df$SiteID, facet_args = list(scales="free_y"))

ppc_bars_grouped(state, state_repf, group = df$ProvenanceID, facet_args = list(scales="free_y")) +
    ggtitle("by Provenance")





