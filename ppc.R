#### posterior predictive checks
## for each phenological observation and model configuration, simulate 1 phenological state

source('phenology_functions.R')

forcingtype = 'scaled_ristos'


library(tidyverse)
library(rstan)
library(bayesplot)
#library(rethinking)


#model
fmod <- readRDS("2019-10-13_slopes_nc_FEMALE_.rds") %>%
    as.data.frame()

mmod <- readRDS("2019-10-13_slopes_nc_MALE_.rds") %>%
    as.data.frame()

state_repf <- t(as.matrix(dplyr::select(fmod, contains("state_rep"))))
state_repm <- t(as.matrix(dplyr::select(mmod, contains("state_rep"))))

reps <- c(1:nrow(fmod))
colnames(state_repm) <- paste("staterep", reps, sep="_")
colnames(state_repf) <- paste("staterep", reps, sep="_")

# original data (this code should match relevant bits in run_stan)
phendf <- read_data()

fdf <- splitsex(phendf, "FEMALE")
# fdf$recordID <- group_indices(fdf, Index, DoY) #create an index for each observation
# 
# dups <-  which(diff(fdf$recordID) < 1)
# fdf <- fdf %>% arrange(recordID)
# whydups <- fdf[sort(c(dups, dups+1)) ,]
mdf <- splitsex(phendf, "MALE")

#merge data and predicted phenological states for females
postcheckf <- cbind(fdf, state_repf)
postcheckm <- cbind(mdf, state_repm)
postrep <- rbind(postcheckf, postcheckm)
postrep$recordID <- group_indices(postrep, Index, DoY, Sex, Year, Site, Orchard, Clone, TreeUnique)

postrep <- gather(postrep, key = "yrep", value = "state", starts_with("staterep")) %>%
    mutate(repid = str_extract(yrep, "\\d+"))
postrep$state <- as.factor(postrep$state)
postrep$Phenophase_Derived <- as.factor(postrep$Phenophase_Derived)
postrep$staterep <- group_indices(postrep, Sex, state, yrep) 

#merge data and predicted phenological states for males

# postcheckm <- cbind(mdf, state_repm)
# postcheckm$recordID <- group_indices(postcheckm, DoY, Sex, Year, Site, Orchard, Clone, TreeUnique)
# 
# postcheckm <- gather(postcheckm, key = "yrep", value = "state", starts_with("staterep")) %>%
#   mutate(repid = str_extract(yrep, "\\d+"))
# postcheckm$state <- as.factor(postcheckm$state)
# postcheckm$Phenophase_Derived <- as.factor(postcheckm$Phenophase_Derived)
# postcheckm$staterep <- group_indices(postcheckm, state, yrep)


# calculate proportion of correct predictions by observation and phenophase
postrep <- postrep %>%
  mutate(correctyn = case_when(state==Phenophase_Derived ~ 1,
                             state!=Phenophase_Derived ~ 0)) 
propstate <- postrep %>%
  group_by(recordID, state, Sex) %>%
  summarise(propstate = n()/length(reps)) %>% # calculate proportion of each state predicted for each observation
  left_join(forcingonly) %>%
  #  %>%
    mutate(correctyn = case_when(state==Phenophase_Derived ~ 1,
                                 state!=Phenophase_Derived ~ 0))
ppc <- postrep %>%
  group_by(recordID) %>%
  summarise(prop_correct = sum(correctyn)/n()) %>% # calculate proportion correct for each observation
  #left_join(propstate) %>%
  left_join(postrep) %>%
  select(recordID, Index, DoY, Sex, Year, Site, Orchard, Clone, TreeUnique, sum_forcing,
         SPU_Name, Phenophase_Derived, prop_correct, state1, state2, state3) %>%
  distinct()

correctfindex <- which(postcheckf$state == postcheckf$Phenophase_Derived)
postcheckf$correct[correctfindex] <- 1
postcheckf$correct[is.na(postcheckf$correct)] <- 0

ppcf <- postcheckf %>%
    group_by(recordID) %>%
    mutate(n = length(recordID)) %>%
  mutate(prop_correct = sum(correct)/length(correct)) %>%
    select(recordID, DoY, Sex, Year, Site, Orchard, Clone, TreeUnique, sum_forcing,
           SPU_Name, Phenophase_Derived, prop_correct) %>%
    distinct()

toomany <- filter(ppcf, n==3600)

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
  group_by(Sex, Phenophase_Derived) %>%
  summarise(prop_correct = sum(correctyn)/n()) %>% # calculate proportion correct for each observation
  left_join(postrep) %>%
  select(recordID, DoY, Sex, Year, Site, Orchard, Clone, TreeUnique, sum_forcing,
         SPU_Name, Phenophase_Derived, prop_correct) %>%
  distinct() %>%
  arrange(Sex, Phenophase_Derived)

ppc_phenophase_site <- postrep %>%
  group_by(Sex, Phenophase_Derived, Site) %>%
  summarise(prop_correct = sum(correctyn)/n()) %>% # calculate proportion correct for each observation
  left_join(postrep) %>%
  select(recordID, DoY, Sex, Year, Site, Orchard, Clone, TreeUnique, sum_forcing,
         SPU_Name, Phenophase_Derived, prop_correct) %>%
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





