# Predict from the model
# Simulate from the model using mean effects and sigma.

# Need parameter values (distributions?)
# Need accumulated forcing time series - use existing for now.

# No generated quantities block. I'm calculating fstart and fend, right? so I can just use the values for b_*_mean and sigma_* that I pulled
# from the posterior when I simulated the model in the first place.

library(tidyverse)
library(gtools)

source('phenology_functions.R')


fmod_raw <- readRDS("2019-10-28phenologyFEMALE.rds") %>%
  as.data.frame() %>%
  select(-contains("state_rep"))

mmod_raw <- readRDS("2019-10-28phenologyMALE.rds") %>%
  as.data.frame() %>%
  select(-contains("state_rep"))

# GLOBALS ###########
forcingtype="scaled_ristos"

# MODELS #############
fmod <- fmod_raw %>%
  select(contains("kappa"), beta, contains("mean"))
fmod$iter <- 1:nrow(fmod)
rm(fmod_raw)
mmod <- mmod_raw %>%
  select(contains("kappa"), beta, contains("mean"))
mmod$iter <- 1:nrow(mmod)
rm(mmod_raw)

# FORCING PREDICTIONS #############
fpredict <- fmod %>%
  mutate(begin = (logit(0.2) + `kappa[1]`)/(beta + b_site_mean + b_clone_mean + b_year_mean + b_prov_mean)) %>%
  mutate(end = (logit(0.8) + `kappa[2]`)/(beta + b_site_mean + b_clone_mean + b_year_mean + b_prov_mean)) %>%
  mutate(Sex="FEMALE") %>%
  select(iter, begin, end, Sex)

mpredict <- mmod %>%
  mutate(begin = (logit(0.2) + `kappa[1]`)/(beta + b_site_mean + b_clone_mean + b_year_mean + b_prov_mean)) %>%
  mutate(end = (logit(0.8) + `kappa[2]`)/(beta + b_site_mean + b_clone_mean + b_year_mean + b_prov_mean)) %>%
  mutate(Sex="MALE") %>%
  select(iter, begin, end, Sex)

rm(fmod, mmod)

predictphen <- rbind(fpredict, mpredict) %>%
  mutate(lengthforcing = end-begin) %>%
  pivot_longer(cols=c("end", "begin"), names_to = "side", values_to = "sum_forcing" )

rm(mpredict, fpredict)
gc()

#scaled density of start and end with accumulated forcing (or flowering period - see commented line)
ggplot(predictphen, aes(x=sum_forcing, y= ..scaled.., fill=Sex, linetype=side)) +
  geom_density(alpha=0.8) +
  #geom_density(data=filter(phendf, Phenophase_Derived==2), aes(x=sum_forcing, y= ..scaled..), inherit.aes = FALSE) +
  stat_ecdf(data=filter(phendf, Phenophase_Derived==2), aes(x=sum_forcing), inherit.aes = FALSE) +
  scale_fill_viridis_d(end=0.8) +
  scale_color_viridis_d(end=0.8) +
  ggtitle("Forcing required to begin and end flowering", subtitle = "with cumulative distribution of accumulated forcing \non observed flowering dates") +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  xlab("Forcing accumulation") +
  ylab("") +
  facet_grid(Sex ~ .)



# Climate data ####
clim <- read.csv("data/all_clim_PCIC.csv", header=TRUE, stringsAsFactors=FALSE) %>%
  filter(forcing_type==forcingtype)
clim$siteyear <- paste(clim$Site, clim$Year, sep='')

# CLIMATE PREDICTIONS ###########
splitclim <- split(clim, clim$siteyear)



ifinder <-  function(x) {
   index <- findInterval(predictphen$sum_forcing, x$sum_forcing)
   doy <- x$DoY[index]
   df <- data.frame(iter=predictphen$iter, DoY=doy, siteyear=x$siteyear[1], Sex=predictphen$Sex, side=predictphen$side, sum_forcing=predictphen$sum_forcing)
   return(df)
 }

DoYpredict <- purrr::map_df(splitclim, ifinder)
assertthat::assert_that(nrow(DoYpredict)==nrow(predictphen) * length(splitclim))

ryears <- sample(unique(clim$siteyear), 30)
ggplot(filter(DoYpredict, siteyear %in% ryears), aes(x=siteyear, y=DoY, color=Sex)) +
  geom_violin() +
  facet_grid(. ~ side)+
  theme(legend.position = "none") +
  scale_color_viridis_d(end=0.8) +
  scale_x_discrete(labels=NULL)

# ggplot(filter(DoYpredict, siteyear %in% ryears), aes(x=DoY, fill=Sex)) +
#   geom_density(alpha=0.5) +
#   facet_grid(. ~ side)+
#   theme(legend.position = "none") +
#   scale_fill_viridis_d(end=0.8) +
#   scale_x_discrete(labels=NULL)

phendf <- read_data(slim = FALSE)

firstandlastrf <- phendf %>%
  group_by(Index, Sex) %>%
  filter(Phenophase_Derived==2) %>%
  summarise(firstrf=min(DoY), lastrf=max(DoY)) %>%
  pivot_longer(cols=contains("rf"), names_to = "beginorend", values_to = "DoY") %>%
  left_join(select(phendf, Index, Sex, DoY, sum_forcing))

ggplot(sample_frac(DoYpredict,0.2), aes(x=sum_forcing, y=DoY)) +
  #stat_binhex(binwidth=c(0.2, 1), alpha=0.7) +
  geom_density_2d(color="gray") +
  #geom_jitter(data=filter(phendf, Phenophase_Derived==2), aes(x=sum_forcing, y=DoY, color=Sex), alpha=0.1, pch=1) +
  #geom_jitter(data=firstandlastrf, aes(x=sum_forcing, y=DoY, color=beginorend), alpha=0.2, pch=1) +
  geom_density_2d(data=firstandlastrf, aes(x=sum_forcing, y=DoY, color=beginorend))+
  facet_grid(Sex ~ .) +
  scale_color_viridis_d(option="A", end=0.8) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  ylab("Day of Year") +
  xlab("accumulated forcing") +
  ggtitle("Flowering start and end date", subtitle="forcing accumulation and day of year")


# plot stage 1 and 3 in above graph, too?
