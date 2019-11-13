# Predict from the model
# Simulate from the model using mean effects and sigma.

# Need parameter values (distributions?)
# Need accumulated forcing time series - use existing for now.

# No generated quantities block. I'm calculating fstart and fend, right? so I can just use the values for b_*_mean and sigma_* that I pulled 
# from the posterior when I simulated the model in the first place.


fmod_raw <- readRDS("2019-10-28phenologyFEMALE.rds") %>%
  as.data.frame() %>%
  select(-contains("state_rep")) 

mmod_raw <- readRDS("2019-10-28phenologyMALE.rds") %>%
  as.data.frame() %>%
  select(-contains("state_rep"))

fmod <- fmod_raw %>%
  select(contains("kappa"), beta, contains("mean"))
fmod$iter <- 1:nrow(fmod)
mmod <- mmod_raw %>%
  select(contains("kappa"), beta, contains("mean"))
mmod$iter <- 1:nrow(mmod)


calc_fstart_means <- 
calc_fstart <- function(df) {
  forcing <- with(df, {
    (logit(0.2) + kappa1)/(beta + b_site + b_prov + b_clone + b_year)
  })
  return(forcing)
}

fpredict <- fmod %>% 
  mutate(start = (logit(0.2) + `kappa[1]`)/(beta + b_site_mean + b_clone_mean + b_year_mean + b_prov_mean)) %>%
  mutate(end = (logit(0.8) + `kappa[2]`)/(beta + b_site_mean + b_clone_mean + b_year_mean + b_prov_mean)) %>%
  mutate(Sex="FEMALE") %>%
  select(iter, start, end, Sex)

mpredict <- mmod %>%
  mutate(start = (logit(0.2) + `kappa[1]`)/(beta + b_site_mean + b_clone_mean + b_year_mean + b_prov_mean)) %>%
  mutate(end = (logit(0.8) + `kappa[2]`)/(beta + b_site_mean + b_clone_mean + b_year_mean + b_prov_mean)) %>%
  mutate(Sex="MALE") %>%
  select(iter, start, end, Sex)

predictphen <- rbind(fpredict, mpredict) %>%
  pivot_longer(cols=c("start", "end"), names_to = "side", values_to = "accumulated_forcing" )

ggplot(predictphen, aes(x=accumulated_forcing, fill=Sex, linetype=side)) +
  geom_density(alpha=0.5) +
  scale_fill_viridis_d(end=0.8) +
  ggtitle("Forcing required to begin and end flowering") +
  theme_bw(base_size=18)

## Climate data ####
clim <- read.csv("data/all_clim_PCIC.csv", header=TRUE, stringsAsFactors=FALSE) %>%
  filter(forcing_type==forcingtype)
clim$siteyear <- paste(clim$Site, clim$Year, sep='')

mpredict

clim
 climset$DoY[findInterval(mpredict$start,
                                 climset$sum_forcing) + 1]
 