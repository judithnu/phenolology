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

lengthframephen <- select(predictphen, iter, Sex, lengthforcing) %>% distinct()

# histograms of phenological period length
ggplot(lengthframephen, aes(x=lengthforcing, fill=Sex)) +
  geom_histogram() +
  scale_fill_viridis_d(end=0.8) +
  theme_bw(base_size=18) +
  ggtitle("Length of phenological period in forcing")

# Climate data ####
clim <- read.csv("data/all_clim_PCIC.csv", header=TRUE, stringsAsFactors=FALSE) %>%
  filter(forcing_type==forcingtype)
clim$siteyear <- paste(clim$Site, clim$Year, sep='')
ryears <- sample(unique(clim$siteyear), 20)
clim <- filter(clim, siteyear %in% ryears)

# CLIMATE PREDICTIONS ###########
splitclim <- split(clim, clim$siteyear)

ifinder <-  function(x) {
   index <- findInterval(predictphen$sum_forcing, x$sum_forcing)
   doy <- x$DoY[index]
   df <- data.frame(iter=predictphen$iter, DoY=doy, siteyear=x$siteyear[1],
                    side=predictphen$side, Sex=predictphen$Sex)
   return(df)
 }

predictdoy <- purrr::map_df(splitclim, ifinder) %>%
  pivot_wider(names_from = side, values_from = DoY) %>% #expensive operation
  mutate(lengthdays = end - begin) %>%
  pivot_longer(cols=c(end, begin), names_to = "side", values_to = "DoY") %>%
  left_join(predictphen)

assertthat::assert_that(nrow(predictdoy)==nrow(predictphen) * length(splitclim))
assertthat::assert_that(length(unique(predictdoy$siteyear))==length(splitclim)) # did you select the right/enough climate time series?

lengthframe <- select(predictdoy, iter, Sex, siteyear, contains("length")) %>%
  distinct()

# violin plot of start and end day of year for every siteyear
ggplot(predictdoy, aes(x=siteyear, y=DoY, color=Sex)) +
  geom_violin() +
  facet_grid(. ~ side)+
  theme(legend.position = "none") +
  scale_color_viridis_d(end=0.8) +
  scale_x_discrete(labels=NULL)

ggplot(lengthframe, aes(x=lengthdays, fill=Sex)) +
  geom_histogram(binwidth=1) +
  facet_grid(Sex ~ .) +
  scale_fill_viridis_d(end=0.8) +
  ggtitle("histograms of phenological period length")

ggplot(lengthframe, aes(x=siteyear, y=lengthdays, color=Sex)) +
  geom_violin() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_viridis_d(end=0.8)

shortid <- c("Tolko2006", "Vernon2001")
longid <- c("Kalamalka2010", "Sorrento2010")

climphen <- filter(clim, DoY < 180)
short <- filter(climphen, siteyear %in% shortid)
long <- filter(climphen, siteyear %in% longid)

forcingonly <- select(predictdoy, sum_forcing, Sex, side)
ggplot(climphen, aes(x=DoY, y=sum_forcing, group=siteyear)) +
  geom_line(color="gray") +
  geom_line(data=short, aes(color="short")) +
  geom_line(data=long, aes(color="long")) +
  scale_color_viridis_d(option="B", end=0.8) +
  geom_hline(data=predictdoy, aes(yintercept = sum_forcing, linetype=Sex), alpha=0.2)


meantemps <- left_join(predictdoy, select(clim, siteyear, DoY, mean_temp), by=c("siteyear", "DoY"))

hpd_lower = function(x, prob) rethinking::HPDI(x, prob)[1]
hpd_upper = function(x, prob) rethinking::HPDI(x, prob)[2]

predictsummary <- meantemps %>% group_by(siteyear, Sex, side) %>%
  summarize(lowerDoY=hpd_lower(DoY, 0.5), upperDoY=hpd_upper(DoY, 0.5), medDoY=median(DoY),
            lowerforcing=hpd_lower(sum_forcing, 0.5), upperforcing=hpd_upper(sum_forcing, 0.5), medforcing=median(sum_forcing),
            lowerlengthforcing=hpd_lower(lengthforcing, 0.5), upperlengthforcing=hpd_upper(lengthforcing, 0.5), medlengthforcing = median(lengthforcing),
            lowerlengthDoY=hpd_lower(lengthdays, 0.5), upperlengthDoY=hpd_upper(lengthdays, 0.5), medlengthDoY=median(lengthdays))

ggplot(predictsummary, aes(x=side, y=medDoY)) +
  geom_point() +
  #geom_errorbar(aes(x=side, ymin=lowerDoY, ymax=upperDoY)) +
  geom_line(aes(group=siteyear)) +
  facet_wrap("Sex")

sumwithclim <- clim %>%
  filter(DoY > 105 & DoY < 170)


ggplot(sumwithclim, aes(x=sum_forcing, y=DoY, group=siteyear, color=siteyear)) +
         geom_line() +
  geom_vline(data=predictsummary, aes(xintercept = medforcing)) +
  facet_grid(Sex ~ .) +
  ggtitle("median day of year for flowering start and end", subtitle="with forcing accumulation") +
  scale_color_viridis_d()+
  theme_bw(base_size=18) +
  xlim(c(8,15))

ggplot(predictdoy, aes(x=lengthdays, y=siteyear, color=Sex) )+
  geom_point(pch=1)

# no
ggplot(meantemps, aes(x=DoY, y=sum_forcing, group=as.factor(iter), color=siteyear)) +
  geom_line(alpha=0.3) +
  facet_wrap("Sex")

ggplot(meantemps, aes(x=DoY, y=mean_temp, group=as.factor(iter)), color=siteyear) +
  geom_line(alpha=0.3) +
  geom_point(aes(x=DoY, y)



# ggplot(filter(predictdoy, siteyear %in% ryears), aes(x=DoY, fill=Sex)) +
#   geom_density(alpha=0.5) +
#   facet_grid(. ~ side)+
#   theme(legend.position = "none") +
#   scale_fill_viridis_d(end=0.8) +
#   scale_x_discrete(labels=NULL)

phendf <- read_data(slim = FALSE)

firstandlastrf <- phendf %>%
  group_by(Index, Sex, Phenophase_Derived) %>%
  filter(Phenophase_Derived!=1) %>%
  summarise(DoY=min(DoY), sum_forcing=min(sum_forcing)) #%>%
  #left_join(select(phendf, Index, Sex, DoY, sum_forcing))
firstandlastrf$side <- NA
firstandlastrf$side[which(firstandlastrf$Phenophase_Derived==2)] <- "begin"
firstandlastrf$side[which(firstandlastrf$Phenophase_Derived==3)] <- "end"

ggplot(sample_frac(predictdoy,0.05), aes(x=sum_forcing, y=DoY)) +
  #stat_binhex(binwidth=c(0.2, 1), alpha=0.7) +
  geom_density_2d(color="gray") +
  #geom_jitter(data=filter(phendf, Phenophase_Derived==2), aes(x=sum_forcing, y=DoY, color=Sex), alpha=0.1, pch=1) +
  #geom_jitter(data=firstandlastrf, aes(x=sum_forcing, y=DoY, color=beginorend), alpha=0.2, pch=1) +
  #geom_density_2d(data=firstandlastrf, aes(x=sum_forcing, y=DoY, color=beginorend))+
  facet_grid(Sex ~ .) +
  scale_color_viridis_d(option="A", end=0.8) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  ylab("Day of Year") +
  xlab("accumulated forcing") +
  ggtitle("Flowering start and end date", subtitle="forcing accumulation and day of year")

# Calculate length of phenological period


# calculate hdpi for each site year, or density?
hpd_lower = function(x, prob) rethinking::HPDI(x, prob)[1]
hpd_upper = function(x, prob) rethinking::HPDI(x, prob)[2]

fifty <- group_by(predictdoy, siteyear, Sex, side) %>%
  summarize(x=hpd_lower(DoY, prob=.99), xend=hpd_upper(DoY, prob=.99), median=median(DoY)) %>%
  mutate(interval="0.5")

ggplot(data=fifty, aes(x=siteyear, y=median)) +
  geom_line(color="gray") +
  geom_point(aes(color=side)) +
  geom_errorbar(aes(ymin=x, ymax=xend)) +
  facet_wrap("Sex") +
  ylab("Day of Year") +
  scale_color_viridis_d(option="B", end=0.9) +
  ggtitle("Flowering periods in all accumulated forcing timeseries")

ggplot(data=predictdoy, aes(x=side, y=DoY, fill=Sex)) +
  geom_violin(alpha=0.8) +
  scale_fill_viridis_d(end=0.8) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  ylab("Day of Year") +
  xlab("") +
  ggtitle("Day of Year flowering starts and ends") +
  geom_violin(data=firstandlastrf, aes(x=side, y=DoY, fill=Sex), alpha=0.5)

ggplot(data=predictdoy, aes(x=side, y=DoY, fill=Sex)) +
  geom_violin(alpha=0.8) +
  scale_fill_viridis_d(end=0.8) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  ylab("Day of Year") +
  xlab("") +
  ggtitle("Day of Year flowering starts and ends") +
  geom_violin(data=firstandlastrf, aes(x=side, y=DoY, fill=Sex), alpha=0.5)

ggplot(data=predictphen, aes(x=side, y=sum_forcing, fill=Sex)) +
  geom_violin(alpha=0.8) +
  scale_fill_viridis_d(end=0.8) +
  theme_bw(base_size=18) +
  theme(legend.position = "none") +
  ylab("Accumulated forcing") +
  xlab("") +
  ggtitle("Forcing required for flowering start and end") +
  geom_violin(data=firstandlastrf, aes(x=side, y=sum_forcing, fill=Sex), alpha=0.5) #+
  #facet_wrap("Phenophase_Derived")

# better time series #################

medforc = 12.3
meddoy = 150
a = c(0.01, 0.1, 1, 10)
b <- ((medforc)/a)^(1/meddoy)

ab <- data.frame(a=a, b=b)

x <- seq(0, meddoy, by=1)

forcingfunction <- function(x, a, b) {
  y <- a*(b^x)
  return(y)
}

forcingfunction(150, a, b)
forcingfunction(100, a, b)
forcingfunction(50, a, b)

ts <- data.frame(a=NA, b=NA, paramset=NA, DoY=NA, sum_forcing=NA)
for (i in 1:nrow(ab)) {
  a <- ab$a[i]
  b <- ab$b[i]
  paramset=i
  DoY=x
  sum_forcing <- forcingfunction(x, a, b)
  tst <- data.frame(a, b, paramset, DoY, sum_forcing)
  ts <- rbind(ts, tst)
}
ts <- ts[-1,]



ggplot(filter(clim, DoY < 152), aes(x=DoY, y=mean_temp, group=siteyear)) +
  geom_line()

climit <- filter(clim, DoY < 152)
fit <- lm(log(climit$sum_forcing) ~ climit$DoY)
climit$sum_forcing_pred <- exp(predict(fit,list(DoY=climit$DoY)))
climit$sum_forcing_t <- exp((climit$DoY*0.04 - 3.48))

test <- exp((climit$DoY*0.029)+-1.48)

b <- seq(from=0, to=5, by =1)
a <- (log(12.3) - b)/150
ab <- data.frame(a=a, b=b)

forcingfunction2 <- function(x, a, b) {
  y <- exp((x^a) + b)
  return(y)
}2b7q8


forcingfunction2(150, a, b)

ts2 <- data.frame(a=NA, b=NA, paramset=NA, DoY=NA, sum_forcing=NA)
for (i in 1:nrow(ab)) {
  a <- ab$a[i]
  b <- ab$b[i]
  paramset <- i
  DoY <- x
  sum_forcing <- forcingfunction2(x, a, b)
  tst <- data.frame(a, b, paramset, DoY, sum_forcing)
  ts2 <- rbind(ts2, tst)
}
ts2 <- ts2[-1,]

ggplot(ts2, aes(x=DoY, y=sum_forcing, color=as.factor(paramset))) +
  geom_line() #+
  geom_line(data=filter(clim, DoY < 152), aes( x=DoY, y=sum_forcing, group=siteyear), color="gray", inherit.aes = FALSE, alpha=0.5) +
  geom_point(data=filter(clim, DoY < 152), aes(x=DoY, y=sum_forcing_pred), inherit.aes=FALSE, color="blue")

ggplot(climit, aes(x=DoY, y=sum_forcing, group=siteyear)) +
  geom_line(color="gray") +
  geom_line(data=climit, aes(x=DoY, y=sum_forcing_pred), color="blue") +
  geom_line(data=climit, aes(x=DoY, y=sum_forcing_t), color="red")

