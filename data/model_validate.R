# validate model
#read in summary of pollen data for Kalamalka

pollenmax <- read.csv("../phd/data/PhenologyAndPollenCounts/maxpollenslidesummary.csv", header=TRUE, stringsAsFactors = FALSE)
pollenmax$DoY <- lubridate::yday(pollenmax$Date)

male_hpd <- filter(doyperiod_hpd, Sex=="MALE") %>%
  select(-doylength) %>%
  pivot_wider(names_from = hpd_end, values_from = c(dbegin, dend))

pollenplus <- left_join(pollenmax, male_hpd) %>%
  mutate(pollenmax_in_interval = (DoY >=dbegin_hpd_low & DoY <=dend_hpd_high))

ggplot(pollenplus, aes(x=pollenmax_in_interval, fill=as.factor(intervalwidth))) +
  geom_bar(position="dodge") +
  ggtitle("Pollen slide max in HPDI?", subtitle = "3 years, 10 obs") +
  scale_fill_viridis_d(option="cividis", end=0.8)
