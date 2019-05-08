#calculate mean different in heatsum between days

diffdata <- ristoframe %>%
    select(Year, DoY, Site, sum_forcing) %>%
    group_by(Site, Year) %>%
    arrange(Site, Year,DoY) %>%
    mutate(diff = c(0, diff(sum_forcing))) %>%
    filter(DoY > 60 & DoY < 180)

hist(diffdata$diff)
sd(diffdata$diff)
hist(rnorm(1000,0,10))
hist(rexp(1000, .25))
