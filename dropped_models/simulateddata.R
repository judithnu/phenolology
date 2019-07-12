# Simulate data
## crap
a <- 0
b <- 1
sigma <- 2
temp = 5
threshold = 3
dd <- max(0, temp - threshold)
cum_temp <- cumsum(temp)

log(.2/(1-.2))
log(.8/(1-0.8))

rnorm(100, mean = a + b * cum_temp, sd = sigma)
##

## i can make a logit
fauxdat <- data.frame(ind = rep(1:10,2), atemp = c(rnorm(10, 3, 5), rnorm(10, 15, 5)), fl = sort(rep(c(.2,.8),10)))

logit <- glm(fl ~ atemp, family = gaussian(link = "logit"), data = fauxdat)

summary(logit)

plot(fauxdat$atemp, fauxdat$fl)
curve(predict(logit, data.frame(atemp = x), type = "resp"), add = TRUE)

##

## tree level toy
### Tree has 10 cones. Over a 10 day period, cones become active for some amount of time, then finish. Cones become active at heatsum = 20 with some error and require 10 additional heatsum to finish.
library(dplyr)
library(tidyr)
library(ggplot2)

get_closest_bigger_number <- function(value, vector) { #gets index of number in a vector that's closest to the value while also being larger than the value
    diffs <- value - vector
    position <- which.min(diffs > 0)
    return(position)
}

calc_event_day <- function(crittemp, acumtemp = tdat$atemp, dayvec = tdat$day) { #determine day that a critical temperature is reached. acumtemp is the accumulated temperature on a day of dayvec.
    pos <- sapply(crittemp, get_closest_bigger_number, vector = acumtemp)
    crit_day <- dayvec[pos] #day tree "flowers"
}

temp = 5
temp_crit = 20
ncones = 100

tdat <- data.frame(day = c(1:10), temp = 5) #daily temp data
tdat$atemp <- cumsum(tdat$temp)

##critical temperatures for individual cones to begin and finish flowering
begin_crit <- rnorm(n = ncones, mean = temp_crit, sd = temp_crit/4)
end_crit <- rnorm(n = ncones, mean = temp_crit/2, sd = temp_crit/8) + begin_crit

cdat <- data.frame(cone = c(1:ncones), begin_crit, end_crit) #cone data

#count up events daily
cdat$begin_day <- calc_event_day(crittemp = cdat$begin_crit)
cdat$end_day <- calc_event_day(crittemp = cdat$end_crit)

#count events by day daily
dailyevents <- cdat %>%
    select(begin_day, end_day) %>%
    gather(key = "event", value = "day") %>%
    group_by(event, day) %>%
    summarize(daily_events = n())

#add in 0 days
dayframe <- data.frame(day = sort(rep(1:10, 2)), event = c("begin_day", "end_day"), daily_events = 0)

fulldaily <- dailyevents %>%
    anti_join(x = dayframe, by = c("day", "event")) %>%
    full_join(dailyevents)

# count cumulative events
cumulevents <- fulldaily %>%
    group_by(event) %>%
    arrange(event, day) %>%
    mutate(cumul_events = cumsum(daily_events))

eventdays <- cumulevents %>%
    filter(daily_events > 0)

phenperiod <- summarise(eventdays, first = min(day), last = max(day)) %>%
    gather(first, last, key = "occurence", value = "day")#get ends of phenological periods. begin_day rows have the first and last day that any cone began flowering on the tree. end_day rows have the first and last day any cone stopped flowering on the tree

sparsephenperiod <- phenperiod %>%
    summarise(start = min(day), end = max(day)) %>%
    gather(start, end, key = "tree_event", value = "day")

ggplot(cumulevents, aes(x = day, y = daily_events, fill = event)) +
    geom_col(alpha = .5, position = "dodge") +
    xlim(1,10) +
    geom_line(aes(x = day, y = cumul_events, color = event)) +
    ylab("cone count") +
    ggtitle("number of cones beginning and finishing flowering daily \n with daily counts and cumulative counts")

#make phenperiod plots
ggplot(phenperiod, aes(x = day, y = 0, color = event) ) +
    geom_point(size = 5)

ggplot(filter(cumulevents, daily_events >0), aes(x = event, y = daily_events)) +
    geom_violin() +
    geom_point()
