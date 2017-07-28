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
### Tree has 10 cones. Over a 10 day period, cones become active for some amount of time, then finish. Cones become active at DD = 20 with some error.

temp = 5
temp_crit = 20

tdat <- data.frame(day = c(1:10), temp = 5)
tdat$atemp <- cumsum(tdat$temp)

temp_crit_ind <- rnorm(n = 10, mean = temp_crit, sd = temp_crit/4) #critical temperatures for individual cones

cdat <- data.frame(cone = c(1:10), temp_crit_ind = temp_crit_ind)

get_closest_bigger_number <- function(value, vector) { #gets index of number in a vector that's closest to the value while also being larger than the value
    diffs <- value - vector
    position <- which.min(diffs > 0)
    return(position)
}

day_crit_pos <- sapply(cdat$temp_crit_ind, get_closest_bigger_number, vector = tdat$atemp)
cdat$day_crit <- tdat$day[day_crit_pos] #day tree "flowers"

library(dplyr)

begin_dat <- left_join(tdat, cdat, by = c("day" ="day_crit"))
table(begin_dat$day, begin_dat$cone)



library(ggplot2)
ggplot(begin_dat, aes(x = day, group=cone)) +
    geom_bar(na.rm = TRUE)
