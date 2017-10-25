# Simulate phenology process

# Assume fixed probability of transition, bernoulli
transition <- function(state1 = 0, state2 = 1, covariate = heatsum) { #simulate the transition between two states and record the states.
    x <- c(state1) # initial state
    transition <- 0 # initialize transition y/n
    b <- -0.1 # heatsum function parameter
    c <- 50 # heatsum function parameter
    #initialize while loop
    pt <- c() #intialize pt vec
    i <- 1
    while (transition < 1) {
        x <- append(x, state1) #record state1
        p <- 1/(1 + exp(b * (covariate[i] - c))) #function from Chuine's universal
        pt <- append(pt, p) #record prob
        transition <- rbinom(1, size = 1, prob = p) # sample for new transition probability
        i <- i + 1
    }
    x <- append(x, state2) # record transition to state2 #this could be a real slow step
    pt <- c(NA, pt, 1)

    return(data.frame(states = x, probabilities = pt, heatsum = covariate[1:length(x)]))
}

#simulate temperature data
temp = rep(5, 50) #constant daily temperature of 5 degrees
heatsum <- cumsum(temp)



seriesbuilder <- function(covariate) {
    firsttrans <- transition(covariate, state1 = 0, state2 = 1)
    secondtrans <- transition(covariate, state1 = 1, state2 = 2)
    series <- rbind(firsttrans, secondtrans)
    return(series)
}



phenofakes <- data.frame(ind = c(), states = c(), probabilities = c(), heatsum = c())
for (i in c(1:100)) {
    x <- seriesbuilder(heatsum)
    #print(c("length", length(x[[1]])))
    indset <- cbind(ind = i, x)
    phenofakes <- rbind(phenofakes, indset)
}

library(ggplot2)

ggplot(phenofakes, aes(x = probabilities, y = states, color = as.factor(ind))) +
    geom_point() +
    geom_jitter()

ggplot(phenofakes[which(phenofakes$states < 2),], aes(x = heatsum, y = probabilities)) +
    geom_jitter() +
    stat_function(fun=function(x) 1/(1 + exp(-0.1 * (x - 50))))

