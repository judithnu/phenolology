# Simulate phenology process

# Assume fixed probability of transition, bernoulli
transition <- function(state1 = 0, state2 = 1, covariate = heatsum) { #simulate the transition between two states and record the states.
    x <- c(state1) # initial state
    pt <- c() # initialize transition probability
    #initialize while loop
    transition <- 0 #transition y/n
    i <- 1
    while (ptrans < 1) {
        x <- append(x, state1) #record state1
        p <- 0.01 * covariate[i] #assume threshold heatsum is 100. this is likely to be a slow step and can be calculated outside the loop or function eventually
        pt <- append(pt, p) #record prob
        transition <- rbinom(1, size = 1, prob = p) # sample for new transition probability
        i <- i + 1
    }
    x <- append(x, state2) # record transition to state2 #this could be a real slow step

    return(list(states = x, probabilities = y))
}

b <- -3
c <- -1

ptry <- 1/(1 + exp(b * (temp - c)))
ptrycum <- cumsum(ptry)
#Scrit <- 242

seriesbuilder <- function(covariate) {
    firsttrans <- transition(covariate, state1 = 0, state2 = 1)
    secondtrans <- transition(covariate, state1 = 1, state2 = 2)
    series <- c(firsttrans, secondtrans)
    return(series)
}

#simulate temperature data
temp = rep(5, 50) #constant daily temperature of 5 degrees
heatsum <- cumsum(temp)

for (i in c(1:10)) {
    x <- seriesbuilder(heatsum)
    print(c("length", length(x[[1]])))
    print(x)
}

for (i in c(1:10)) {
    x <- seriesbuilder(heatsum)
    y <- c(1:10)
    print(c("length", length(x[[1]])))
    print(x)
}



