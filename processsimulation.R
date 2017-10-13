# Simulate phenology process

# Assume fixed probability of transition, bernoulli
transition <- function(state1 = 0, state2 = 1, covariate = heatsum) { #simulate the transition between two states and record the states.
    ptrans <- 0 #initial probability of transition
    x <- c(state1)
    y <- c(ptrans)
    i <- 1
    while (ptrans < 1) {
        x <- append(x, state1) #record state1
        p <- 0.01 * covariate[i] #assume threshold heatsum is 100
        y <- append(y, p) #record prob
        ptrans <- rbinom(1, size = 1, prob = p) # sample for new transition probability
        i <- i + 1
    }
    x <- append(x, state2) # record transition to state2

    return(list(states = x, probabilities = y))
}

seriesbuilder <- function(covariate) {
    firsttrans <- transition(covariate, state1 = 0, state2 = 1)
    secondtrans <- transition(covariate, state1 = 1, state2 = 2)
    series <- c(firsttrans, secondtrans)
}

#simulate temperature data
temp = rep(5, 50) #constant daily temperature of 5 degrees
heatsum <- cumsum(temp)

for (i in c(1:10)) {
    x <- seriesbuilder(heatsum)
    print(c("length", length(x[[1]])))
    print(x)
}


