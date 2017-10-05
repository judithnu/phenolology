# Simulate phenology process


# Assume fixed probability of transition, bernoulli
transition <- function(p = 0.1, state1 = 0, state2 = 1) { #simulate the transition between two states and record the states.
    ptrans <- rbinom(1, size = 1, prob = p) #y/n transition from initial to next stage
    x <- c("pre")
    while (ptrans == 0) {
        x <- append(x, state1)
        ptrans <- rbinom(1, size = 1, prob = p)
    }
    x <- append(x, state2)
    return(x)
}

seriesbuilder <- function(p = 0.1, q = 0.1) {
    firsttrans <- transition(p, state1 = 0, state2 = 1)
    secondtrans <- transition(q, state1 = 1, state2 = 2)
    series <- c(firsttrans, secondtrans)
}

for (i in c(1:10)) {
    x <- seriesbuilder()
    print(x)
}

