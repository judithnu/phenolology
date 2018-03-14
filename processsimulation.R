# Simulate phenology process

# Phenology transition function
## Logistic function with asymptote at 1, steepness of curve determined by k and midpoint determined by h.
calc_probability <- function(x, k = steepness, h = midpoint, i = individual_effect) { # a logistic function parameterized to estimate steepness and transition points with individual effects
    1/(1 + exp(-(k + i) * (x - h)))
}

# Assume fixed probability of transition, bernoulli
transition <- function(state1 = 0, state2 = 1, p) { #simulate the transition between two states and record the states. p is the probability of transition. Can be a vector
    x <- c(state1) # initial state
    do_trans <- 0 # initialize transition y/n
    #initialize while loop
    i <- 1
    while (do_trans < 1) {
        x <- append(x, state1) #record state1
        do_trans <- rbinom(1, size = 1, prob = p[i]) # sample for new transition probability based on p at timestep i
        i <- i + 1
    }
    x <- append(x, state2) # record transition to state2 #this could be a real slow step

    return(x)
}

build_series <- function(p1, p2) { #This function returns a series of states based on transition probability vectors. p1 is a vector of transition probabilities for transitioning from 0 to 1, p2 is a vector of transition probabilities for transitioning from 1 to 2
    state0 <- transition(state1 = 0, state2 = 1, p = p1)
    state1 <- transition(state1 = 1, state2 = 2, p = p2)
    len0 <- length(state0)
    len1 <- length(state1)
    step = c(1:(len0+len1))
    firsttrans <- data.frame(trans_prob = p1[1:length(state0)], state = state0)
    secondtrans <- data.frame(trans_prob = p2[1:length(state1)], state = state1)
    series <- rbind(firsttrans, secondtrans)
    series <- cbind(step, series)
    return(series)
}

# individuals
n_ind <- 100 #number of individuals
ind_effect <- rnorm(n_ind, mean = .05, sd = .005)

simulate_phenodata <- function(n_ind, ind_effect, year) {
    #simulate temperature data
    #temp = rep(5, 50) #constant daily temperature of 5 degrees
    temp <- rnorm(50, mean = 5, sd = 2) #simulate more complex temperature dataset - no trend, but variability
    heatsum <- cumsum(temp)

    #calculate probabilities of transition at different heatsums
    probs1 <- lapply(ind_effect, calc_probability, x=heatsum, k = 0.1, h = 60)
    probs2 <- lapply(ind_effect, calc_probability, x = heatsum, k = 0.1, h = 60)
    #probs1 <- calc_probability(heatsum, k = 0.1, h = 60)
    #probs2 <- calc_probability(heatsum, k = 0.1, h = 120)

    phenofakes <- data.frame()
    for (i in c(1:n_ind)) {
        x <- build_series(probs1[[i]], probs2[[i]])
        #print(c("length", length(x[[1]])))
        indset <- cbind(ind = i, x)
        phenofakes <- rbind(phenofakes, indset)
    }

    phenofakes$heatsum <- heatsum[phenofakes$step]
    phenofakes$year <- year
    return(phenofakes)

}

phenofakes <- lapply(1:15, simulate_phenodata, n_ind = n_ind, ind_effect = ind_effect)
phenofakes <- do.call(rbind, phenofakes)


library(ggplot2)

# ggplot(phenofakes, aes(x = trans_prob, fill = as.factor(state))) +
#     geom_histogram(position = "dodge")


ggplot(phenofakes[which(phenofakes$state < 2),], aes(x = heatsum, y = trans_prob)) +
    geom_jitter() #+
   # stat_function(fun=function(x) 1/(1 + exp(-0.1 * (x - 50))))

phenofakes %>%
    filter(state < 2) %>%
ggplot(aes(x = heatsum[step], y = state)) +
    geom_count() +
    stat_function(fun=function(x) 1/(1 + exp(-0.1 * (x - 60)))) +
    facet_grid(year ~ .) +
    #stat_function(fun=function(x) 1+(1/(1 + exp(-0.1 * (x - 120))))) +
    ggtitle("Simulation of 100 individuals transitioning between phenological states 0, 1, and 2") +
    ylab("Phenological State") +
    xlab("Heat Sum (arbitrary units)")+
    theme_set(theme_gray(base_size = 25))

curve(1/(1 + exp(-0.1 * (x - 60))), from = 0, to = 120, ylab = "p", xlab = "x")
curve(1/(1 + exp(-0.5 * (x - 60))), from = 0, to = 120, add = TRUE)
abline(v = 60, lty = 2)
