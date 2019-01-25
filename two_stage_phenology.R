library(rethinking)
library(dplyr)
library(tidyr)
library(viridis)


# Ordered logit -----------------------------------------------------------


# Simulate phenology data from an ordered logit model for phenology with 3 states, transitions at heatsums of 60 and 120 and a transition speed of 0.1 (in a model logstic model with the exponent parameterized k(x-h), h is at 60 and 120 and k is 0.1 info and equations used for earlier process simulation
maxtemp = 150
heatsums <- seq(from=1, to = maxtemp, length.out = 1000)
h = c(60,120) # temp thresholds
k = 0.1 # transition speed
phi <- function(x) k*x #x is heatsum
a = k*h

#s1 <- rordlogistic(n=1e4, phi = heatsums , a = a)
s1 <- rordlogit(1e4, phi = phi(heatsums), a = a)
hist(s1)
#p1 <- dordlogistic(s1, phi = heatsums, a = a)
l1 <- dordlogit(s1, phi = phi(heatsums), a = a) # probability of being in a given state at a given heatsum
dens(l1)

p1 <- pordlogit(1:3, a = a, phi = phi(heatsums)) #cum probability of _being_ in a given state
colnames(p1) <- c('p1_s1', 'p1_s2', 'p1_s3')
p1 <- data.frame(heatsums, p1)
dens(p1)

d <- data.frame(heatsums, s1, l1) %>%
    dplyr::inner_join(p1, by = 'heatsums')
d$s1 <- as.factor(d$s1)

dp <- d %>% gather('p1_s1','p1_s2', 'p1_s3', key='state', value='probability')


ggplot(d, aes(x = heatsums, y = s1, color = s1)) +
    geom_point() +
    scale_color_viridis(discrete = TRUE)

ggplot(d, aes(x = heatsums, y = l1, colour = s1)) +
    geom_point() +
    scale_color_viridis(discrete = TRUE)

ggplot(dp, aes(x = heatsums, y = probability, color = state)) +
    geom_point() +
    scale_color_viridis(discrete = TRUE)

ggplot(dp, aes(x = probability, y = s1, color = state)) +
    geom_point() +
    scale_color_viridis(discrete = TRUE)

## Now fit a model to the simulated data

m1 <- map(
    alist(
        s1 ~ dordlogit(phi, c(a1, a2)),
        phi <- b * heatsums,
        b ~ dnorm(0,10), # huge considering b should be equal to k
        c(a1, a2) ~ dnorm(0,10)
    ),
    data = d,
    start= list(a1=1, a2=3, b = 2)
)

post_m1 <- extract.samples(m1)
precis(m1)

plot(200,1, type = "n", xlab='heatsum', ylab='probability', xlim = c(0,200), ylim = c(0, 1))

shortheat <- seq(from = 0, to = 150, length.out = 100)
s <- 1e4
post_m1_samples <- post_m1[sample(1:s, s, replace = FALSE),]
for (s in 1:1000) {
    ak <- as.numeric(p[1:2])
    phi <- post_m1_samples$b[s] * shortheat
    pk <- pordlogit(1:2, a = ak, phi = phi)
    for (i in 1:2) {
        lines(shortheat, pk[,i], col = col.alpha(rangi2,0.1))
        abline(v = h)
    }
}

head(post_m1)
dens(post_m1, show.HPDI = TRUE)

# Now simulate data from the fitted model (sample to simulate predictions)

## get samples from post
samples_m1 <- sample_n(post_m1, 1e3)


#for each set of parameter values, simulate 100 phenology time series


s <- numeric(0)
for (i in 1:length(samples_m1)) {
    si <- apply(as.data.frame(shortheat), 1,
          function(y) rordlogit(100, phi = samples_m1[i,3]*y, a = samples_m1[i,1:2]))# simulate data at all input values (shortheat) at a given parameter set
    si <- cbind(as.numeric(rownames(samples_m1)[i]), si)
    s <- rbind(s, si)
}

s <- data.frame(s)
colnames(s) <- c("param_set", shortheat)
post_pred_m1 <- gather(s, key = heatsum, value = state, -param_set)
post_pred_m1$heatsum <- as.numeric(post_pred_m1$heatsum)

ggplot(post_pred_m1, aes(x = heatsum, y = state, color = as.factor(state))) +
    geom_point() +
    scale_color_viridis(discrete = TRUE)

ggplot(post_pred_m1, aes(factor(state), y = heatsum, fill = as.factor(state))) +
    geom_violin(trim=FALSE) +
    scale_fill_viridis(discrete=TRUE) +
    ggtitle("Data simulated from parameters from fitted model")

ggplot(d, aes(factor(s1), y = heatsums, fill = as.factor(s1))) +
    geom_violin(trim=FALSE) +
    scale_fill_viridis(discrete=TRUE) +
    ggtitle("Data simulated from model")

# Attempt to add levels

# Hierarchical Ordered Logit ----------------------------------------------

maxtemp = 150
heatsums <- seq(from=1, to = maxtemp, length.out = 100)
h = c(60,120) # temp thresholds
k = 0.1 # transition speed
phi <- function(x, ind) k*x + ind #x is heatsum
a = k*h

ind <- rnorm(10, mean=5, sd = 2) #individual effects for 10 individuals

#s1 <- rordlogistic(n=1e4, phi = heatsums , a = a)
s1 <- list()
for (i in 1:length(ind)) {
s1[[i]] <- rordlogit(15, phi = phi(heatsums, ind[i]), a = a)
}
hist(s1)
s1 <- rordlogit(1e4, phi = phi(heatsums, ind = 0), a = a)
#p1 <- dordlogistic(s1, phi = heatsums, a = a)
l1 <- dordlogit(s1, phi = phi(heatsums, ind = 0), a = a) # probability of being in a given state at a given heatsum
dens(l1)

p1 <- pordlogit(1:3, a = a, phi = phi(heatsums)) #cum probability of _being_ in a given state
colnames(p1) <- c('p1_s1', 'p1_s2', 'p1_s3')
p1 <- data.frame(heatsums, p1)
dens(p1)

d <- data.frame(heatsums, s1, l1) %>%
    dplyr::inner_join(p1, by = 'heatsums')
d$s1 <- as.factor(d$s1)

dp <- d %>% gather('p1_s1','p1_s2', 'p1_s3', key='state', value='probability')


ggplot(d, aes(x = heatsums, y = s1, color = s1)) +
    geom_point() +
    scale_color_viridis(discrete = TRUE)

ggplot(d, aes(x = heatsums, y = l1, colour = s1)) +
    geom_point() +
    scale_color_viridis(discrete = TRUE)

ggplot(dp, aes(x = heatsums, y = probability, color = state)) +
    geom_point() +
    scale_color_viridis(discrete = TRUE)

ggplot(dp, aes(x = probability, y = s1, color = state)) +
    geom_point() +
    scale_color_viridis(discrete = TRUE)

m1 <- map(
    alist(
        s1 ~ dordlogit(phi, c(a1, a2)),
        phi <- b * heatsums,
        b ~ dnorm(0,10), # huge considering b should be equal to k
        c(a1, a2) ~ dnorm(0,10)
    ),
    data = d,
    start= list(a1=1, a2=3, b = 2)
)
