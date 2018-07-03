library(rethinking)

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

library(viridis)
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
        b ~ dnorm(0,10),
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

precis(m_hs, depth = 2)
logistic(coef(m_hs))
