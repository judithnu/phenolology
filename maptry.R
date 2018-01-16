#map experiment

library(rethinking)

pf <- phenofakes[phenofakes$state < 2, ]

## toy

flist <- alist(
    state ~ dbinom(1,  prob = p),
    p ~ dunif(0,1)
)

mtoy <- map(flist, data = pf)
precis(mtoy)
post <- extract.samples(mtoy)
for (i in 1:length(post)) {
    dens(post[,i])
}

## very simple version of the model. bad because it doesn't "see" individuals
flist <- alist(
    state ~ dbinom(1,  prob = p),
    p <- 1/(1 + exp(-k * (heatsum - h))),
    k ~ dunif(min = .02, max = 0.5),
    h ~ dnorm(mean = 60, sd = 10)
)

m_ex <- map2stan(flist,
                data = pf,
                iter = 4000,
                chains = 4,
                start = list(k = .2, h = 65)
                  )

post <- extract.samples(m_ex)
for (i in 1:length(post)) {
    dens(post[[i]])
}

## logit style

flist <- alist(
    state ~ dbinom(1,  prob = p),
    logit(p) <- k * (heatsum - (h + h_ind[ind])),
    h_ind[ind] ~ dnorm(0, sigma_ind),
    k ~ dunif(min = .01, max = 0.5),
    h ~ dnorm(mean = 55, sd = 10),
    sigma_ind ~ dcauchy(0,1)
)

m_bin <- map2stan(flist,
                 data = pf,
                 iter = 4000,
                 chains = 2
)

post <- extract.samples(m_bin)
total_h_ind <- sapply(1:100, function(ind) post$h + post$h_ind[ , ind])
round(apply(total_h_ind, 2, mean), 2)

dens(post$k)
dens(post$h)
dens(unlist(total_h_ind), show.HPDI = .80)

for (i in 1:length(post)) {
    dens(post[[i]])
}

## slightly more complex version of the model that "sees" individuals



# #attempt to represent model as logit
# flist <- alist(
#     state ~ dbinom(1,p),
#     logit(p) <- a + b*heatsum,
#     a ~ dnorm(mu, sigma),
#     b ~ dnorm(0, 10),
#     mu ~ normal(-5, 10),
#     sigma ~ dnorm(0,10)
# )
#
# m_simple <- map2stan( flist, data = pf, iter = 4000, chains = 2 )
# post <- extract.samples(m_simple)
# a_prob <- logistic(post$a + b*heatsum)
# dens(a_prob)
# dens(b_prob)
# plot(unique(pf$ind), s)
