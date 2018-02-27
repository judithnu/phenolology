library(rethinking)
library(dplyr)
library(lme4)
library(scales)
library(nlme)

#data variants
pf <- phenofakes[phenofakes$state < 2 & phenofakes$heatsum < 120, ]
pf_freq <- pf %>%
    group_by(heatsum) %>%
    summarize(flower_freq = rnorm(n = 1, mean = sum(state)/100, sd = 0.05)) %>%
    filter(flower_freq > 0, flower_freq < 1)

pf_ic3 <- pf[pf$step %% 3 == 0, ] # interval censoring - every 3 days

# get pre and post transition only

pf_temp <- pf %>%
    group_by(ind, state) %>%
    mutate(max_index = max(step)) %>%
    mutate(min_index = min(step))

pre <- pf_temp %>%
    filter(state == 0, step == max_index)

active <- pf_temp %>%
    filter(state ==1, step == min_index)

pf_trans <- rbind(pre, active) %>%
    select(-max_index, - min_index)

#basic logit model

pred_plotter <- function(modeldat, model) { #function to plot data and model predictions from logit
    plot(state~heatsum, data = modeldat)
    points(modeldat$heatsum, fitted(model), col = "red", lwd = 2)
    curve(calc_probability(x, h = 60, k = 0.1), add = TRUE, col = "green")
    #lines(arm::invlogit(state)~heatsum, data = newdat)
    title("basic logit model \n red = model curve, green = source curve")
}

logit <- glm(state ~ heatsum, family = binomial(link = 'logit'), data = pf)
pred_plotter(pf, logit)

# mixed model with random individual effects

pf$indfac <- as.factor(pf$ind)
pf$heatsum_scaled <- pf$heatsum/100


scaled <- scale_params(par = c(.1, 60)) #scaled params

logit2 <- glmer(state ~ heatsum_scaled + (1|indfac), family = binomial, data = pf) #mixed effect model with random individual effect
indcoefs <- fitted(logit2, pf, type = 'response')


pred_plotter2 <- function(modeldat, model, individualeffects) { #function to plot data and model predictions for glmer model
    #newdat <- data.frame(heatsum_scaled = seq(min(modeldat$heatsum_scaled), max(modeldat$heatsum_scaled)), len = 100)
    #newdat$state <- predict(model, newdata = newdat, type = "response")
    plot(jitter(state, amount = 0.02) ~ heatsum_scaled, data = modeldat, col = "red4")
    #lines(state~heatsum_scaled, data = newdat, col = "green", lwd = 2)
    # for (i in 1:length(indcoefs)) {
    #     curve(arm::invlogit(cbind(1,x) %*% fixef(model) + individualeffects[i]), add = TRUE, col = alpha("black", 0.1))
    # }
    curve(arm::invlogit(cbind(1,x) %*% fixef(model)), add = TRUE, col = "red")
    curve(1/(1 + exp(-.1*100 * (x - 60/100))), add = TRUE, col = "green")
    title("binomial mixed model")
}

pred_plotter2(pf, logit2, individualeffects = indcoefs)

calc_probability2 <- function(x, k = steepness, h = midpoint) {
    1/(1 + exp(-k * (x - h + hi)))
}

pf_grouped <- groupedData(state ~ heatsum | ind, pf_trans)

logit3 <- nlme(state ~ 1/(1 + exp(-k * (heatsum - h))),
               data = pf_grouped,
               fixed = h + k ~ 1,
               random = k ~ 1,
               start = c(k = .1, h = 60))

                   nlme(lgcopy ~ logexp2(p1,b1,p2,b2,day),
                        fixed = p1+b1+p2+b2 ~1,random = p1+b1+p2+b2 ~1,
                        data =aids.dat,start=c(start))#multilevel

#map experiment





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


## logit style with individual effects

flist <- alist(
    state ~ dbinom(1,  prob = p),
    logit(p) <- k * (heatsum - (h + h_ind[ind])),
    h_ind[ind] ~ dnorm(0, sigma_ind),
    k ~ dnorm(mean = -.21, sd = 0.1),
    h ~ dnorm(mean = 55, sd = 10),
    sigma_ind ~ dnorm(0,5)
)

m_bin <- map2stan(flist,
                 data = pf,
                 iter = 5000,
                 chains = 4
)

post <- extract.samples(m_bin)
total_h_ind <- sapply(1:100, function(ind) post$h + post$h_ind[ , ind])
round(apply(total_h_ind, 2, mean), 2)

dens(post$k)
dens(post$h)
dens(unlist(total_h_ind), show.HPDI = .80)

probs <- logistic(.45 * (clim$Heatsum - 60))
plot(clim$Heatsum, probs, xlim = c(30,100))

## logit style with individual effects for both h and k

flist <- alist(
    state ~ dbinom(1,  prob = p),
    logit(p) <- (k + k_ind[ind]) * (heatsum - (h + h_ind[ind])),
    h_ind[ind] ~ dnorm(0, sigma_ind),
    k_ind[ind] ~ dnorm(0, sigmak_ind),
    k ~ dnorm(mean = 0.1, sd = 0.1),
    h ~ dnorm(mean = 55, sd = 10),
    sigma_ind ~ dnorm(0,5),
    sigmak_ind ~ dnorm(0, .1)
)

m_bin <- map2stan(flist,
                  data = pf_ic3,
                  iter = 4000,
                  chains = 3
)

post <- extract.samples(m_bin)
total_h_ind <- sapply(1:100, function(ind) post$h + post$h_ind[ , ind])
total_k_ind <- sapply(1:100, function(ind) post$k + post$k_ind[ , ind])
h_sum <- round(apply(total_h_ind, 2, mean), 2)
k_sum <- round(apply(total_k_ind, 2, mean), 2)

dens(post$k)
dens(post$h)
dens(unlist(total_h_ind), show.HPDI = .80)
title("Posterior for h with 80% HDPI")
abline(v = 60, col = "red")
dens(unlist(total_k_ind), show.HPDI = .80)
title("Posterior for k with 80% HDPI")
abline(v = .1, col = "red")

HPDI(post$k, prob = .8)
HPDI(post$h)
probs <- sapply(unique(pf$heatsum), function(x) logistic(k_sum * (x - h_sum)))
pframe <- data.frame(probs)
colnames(pframe) <- unique(pf$heatsum)
pframe$ind <- c(1:100)

library(tidyr)
ppframe <- gather(pframe, key = heatsum, value = prob, -ind)
#probs <- logistic(.45 * (clim$Heatsum - 60))
#plot(clim$Heatsum, probs, xlim = c(30,100))

curve(1/(1 + exp(-0.1 * (x - 60))), from = 0, to = 120, ylab = "p", xlab = "x", col = "red")
for (i in 1:100) {
    curve(1/(1 + exp(-k_sum[i] * (x - h_sum[i]))), add = TRUE, lty = 3)
}

ggplot(ppframe, aes(x = heatsum, y = prob, color = ind)) +
    geom_point() +
    stat_function(fun=function(x) 1/(1 + exp(-k_sum * (x - h_sum))))

plot(unique(pf$heatsum), pframe[,1], xlim = c(0, 120))

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
