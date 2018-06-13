library(rethinking)
phenofakes$state <- phenofakes$state+1
simplehist(phenofakes$state)


pr_s <- table(phenofakes$state)/nrow(phenofakes)
cum_pr_s <- cumsum(pr_s)
plot(1:3, cum_pr_s, type = "b", xlab="response", ylab="cumulative proportion", ylim=c(0,1))

logit <- function(x) log(x/(1-x))
lco <- logit(cum_pr_s)

m_intercept_only <- map2stan(
    alist(
        state ~ dordlogit(phi, cutpoints),
        phi <- 0,
        cutpoints ~ dnorm(0,10)
          ),
    data=list(state=phenofakes$state),
    start=list(cutpoints=c(-1,0)),
    chains = 2, cores = 2

)

precis(m_intercept_only, depth=2)
logistic(coef(m_intercept_only))

m_hs <- map2stan(
    alist(
        state ~ dordlogit(phi, cutpoints),
        phi <- b*heatsum,
        b ~ dnorm(0,10),
        cutpoints ~ dnorm(0,10)
    ),
    data=phenofakes,
    start=list(cutpoints=c(-1,0)),
    chains = 2, cores = 2
    )

precis(m_hs, depth = 2)
logistic(coef(m_hs))
