# Plot some priors
library(tidyr)

beta_pop <- rexp(1000,2)
beta_eff <- rnorm(1000, mean=0, rexp(100, 6))

kappas <- rgamma(1000, 7.5, 1)

fus <- kappas/beta_pop
fu_eff <- kappas/(beta_pop+beta_eff)


priors <- data.frame(beta_pop, kappas, fu_pop = fus, fu_eff = fu_eff) %>%
    gather(key=param, value=value) %>%
    filter(value > -10 & value < 40)
ggplot(priors, aes(x=value)) +
    geom_histogram() +
    facet_wrap("param", scales="free") +
    theme_bw(base_size=18)

min(fu_eff)

length(which(fus > 30))/10000
(length(which(fu_eff >30)) +length(which(fu_eff < -30)))/10000
