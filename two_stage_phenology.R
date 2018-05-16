## This is code from dave
library(brms, quietly = TRUE)
set.seed(1)

phenofakes$state <- phenofakes$state + 1
fit = brm(state ~ heatsum,
    family = cumulative(),
   # prior = set_prior("student_t(3, 0, 25)", class = "b"),
    data = phenofakes,
   chains = 4,
    cores = 4,
    #refresh = 50,
    save_model = "brm.stan"
)

saveRDS(fit, file = "fit.rds")
fit <- readr::read_rds("fit.rds")


plot(fit, pars = c("heatsum"))

graph_probs <- posterior_linpred(fit, transform = TRUE, nsamples = 2000)
