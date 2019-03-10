## This is code from dave
library(brms, quietly = TRUE)
set.seed(1)

fit = brm(
    bf(
        cbind("Speech", "Salivation", "Swallowing", "Handwriting", "Cutting", "Dressing",
              "Turning", "Walking", "Climbing", "Respiratory") ~
            elapsed +
            (elapsed | s | subject)
    ) +
        set_rescor(FALSE),
    family = cumulative(),
    prior = set_prior("student_t(3, 0, 25)", class = "b"),
    data = training_data,
    chains = 2,
    cores = 2,
    refresh = 50,
    save_model = "brm.stan"
)

saveRDS(fit, file = "fit.rds")

# me trying to copy Dave's code


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

##########

library(brms, quietly = TRUE)
set.seed(1)
df$Phenophase_Simp <- as.integer(df$Phenophase_Simp)
fit = brm(
    bf(
        Phenophase_Simp ~
            Heatsum + (1|Clone)),
    family = cumulative(),
    prior = set_prior("beta(2,5)", class = "b", coef="Heatsum"),
    data = df,
    chains = 2,
    cores = 2,
    refresh = 50,
    save_model = "brm.stan"
)

saveRDS(fit, file = "fit.rds")
plot(fit)

# add nested tree in clone #USE
set.seed(1)
fit2 = brm(
    bf(
        Phenophase_Simp ~
            Heatsum + (1|Tree/Clone/Orchard)),
    family = cumulative(),
    prior = set_prior("beta(2,5)", class = "b", coef="Heatsum"),
    data = df,
    chains = 5,
    cores = 5,
    refresh = 50,
    iter = 1e4,
    warmup = 2000,
    save_model = "brm.stan"
)

plot(fit2)

set.seed(1)
fit3 = brm(
    bf(
        Phenophase_Simp ~
            Heatsum + (Heatsum+ 1|Tree/Clone/Orchard)),
    family = cumulative(),
    prior = set_prior("beta(2,5)", class = "b", coef="Heatsum"),
    data = df,
    chains = 5,
    cores = 5,
    refresh = 50,
    iter = 3e3,
    warmup = 1000,
    save_model = "brm.stan"
)

plot(fit3)

set.seed(1)
fit4 = brm(
    bf(
        Phenophase_Simp ~
            Heatsum + (Heatsum+ 1|Tree/Clone)),
    family = cumulative(),
    prior = set_prior("beta(2,5)", class = "b", coef="Heatsum"),
    data = df,
    chains = 5,
    cores = 5,
    refresh = 50,
    iter = 3e3,
    warmup = 1000,
    save_model = "brm.stan"
)

plot(fit4)