# functions supporting priorAnalysisGamma.Rmd




# simulate data from an ordinal logistic model and format it as input for stan. input is a list for stan, as is output. groups is TRUE or FALSE
simulate_data <- function(input, groups) {
    # simulate data
    if (isTRUE(groups)) {
        simu <- rstan::stan(file='dirichlet prior/covar_group_sim.stan', iter=1, chains=1, algorithm="Fixed_param", data=input)
    } else {
        simu <- rstan::stan(file='dirichlet prior/covar_sim.stan', iter=1, chains=1, algorithm="Fixed_param", data=input)
    }
    # extract data from stan model
    simu_params <- rstan::extract(simu)
    # format data as input to another stan model
    input_data_for_model <- list("N" = input$N, "K" = input$K, "x" = input$x, "y" = array(simu_params$y[1,]))

    if (groups==TRUE) {
        append(input_data_for_model, "G"=input$G, "GID" = input$GID)
    }

    return(input_data_for_model)
}

# fit a model with a gamma prior in stan. simdatlist is a list of simulated data (1 simulated dataset per list entry), pars is a list of parameter values (1 set of parameter values per list entry) and groups is TRUE or FALSE indicating whether you're trying to fit groups.
fit_gamma_model <- function(simdatlist, pars, groups) {
    #choose whether to use data simulated with a rapid or slow transition
    if (pars$transition == "slow") {
        simdat <- simdatlist$slow
    }
    if (pars$transition == "medium") {
        simdat <- simdatlist$medium
    }
    if (pars$transition == "fast") {
        simdat <- simdatlist$fast
    }
    #extract parameters for prior distribtuions
    simdat$shape <- pars$shape
    simdat$cut_rate <- pars$cut_rate
    simdat$beta_rate <- pars$beta_rate


    #fit the model
    if (isTRUE(groups)) {
        fitgam <- stan(file='dirichlet prior/gamma/gamma_covar_group.stan', data=simdat, chains=4, iter = 3500, warmup=1000)
    } else {
        fitgam <- stan(file='dirichlet prior/gamma/gamma_covar.stan', data=simdat, chains=4)
    }
    return(fitgam)
}

# fit model with induced dirichlet prior
fit_indir_model <- function(simdatlist, pars, groups) {
    #choose whether to use data simulated with a rapid or slow transition
    if (pars$transition == "slow") {
        simdat <- simdatlist$slow
    }
    if (pars$transition == "medium") {
        simdat <- simdatlist$medium
    }
    if (pars$transition == "fast") {
        simdat <- simdatlist$fast
    }

    #extract parameters for prior distribtuions
    simdat$beta_rate <- pars$beta_rate
    simdat$anchor <- pars$anchor

    #fit the model

    fit <- stan(file='dirichlet prior/dirichlet/dirichlet_covar.stan', data=simdat, chains=4)
    return(fit)
}

## append a label (string) to all columnnames in a dataframe (x)
label_names <- function(x, label) {
    colnames(x) <- paste0(colnames(x), "_", label)
    return(x)
}

#bind true parameters (in list parlist) and model parameters (in list fit) even tho it will make a giant df.
bind_true_model_pars <- function(fits, parlist) {
    # extract params from model object
    params <- lapply(fits, function(x) {data.frame(rstan::extract(x) ) } )

    # label params as coming from the model or as true params used to
    params <- map(params, label_names, label="model")
    parlist <- map(parlist, label_names, label="true")

    # combine model and true params in a big list of dataframes - each list entry is a dataframe for a single model
    params <- map2(params, parlist, cbind)

    return(params)
}

## append a label (string) to all columnnames in a dataframe (x)
label_names <- function(x, label) {
    colnames(x) <- paste0(colnames(x), "_", label)
    return(x)
}

# # function to plot modeled parameters with label being a string to label the graph (usually prior and whether groups were included and pardf is the modeled parameter dataframe) trueparams is a one row dataframe of true parameters. h is a global parameter have fun!
# parplot <- function(pars) {
#     cutplot <- mcmc_intervals(pars, regex_pars="c.\\d_model") +
#         geom_vline(xintercept=c(unique(pars$c.1_true), unique(pars$c.2_true)))
#
#     betaplot <- mcmc_intervals(pars, pars="beta_model") +
#         geom_vline(xintercept=unique(pars$beta_true))
#
#     h1plot <- mcmc_intervals(pars, regex_pars = "h1_model") +
#         geom_vline(xintercept=unique(pars$h1_true))
#     h2plot <- mcmc_intervals(pars, regex_pars = "h2_model") +
#         geom_vline(xintercept=unique(pars$h2_true))
#     allplot <- cowplot::plot_grid(cutplot, betaplot, h1plot, h2plot,
#                                   nrow=1, ncol=4,
#                                   labels=paste("model", unique(pars$modelid_true),
#                                                "rate=", unique(pars$beta_rate_true),
#                                                "shape=", unique(pars$shape_true)))
#     print(allplot)
# }

# function to calculate the difference between modeled parameters (in the model_params dataframe) and true parameters (h is globally declared). Where true params is a dataframe
posterior_differencer <- function(pars, h) {
    c1_diff <- pars$c.1_model - pars$c.1_true
    c2_diff <- pars$c.2_model - pars$c.1_true
    h1_diff <- pars$h1_model - pars$h1_true
    h2_diff <- pars$h2_model - pars$h2_true
    beta_diff <- pars$beta_model - pars$beta_true
    diffframe <- data.frame(c1_diff, c2_diff, h1_diff, h2_diff, beta_diff)
    return(diffframe)
}

# function to plot histograms of differences between true params and modeled params (model_params dataframe). MAKE IT WORK WITH BOTH GAM AND INDIR
diffplotter <- function(diffs, pars) {
    #diffs <- posterior_differencer(pars)
    cuts <- mcmc_intervals(diffs, regex_pars = "c") +
        ggtitle("", subtitle = "differences between modeled and true params")
    opars <- mcmc_intervals(diffs, pars=c("h1_diff", "h2_diff", "beta_diff"))
    cowplot::plot_grid(cuts, opars, labels=paste("model", unique(pars$modelid_true),
                                                 "rate=", unique(pars$beta_rate_true),
                                                 "anchor=", unique(pars$anchor_true)))
}

HPDIlow <- function(x, prob) {
    HPDI <- rethinking::HPDI(x, prob=prob)
    return(HPDI[1])
}

HPDIhigh <- function(x, prob) {
    HPDI <- rethinking::HPDI(x, prob=prob)
    return(HPDI[2])
}


calc_HPDI <- function(params, prob) {
    low <- params %>% dplyr::summarise_at(vars(ends_with("model")), HPDIlow, prob=prob)
    high <- params %>% dplyr::summarise_at(vars(ends_with("model")), HPDIhigh, prob=prob)

    # awkward formatting
    hdpis <- dplyr::full_join(low, high) %>%
        select(-contains("lp"))
    colnames(hdpis) <- stringr::str_replace(colnames(hdpis), "_model", "")

    # true param
    true <- params %>% dplyr::summarise_at(vars(ends_with("true")), unique)
    colnames(true) <- stringr::str_replace(colnames(true), "_true", "")
    true <- select(true, colnames(hdpis))

    # more awkward formatting
    compframe <- dplyr::full_join(hdpis, true) %>%
        t(.) %>%
        data.frame()
    colnames(compframe) <- c("low", "high", "true")
    compframe$params <- rownames(compframe)

    # true param in interval?
    tf <- compframe %>% mutate(inint = true > low & true < high)
    return(tf)

}


