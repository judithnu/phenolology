# functions for prior analysis

# generate samples from a truncated normal distribution.
# n = how many samples, mean = mean, sd = sd. min and max are the limits of the distribution/truncation points.
rtnorm <- function(n, mean, sd, min, max) {
    x <- rnorm(n, mean=mean, sd=sd)
    x <- x[x >= min & x <= max]
    while(length(x) < n) {
        newx <- rnorm(1, mean=mean, sd=sd)
        while(newx <= min | newx >= max) {
            newx <- rnorm(1, mean=mean, sd=sd)
        }
        x <- c(x, newx)
    }
    length(x)==n
    return(x)
}

# Build a list of objects that are used for simulating data in stan. simu_pars is a dataframe of beta and cutpoint parameters, inputlist is what gets fed directly to stan, and h is the transition points and is useful for plotting later
set_simulation_parameters <- function() {
    N <- 50
    K <- 3

    beta <- data.frame(transition = c("slow", "medium", "fast"), beta=c(0.5, 1, 2))
    cutpoints <- data.frame(c.1= c(4,8,16), c.2=c(6, 12, 24), transition=c("slow", "medium", "fast"))

    simu_pars <- merge(beta, cutpoints)

    # half transition points, engineered to be identical all transitions
    h1 <- unique(simu_pars$c.1/simu_pars$beta )
    h2 <- unique(simu_pars$c.2/simu_pars$beta )
    h <- c(h1, h2)

    # covariate over full range of heat accumulation (risto scale) Jan-Julyish
    #x <- rtnorm(n=N, mean=mean(h), sd=2, min=0, max=20) #covariate
    x <- runif(n=N, min=0, max=20)

    inputs_for_sim <- split(simu_pars, simu_pars$transition) %>%
        purrr::map(.f=function(y) {list("N" = N, "K" = K, "c" = c(y$c.1, y$c.2), "beta"=y$beta, "h" = h, "x" = x)})

    return(list(pars=simu_pars, inputlist=inputs_for_sim, h=h))
}

# Plot simulated data to make sure it's reasonable
plot_simulated_data <- function(simdat, simulation_input) {
    simdf <- purrr::map(simdat, .f = function(x) {x[c("x", "y")]}) %>%
        purrr::map_dfr(.f = bind_rows, .id=".id")

    p1 <- ggplot(simdf, aes(x=x, y=y)) +
        geom_jitter(shape=1, height=0.1, alpha=0.5) +
        geom_vline(xintercept = simulation_input[[1]]$h) +
        ggtitle("Simulated data with cutpoints") +
        facet_grid(.id ~ .)

    p2 <- ggplot(simdf, aes(x=x, colour=as.factor(y))) +
        stat_ecdf() +
        geom_vline(xintercept=simulation_input[[1]]$h) +
        theme(legend.position = "none") +
        ggtitle("Cumulative x for 3 states") +
        facet_grid(.id ~ .)

    cowplot::plot_grid(p1, p2, ncol=2)
}

# turn a dataframe of parameter values into a list of parameters for each model run. Each row of the dataframe becomes a dataframe in a list
make_parframe_list <- function(parframe) {
    parlist <- split(parframe, seq(nrow(parframe)))
    names(parlist) <- parframe$modelid

    return(parlist)
}
