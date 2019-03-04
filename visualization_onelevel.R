#mcmc plots for lab meeting

library(bayesplot)
library(tidyr)
library(dplyr)

#functions #################
#logistic function
logistic_b <- function(x, alpha, beta) {
    y <- 1/(1 + exp(alpha - beta * x))
    return(y)
}
# take a dataframe of parameters for the logistic and calculate curves
build_model_curve <- function(param_frame, alpha, heatsum=x){
    y <- lapply(heatsum, logistic_b, alpha=alpha, beta=param_frame$beta)
    names(y) <- x
    y <- t(data.frame(y))
    y <- gather(as.data.frame(y), key=sample, value=value)
    y <- cbind(y, heatsum=x)
}

#calculate transition thresholds for parameter sets
generate_transitions50s <- function(parameter_df, sex) {
    df <- parameter_df %>%
        mutate(h1 = c.1/beta) %>%
        mutate(h2 = c.2/beta) %>%
        mutate(sex = sex) %>%
        gather(key="transition", value ="h", h1, h2)
    return(df)
}

#Male
#--------------
mfit <- readRDS('2019-02-24_onelevel.rds')
bayesplot_theme_update(text=element_text(size=18))
x <- seq(from=0, to=450)
samples <- 200

show(mfit)
mpost <- as.array(mfit)
mparams <- data.frame(rstan::extract(mfit))
mparams_samp <- sample_n(mparams, size=samples)

color_scheme_set("blue")

#plot params
mcmc_dens_overlay(mpost, pars=c("beta", "c[1]", "c[2]")) + ggtitle("MCMC draws - male") + legend_none()



#Female
#-----------
show(ffit)
fpost <- as.array(ffit)
fparams <- data.frame(rstan::extract(ffit))
fparams_samp <- sample_n(fparams, size=samples)

color_scheme_set("green")

#plot params
mcmc_dens_overlay(fpost, pars=c("beta", "c[1]", "c[2]")) + ggtitle("MCMC draws - female") + legend_none()

#model preds
############################
    #calc_female
ft1 <- build_model_curve(fparams_samp, alpha=fparams_samp$c.1) #female transition 1
ft1$phenophase <- 2

ft2 <- build_model_curve(fparams_samp, alpha=fparams_samp$c.2) #female transition 2
ft2$phenophase <- 3

ft <- full_join(ft1, ft2)
ft$index <- group_indices(ft, sample, phenophase)
ft$sex <- "female"



fparams_samp <- generate_transitions50s(fparams_samp, sex="female")

mparams_samp <- generate_transitions50s(mparams_samp, sex="male")

params_samp <- full_join(fparams_samp, mparams_samp)


    #calc male
mt1 <- build_model_curve(mparams_samp, alpha=mparams_samp$c.1) #female transition 1
mt1$phenophase <- 2

mt2 <- build_model_curve(mparams_samp, alpha=mparams_samp$c.2) #female transition 2
mt2$phenophase <- 3

mt <- full_join(mt1, mt2)
mt$index <- group_indices(mt, sample, phenophase)
mt$sex <- "male"

transitions <- rbind(ft, mt)




ggplot(transitions, aes(x=heatsum, y=value, group=index, color=as.factor(phenophase))) +
    geom_vline(data=params_samp, aes(xintercept=h, group=transition), color="grey") +
    geom_line(size=.5) +
    theme_bw(base_size=20) +
    theme(legend.position = "top") +
    scale_color_viridis_d(name="Phenophase") +
    facet_grid(sex ~ .) +
    ggtitle("Transitions for 200 draws from \n MCMC draws for parameters")


#percentiles
