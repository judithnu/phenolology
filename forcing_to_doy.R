# Convert from forcing to DoY

library(rethinking)

ggplot(phendf, aes(x=sum_forcing, y=DoY, color=as.factor(Year))) +
    geom_point() +
    scale_color_viridis_d()

forcingdat <- dplyr::select(phendf, sum_forcing, DoY, Year, Site)
forcingdat$Site <- group_indices(forcingdat, Site)


dayfit <- ulam(
    alist(
        DoY ~ normal(mu, sigma_day),
        mu = (b + by[Year] + bs[Site]) * sum_forcing + alpha,
        b ~ exponential(1.5),
        alpha ~ exponential(1),
        by[Year] ~ normal(0,sigma_year),
        bs[Site] ~ normal(0, sigma_site),
        sigma_site ~ exponential(1.5),
        sigma_year ~ exponential(1.5),
        sigma_day ~ exponential(1.5)
    ),
    data = forcingdat, declare_all_data = FALSE, chains =1, iter=10
)

lm(phendf.DoY ~ phendf.sum_forcing, data=forcingdat)


