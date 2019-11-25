# demo for 11/29 lab meeting presentation
source('phenology_functions.R')
kappa1 <- 20
kappa2 <- 30

beta = 2
x <- seq(from=3, to = 20, length.out = 100)
t1 <- logistic2(x=x, b=beta, c=kappa1)
t2 <- logistic2(x=x, b=beta, c=kappa2)

fauxdat <- data.frame(x, t1, t2) %>%
    pivot_longer(cols=starts_with("t"), names_to = "transition", values_to = "prob")

start <- (logit(0.2) + kappa1)/beta
end <- (logit(0.8) + kappa2)/beta

baseplot <- ggplot(fauxdat, aes(x=x, y=prob, color=transition)) +
    geom_line(size=1, alpha=0.9) +
    scale_color_viridis_d(option="cividis", end=0.8)+
    theme_bw(base_size = 18) +
    theme(legend.position = "none") +
    xlab("Accumulated forcing units") +
    ylab("Transitioned probability")

startandend <- baseplot +
    geom_point(aes(x=start, y=0.2, color="t1"), size=3) +
    geom_point(aes(x=end, y=0.8, color="t2"), size=3)

startandend +
    annotate("rect", xmin=start, xmax=end, ymin=0, ymax=1, alpha=0.2)




