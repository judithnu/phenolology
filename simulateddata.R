# Simulate data

a <- 0
b <- 1
sigma <- 2
temp = 5
dd <- max(0, temp - threshold)
cum_temp <- cumsum(temp)

log(.2/(1-.2))
log(.8/(1-0.8))

rnorm(100, mean = a + b * cum_temp, sd = sigma)

fauxdat <- data.frame(ind = rep(1:10,2), atemp = c(rnorm(10, 3, 5), rnorm(10, 15, 5)), fl = sort(rep(c(.2,.8),10)))

logit <- glm(fl ~ atemp, family = gaussian(link = "logit"), data = fauxdat)

summary(logit)

plot(fauxdat$atemp, fauxdat$fl)
curve(predict(logit, data.frame(atemp = x), type = "resp"), add = TRUE)

dat <- structure(list(Response = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L,
                                   0L, 0L), Temperature = c(29.33, 30.37, 29.52, 29.66, 29.57, 30.04,
                                                            30.58, 30.41, 29.61, 30.51, 30.91, 30.74, 29.91, 29.99, 29.99,
                                                            29.99, 29.99, 29.99, 29.99, 30.71, 29.56, 29.56, 29.56, 29.56,
                                                            29.56, 29.57, 29.51)), .Names = c("Response", "Temperature"),
                 class = "data.frame", row.names = c(NA, -27L))

temperature.glm <- glm(Response ~ Temperature, data=dat, family=binomial)

plot(dat$Temperature, dat$Response, xlab="Temperature",
     ylab="Probability of Response")
curve(predict(temperature.glm, data.frame(Temperature=x), type="resp"),
      add=TRUE, col="red")
