# From Betancourt's case study https://betanalpha.github.io/assets/case_studies/ordinal_regression.html

logit <- function(x) {
    y <- log(x/(1-x))
    return(y)
}


library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())

#Simulate from an ordered logistic model with induced dirichlet prior
simu <- stan(file='dirichlet prior/simulate_ordered_gamma.stan', iter=1, chains=1,
             seed=4838282, algorithm="Fixed_param")

simu_params <- extract(simu)
round(unique(simu_params$c), 5)

input_data <- list("N" = 50, "K" = 5, "y" = array(simu_params$y[1,]))

table(input_data$y)


#Fit the simulated data using the model with a dirichlet prior

fit <- stan(file='dirichlet prior/ordered_logistic_induced.stan', data=input_data, seed=4938483)

#Let's check the model

params <- extract(fit)

#get cutpoints and "beta"
summary(params$c)
summary(params$gamma)

## params
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")



par(mfrow=c(1, 1))

hist(params$gamma, main="", xlab="gamma", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=simu_params$gamma, col="white", lty=1, lw=3)
abline(v=simu_params$gamma, col="black", lty=1, lw=2)

hist(params$c[, 1], breaks=seq(-11, 6, 0.25), main="",
     xlab="Internal Cut Points", xlim=c(-11, 6),
     yaxt='n', ylab="", ylim=c(0, 600),
     col=c_dark, border=c_dark_highlight)
abline(v=simu_params$c[1], col="white", lty=1, lw=3)
abline(v=simu_params$c[1], col="black", lty=1, lw=2)

hist(params$c[, 2], breaks=seq(-11, 6, 0.25),
     col=c_mid_highlight, border=c_dark_highlight, add=T)
abline(v=simu_params$c[2], col="white", lty=1, lw=3)
abline(v=simu_params$c[2], col="black", lty=1, lw=2)

hist(params$c[, 3], breaks=seq(-11, 6, 0.25),
     col=c_mid, border=c_dark_highlight, add=T)
abline(v=simu_params$c[3], col="white", lty=1, lw=3)
abline(v=simu_params$c[3], col="black", lty=1, lw=2)

hist(params$c[, 4], breaks=seq(-11, 6, 0.25),
     col=c_light_highlight, border=c_dark_highlight, add=T)
abline(v=simu_params$c[4], col="white", lty=1, lw=3)
abline(v=simu_params$c[4], col="black", lty=1, lw=2)

## ppc

B <- 5
idx <- rep(1:B, each=2)
x <- sapply(1:length(idx), function(b) if(b %% 2 == 0) idx[b] + 0.5 else idx[b] - 0.5)

obs_counts <- hist(input_data$y, breaks=(1:(B + 1)) - 0.5, plot=FALSE)$counts
pad_obs_counts <- sapply(idx, function(n) obs_counts[n])

pred_counts <- sapply(1:4000, function(n)
    hist(params$y_ppc[n,], breaks=(1:(B + 1)) - 0.5, plot=FALSE)$counts)
probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:B, function(b) quantile(pred_counts[b,], probs=probs))
pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))

plot(1, type="n", main="Posterior Predictive Distribution",
     xlim=c(0.5, B + 0.5), xlab="y",
     ylim=c(0, max(c(obs_counts, cred[9,]))), ylab="")

polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
        col = c_light, border = NA)
polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
        col = c_mid, border = NA)
polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
        col = c_mid_highlight, border = NA)
lines(x, pad_cred[5,], col=c_dark, lwd=2)

lines(x, pad_obs_counts, col="white", lty=1, lw=2.5)
lines(x, pad_obs_counts, col="black", lty=1, lw=2)

