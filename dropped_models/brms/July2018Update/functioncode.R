dordlogit <- function (x, phi, a, log = FALSE)
{
    a <- c(as.numeric(a), Inf)
    p <- logistic(a[x] - )
        logistic(a[x] - phi)
    na <- c(-Inf, a)
    np <- logistic(na[x] - phi)
    p <- p - np
    if (log == TRUE)
        p <- log(p)
    p
}
dordlogistic <- function (x, phi, a, log = FALSE)
{
    a <- c(as.numeric(a), Inf)
    p <- logistic(a[x] - phi)
    na <- c(-Inf, a)
    np <- logistic(na[x] - phi)
    p <- p - np
    if (log == TRUE)
        p <- log(p)
    p
}


pordlogistic <- function (x, phi, a, log = FALSE)
{
    a <- c(as.numeric(a), Inf)
    if (length(phi) == 1) {
        p <- logistic(a[x] - phi)
    }
    else {
        p <- matrix(NA, ncol = length(x), nrow = length(phi))
        for (i in 1:length(phi)) {
            p[i, ] <- logistic(a[x] - phi[i])
        }
    }
    if (log == TRUE)
        p <- log(p)
    p
}

rordlogistic <- function (n, phi = 0, a)
{
    a <- c(as.numeric(a), Inf)
    k <- 1:length(a)
    if (length(phi) == 1) {
        p <- dordlogit(k, a = a, phi = phi, log = FALSE)
        y <- sample(k, size = n, replace = TRUE, prob = p)
    }
    else {
        y <- rep(NA, n)
        if (n > length(phi)) {
            phi <- rep(phi, ceiling(n/length(phi)))
        }
        for (i in 1:n) {
            p <- dordlogit(k, a = a, phi = phi[i], log = FALSE)
            y[i] <- sample(k, size = 1, replace = TRUE, prob = p)
        }
    }
    y
}