#This is code sent to C. Tysor from W. Pearce with J. Davies knowledge. The code was sent in confidence and should NOT be shared (at least before publication)

k.check <- function(x, k=NULL){
    if(is.null(k))
        k <- length(x)
    if(k < 3 )
        stop("Insufficient data or insufficient time window")
    if(k > length(x)){
        warning("'k' longer than 'x'; using all of 'x'")
        k <- length(x)
    }
    return(k)
}

theta.hat <- function(x, k=NULL, late=FALSE){
    k <- k.check(x, k)
    x <- sort(x, FALSE, decreasing=late)
    x <- x[seq_len(k)]
    tryCatch(sum(x * weights(x)), error=function(x) NA)
}

ci.theta.hat <- function(x, alpha=0.05, k=NULL){
    Sl <- function(x,k,alpha) (-log(1 - alpha/2)/k)^(-v.hat(x))
    Su <- function(x,k,alpha) (-log(alpha/2)/k)^(-v.hat(x))

    if(is.null(k))
        k <- length(x)
    if(k > length(x)){
        warning("'k' longer than 'x'; using all of 'x'")
        k <- length(x)
    }
    x <- sort(x, FALSE)
    x <- x[seq_len(k)]


    lower <- tryCatch(x[1] + ((x[1]-x[k]) / (Sl(x,k,alpha) -1)), error=function(x) NA)
    upper <- tryCatch(x[1] + ((x[1] - x[k]) / (Su(x,k,alpha)-1)), error=function(x) NA)
    return(c(lower, upper))
}

#Internals functions for the above
lam <- function(i, j, v){
    std <- function(x,y) (gamma(2*v+x) * gamma(v+y)) / (gamma(v+x) * gamma(y))
    inv <- function(x,y) (gamma(2*v+y) * gamma(v+x)) / (gamma(v+y) * gamma(x))
    return(ifelse(j <= i, std(i,j), inv(i,j)))
}
weights <- function(x){
    k <- length(x)
    v <- v.hat(x)
    lam.mat <- outer(seq_along(x), seq_along(x), lam, v)
    e <- matrix(rep(1,k), ncol=1)
    alpha <- as.vector(solve(t(e) %*% solve(lam.mat) %*% e)) * as.vector(solve(lam.mat) %*% e)
    return(alpha)
}
v.hat <- function(x){
    x <- sort(x, FALSE)
    if(x[1] == x[2]){
        warning("Repeated earliest measurements; applying correction")
        x <- c(x[1], x[x != x[1]])
    }
    if(length(unique(x)) == 2){
        warning("Only two unique measurement dates; unable to compute")
        return(NA)
    }

    k <- length(x)
    return((1/(k-1)) * sum(log((x[1] - x[k]) / (x[1] - x[seq(2,k-1)]))))
}
sim.theta.se <- function(x, k=NULL, n=1000, max.iter=10, late=FALSE){
    #Estimate parameters
    k <- k.check(x, k)
    theta <- theta.hat(x, k=k, late=late)
    shape <- v.hat(x)
    if(is.na(theta))
        return(c(NA,NA))

    #Handle values less than 0 (through recursion)
    if(theta < 0){
        warning("theta < 0; scaling (will unscale before return)")
        output <- Recall(x + max(abs(c(theta, min(x)))), k, n, max.iter)
        return(c(theta, output[2]))
    }

    #Simulate
    sims <- replicate(n, rpois(k, theta^shape))
    s.thetas <- apply(sims, 2, theta.hat, late=late)
    for(i in seq_len(max.iter)){
        s.thetas <- s.thetas[!is.na(s.thetas)]
        if(length(s.thetas) >= n){
            s.thetas <- s.thetas[seq_len(n)]
            break
        }
        sims <- replicate(n, rpois(k, theta^shape))
        s.thetas <- append(s.thetas, apply(sims, 2, theta.hat, late=late))
    }

    if(length(s.thetas) != n){
        warning("Iteration limit exceeded")
        return(NA)
    }
    return(c(theta,sd(s.thetas)))
}
