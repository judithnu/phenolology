z = DiscreteWeibull::rdweibull(1000, .9995, 3) %>%
    map(~.x < seq(30)) %>%
    as.data.frame() %>%
    bind_cols()

matplot(z, type = "l", col = "#00000010", lty = 1, lwd = 2)

dodo <- c(1662,1638, 1631, 1628, 1628, 1611, 1607, 1602, 1601, 1598)
dodosurvey <- seq(from = 1598, to = 1662, by = 1)

tdat1 <- c(1:10)
tdat2 <- seq(from = 1, to = 10, by = 2)
tdat3 <- seq(from = 1, to = 10, by = 3)
tdat4 <- seq(from = 1, to = 10, by = 4)

name_seq <- function(interval, label) { #create a dataframe with a sequence of numbers from one to ten at interval defined by "interval" (numeric) and labeled with "label" (character)
    a <- seq(from = 1, to = 10, by = interval)
    b <- rep(label, length(a))
    data.frame(interval = b, sequence = a )
}
intervalcensor <- data.frame(interval = character(0), sequence = numeric(0))

interval <- c(1:4)
label <- as.character(interval)


for (i in c(1:length(interval))) {
    c <- name_seq(interval = interval[i], label = label[i])
    intervalcensor <- rbind(intervalcensor, c)
}

library(dplyr)
tt <- intervalcensor %>%
    dplyr::group_by(interval) %>%
    dplyr::mutate(starttime = theta.hat(sequence))

library(ggplot2)
ggplot(tt, aes(x = interval, y = sequence)) +
    geom_point() +
    geom_point(aes(x = tt$interval, y = tt$starttime), color = "red")

theta.hat(tdat1)
theta.hat(tdat2)
theta.hat(tdat3)
theta.hat(tdat4)

theta.hat(dodo)
theta.hat(dodosurvey)

fauxdat <- c(104,106,108,111,113,115,118)
fauxdat2 <- fauxdat-50

theta.hat(fauxdat)
theta.hat(fauxdat[1:3])


theta.hat(fauxdat2)
theta.hat(fauxdat2[1:3])

#############

p <- 0.5
q <- 0.5
t <- 10 #number of days

transition <- function(p = 0.1, initialstate = "pre", finalstate = "flowering") { #simulate the transition between two states and record the states. State arguments are characters.
    ptrans <- rbinom(1, size = 1, prob = p) #y/n transition from pre to flowering
    x <- c()
    while (ptrans == 0) {
        print(initialstate)
        x <- append(x, initialstate)
        pPF <- rbinom(1, size = 1, prob = p)
    }
    x <- append(x, finalstate)
    return(x)
}
