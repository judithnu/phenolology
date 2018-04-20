# try logit model with subset of real data
library(dplyr)

rawdat <- read.csv('/home/sus/Documents/research_phd/data/PhenologyAndPollenCounts/from Rita Wagner/data_cleaned/PGTIS_pheno_1997_2012_cleaned.csv', stringsAsFactors = FALSE)
rclim <- read.csv('/home/sus/Documents/research_phd/data/Climate/formatted/PrinceGeorgeSTP.csv', header = TRUE)

# phenology data: one year, male

pdat <- subset(rawdat, Sex == "MALE")
pdat$DayofYear <- lubridate::yday(pdat$Date)

rdat <- subset(rdat, Year == 1997 & Sex == "FEMALE")
rdat$DayofYear <- lubridate::yday(rdat$Date)

# climate data
clim <- subset(rclim, Year %in% unique(pdat$Year))

#calculate amount of heat per day assume no heating below 5 degrees and linear heating starting at 5
no_heat <- clim %>%
    filter(MeanTempC < 5) %>%
    mutate(Heat = 0)

heat <- clim %>%
    filter(MeanTempC >= 5) %>%
    mutate(Heat = MeanTempC - 5)

clim <- rbind(no_heat, heat) %>%
    arrange(Year, DayofYear) %>%
    group_by(Year) %>%
    mutate(Heatsum = cumsum(Heat)) %>% # add heatsum
    select(Year, DayofYear, Heat, Heatsum)

# combine climate and phenology data

df <- merge(pdat, clim) %>%
    select(Year, DayofYear, Clone, Tree, Phenophase, Heat, Heatsum) %>%
    arrange(Year, Tree, DayofYear) %>%
    filter(!Phenophase==0)# drop unexplained 0s

rdf <- merge(rdat, clim) %>%
    select(DayofYear, Clone, Tree, Phenophase, Heat, Heatsum) %>%
    arrange(Tree, DayofYear)

# transform phenology data into 1 (not started), 2 (active), (finished)


by_tree <- group_by(df, Tree, Year)
rby_tree <- group_by(rdf, Tree)

fo <- by_tree %>%
    filter(Phenophase == '4') %>%
    summarise(First_Occurence = min(DayofYear))

rfo <- rby_tree %>%
    filter(Phenophase == '4') %>%
    summarise(First_Occurence = min(DayofYear))

intermed <- merge(df, fo) %>%# first occurance recorded
    arrange(Tree, DayofYear)
rintermed <- merge(rdf, rfo) %>%# first occurance recorded
    arrange(Tree, DayofYear)



intermed_ind <- which(intermed$DayofYear < intermed$First_Occurence) # not started
intermed$Phenophase_Simp <- NA
intermed$Phenophase_Simp[intermed_ind] <- 0
intermed$Phenophase_Simp[is.na(intermed$Phenophase_Simp) == TRUE & intermed$Phenophase == '4'] <- 1
intermed$Phenophase_Simp[is.na(intermed$Phenophase_Simp) == TRUE] <- 2
intermed

rintermed_ind <- which(rintermed$DayofYear < rintermed$First_Occurence) # not started
rintermed$Phenophase_Simp <- NA
rintermed$Phenophase_Simp[rintermed_ind] <- 0
rintermed$Phenophase_Simp[is.na(rintermed$Phenophase_Simp) == TRUE & rintermed$Phenophase == '4'] <- 1
rintermed$Phenophase_Simp[is.na(rintermed$Phenophase_Simp) == TRUE] <- 2
rintermed

df <- intermed %>%
    filter(Phenophase_Simp < 2)
rdf <- rintermed %>%
    filter(Phenophase_Simp < 2)

## males logit
logitm <- glm(Phenophase_Simp ~ Heatsum, family = binomial(link = 'logit'), data = df)

pred_plotter <- function(modeldat, model) { #function to plot data and model predictions from logit
    #newdat <- data.frame(Heatsum = seq(min(modeldat$Heatsum), max(modeldat$Heatsum)), len = 200)
    #newdat$Phenophase_Simp <- predict(model, newdata = newdat, type = "response")
    plot(Phenophase_Simp ~ Heatsum, data = modeldat, col = "red4")
    #lines(Phenophase_Simp ~ Heatsum, data = newdat, col = "red", lwd = 2)
    curve(1/(1 + exp(-(model$coefficients[2]*x + model$coefficients[1]))), add = TRUE, col = "green")
    curve(1/(1 + exp(-0.0926 * (x - 187.174))), col = "red", add = TRUE)
    #lines(arm::invlogit(state)~heatsum, data = newdat)
    title("Prince George 1997 Males \n red = stan binomial, green = glm logit")
}

pred_plotter(df, logitm)

#males

# no individual effects
flist <- alist(
    Phenophase_Simp ~ dbinom(1,  prob = p),
    p <- 1/(1 + exp(-k * (Heatsum - h))),
    k ~ dunif(min = 0, max = 1),
    h ~ dnorm(mean = 150, sd = 50)
)

m_bin <- map2stan(flist,
                  data = df,
                  iter = 4000,
                  chains = 4,
                  start = list(h = 200, k = 0.1)
)


post <- extract.samples(m_bin)
total_h_ind <- sapply(1:100, function(ind) post$h + post$h_ind[ , ind])
round(apply(total_h_ind, 2, mean), 2)

dens(post$k)
dens(post$h)


dens(unlist(total_h_ind), show.HPDI = .80)

for (i in 1:length(post)) {
    dens(post[[i]])
}

ggplot(df, aes(x = Heatsum, y = as.factor(Phenophase_Simp))) +
    geom_count() +
    xlab("Heatsum (Celsius)") +
    ylab("Shedding Pollen?") +
    ggtitle("Pollen shed vs heatsum \n 1997 at Prince George")
#individual effects for h and k

## calculate some priors

priors <- filter(df, DayofYear %in% unique(First_Occurence)) %>%
    summarise(mean(Heatsum), sd(Heatsum))

flist <- alist(
    state ~ dbinom(1,  prob = p),
    logit(p) <- (k + k_ind[ind]) * (heatsum - (h + h_ind[ind])),
    h_ind[ind] ~ dnorm(0, sigmah_ind),
    k_ind[ind] ~ dnorm(0, sigmak_ind),
    k ~ dnorm(mean = 0.1, sd = 0.3),
    h ~ dnorm(mean = priors[1], sd = priors[2]),
    sigmah_ind ~ dunif(0,50),
    sigmak_ind ~ dunif(0,.1)
)

m_bin <- map2stan(flist,
                  data = df,
                  iter = 1e4,
                  warmup = 2e3,
                  chains = 5,
                  start = list(k = .12, h = priors[1]),
                  cores = parallel::detectCores()
)

#females

ggplot(rdf, aes(x = Heatsum, y = as.factor(Phenophase_Simp))) +
    geom_count() +
    xlab("Heatsum (Celsius)") +
    ylab("Receptive?") +
    ggtitle("Receptivity vs heatsum \n 1997 at Prince George")

flist <- alist(
    Phenophase_Simp ~ dbinom(1,  prob = p),
    logit(p) <- k * (Heatsum - (h + h_ind[Tree])),
    h_ind[Tree] ~ dnorm(0, sigma_ind),
    k ~ dnorm(mean = 0.5, sd = 0.25),
    h ~ dnorm(mean = 150, sd = 25),
    sigma_ind ~ dnorm(0,10)
)

m_bin <- map2stan(flist,
                  data = rdf,
                  iter = 100,
                  chains = 1,
                  start = list(h = 150, k = .4, sigma_ind = 10)
)

