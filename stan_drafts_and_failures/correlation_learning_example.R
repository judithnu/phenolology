#mcelreath
a <- 3.5
b <- (-1)
sigma_a <- 1
sigma_b <- 0.5
rho <- (-0.7)
Mu <- c( a , b )
cov_ab <- sigma_a*sigma_b*rho
Sigma <- matrix( c(sigma_a^2,cov_ab,cov_ab,sigma_b^2) , ncol=2 )
matrix( c(1,2,3,4) , nrow=2 , ncol=2 )
sigmas <- c(sigma_a,sigma_b) # standard deviations
Rho <- matrix( c(1,rho,rho,1) , nrow=2 ) # correlation matrix
# now matrix multiply to get covariance matrix
Sigma <- diag(sigmas) %*% Rho %*% diag(sigmas)
N_cafes <- 20
library(MASS)

vary_effects <- mvrnorm( N_cafes , Mu , Sigma )
a_cafe <- vary_effects[,1]
b_cafe <- vary_effects[,2]

N_visits <- 10
afternoon <- rep(0:1,N_visits*N_cafes/2)
cafe_id <- rep( 1:N_cafes , each=N_visits )
mu <- a_cafe[cafe_id] + b_cafe[cafe_id]*afternoon
sigma <- 0.5 # std dev within cafes
wait <- rnorm( N_visits*N_cafes , mu , sigma )
d <- data.frame( cafe=cafe_id , afternoon=afternoon , wait=wait )
R <- rlkjcorr( 1e4 , K=2 , eta=2 )
dens( R[,1,2] , xlab="correlation" )

m14.1 <- ulam(
    alist(
        wait ~ normal( mu , sigma ),
        mu <- a_cafe[cafe] + b_cafe[cafe]*afternoon,
        c(a_cafe,b_cafe)[cafe] ~ multi_normal( c(a,b) , Rho , sigma_cafe ),
        a ~ normal(5,2),
        b ~ normal(-1,0.5),
        sigma_cafe ~ exponential(1),
        sigma ~ exponential(1),
        Rho ~ lkj_corr(2)
    ) , data=d , sample=FALSE )

stancode(m14.1)

#mine

helper <- ulam(
    alist(
        Phenophase_Derived ~ ordered_logistic(phi, cutpoints),
        phi <- beta[SexID]*forcing_accum + a_sex[SexID],
        # adaptive prior
        c(beta, a_sex)[SexID] ~ multi_normal(c(a,b), Rho, sigma_sex),
        # fixed prior
        a ~ dbeta(.5,5),
        b ~ normal(0,1),
        sigma_sex ~ exponential(1),
        cutpoints ~ dnorm(15,5),
        Rho ~ lkj_corr(2)
    ),
    data=phendf , chains=1 , cores=1, iter=20 )

#can't get this model to work with other effects

# Too many indexes, expression dimensions=0, indexes found=1
# error in 'model45919501e97_21f79605f2e4ee2f4ceb12392211c13a' at line 42, column 95
# -------------------------------------------------
#     40:     }
# 41:     for ( i in 1:22416 ) {
#     42:         phi[i] = beta[SexID[i]] * forcing_accum[i] + a_sex[SexID[i]] + a_prov[ProvenanceID[i]];
#     ^
#         43:     }