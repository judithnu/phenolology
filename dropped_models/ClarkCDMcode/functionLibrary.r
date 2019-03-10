

library(mvtnorm)

############################### map species

mapSpecies <- function(specs,x,y,z,mapx=range(x),mapy=range(y),
                       scale=0,add=F){
	
	specNames <- sort(unique(specs))
	ns        <- length(specNames)
	cc        <- as.vector(sample(c(1:100),ns))
	cvec      <- 1
	if(ns > 1)cvec <- match(specs,specNames)
	
	if(scale > 0)mapsetup(mapx,mapy,scale)
	symbols(x,y,circles=z/10,inches=F,xlim=mapx,ylim=mapy,fg=cc[cvec],add=add)
}

appendData <- function(oldfile,newfile,oldDates,newDates){  #append new plot data
    	
    rold <- rownames(oldfile)
    rnew <- rownames(newfile)

    	if(!is.matrix(newfile)){
    		newfile <- matrix(newfile,length(newfile),1)
    		colnames(newfile) <- newDates
    	}

    	if(length(oldfile) == 0)return(newfile)
    	
    	if(!is.matrix(oldfile)){
    		oldfile <- matrix(oldfile,length(oldfile),1)
    		colnames(oldfile) <- oldDates
    	}
    	
      allDates <- unique(c(oldDates,newDates))
    	
      newMat  <- matrix(NA,nrow(newfile),length(allDates))
      colnames(newMat) <- allDates

      wwc     <- match(colnames(newfile), allDates)
      newMat[,wwc] <- newfile
      newMat  <- matrix(newMat,nrow(newfile),length(allDates))
      colnames(newMat) <- allDates
      rownames(newMat) <- rnew

      oldMat <- matrix(NA,nrow(oldfile),length(allDates))
      colnames(oldMat) <- allDates
      wwc    <- match(colnames(oldfile), allDates)
      oldMat[,wwc] <- oldfile
      oldMat <- matrix(oldMat,nrow(oldfile),length(allDates))
      colnames(oldMat) <- allDates
      rownames(oldMat) <- rold
      
      rbind(oldMat,newMat)
}
#######################################
appendMatrix <- function(c1,c2){

   if(length(c1) == 0)return(c2)
   w1 <- match(colnames(c1),colnames(c2))
   wn1 <- which(is.na(wc))
   w2 <- match(colnames(c2),colnames(c1))
   wn2 <- which(is.na(wc))

   c22 <- matrix(NA,nrow(c2),length(w1))
   c22[,w2] <- c2
   rownames(c22) <- rownames(c2)
   colnames(c22) <- colnames(c1)
   rbind(c1,c22)
}
############################################
row2Mat <- function(vec){

  if(is.matrix(vec))return(vec)
  vn  <- names(vec)
  vec <- matrix(vec,1)
  colnames(vec) <- vn
  vec
}
############################################
col2Mat <- function(vec,namecol=NULL){

  if(is.matrix(vec))return(vec)
  vn  <- names(vec)
  vec <- matrix(vec,ncol=1)
  rownames(vec) <- vn
  colnames(vec) <- namecol
  vec
}
##############################################
rowBind <- function(matnow,row2add,rowName){

  if(length(matnow) == 0){
    matnow <- row2add
    if(!is.matrix(row2add)){
        matnow <- matrix(matnow,nrow=1)
        if(!is.null(names(row2add)))colnames(matnow) <- names(row2add)
    }
    rownames(matnow) <- rowName
    return(matnow)
  }

  matnow <- rbind(matnow,row2add)
  rownames(matnow)[nrow(matnow)] <- rowName
  matnow
}
#############################################

values2contour <- function(x,y,z,q,col='black',lwd=1,zlevs=NULL,add=F){    

  #contours where x,y is not a uniform grid, requires 'spatial' library

  require(spatial)

  xr    <- range(x,na.rm=T)
  yr    <- range(y,na.rm=T)
  xd <- (xr[2] - xr[1])/q
  yd <- (yr[2] - yr[1])/q

  xgrid <- seq(xr[1],xr[2],by=xd)
  ygrid <- seq(yr[1],yr[2],by=yd)

  zsurf <- surf.gls(2,expcov,x,y,z,nx=1000,d=10)

  prsurf <- prmat(zsurf, xr[1],xr[2],yr[1],yr[2], xd)

  
  if(is.null(zlevs))zlevs <- signif(seq(min(z), max(z), length=nlevs),1)
  contour(prsurf, levels=zlevs,lwd=lwd,col=col,add=T)
}
#####################################################

crosscorByRow <- function(xmat,ymat,lag=ncol(xmat),BOOTSTRAP=F,PLOT=F){  
  #cross correlation for each row of xmat[i,] vs ymat[i,]

  nn <- nrow(xmat)
  xx <- c(-lag:lag)
  nc <- length(xx)
  yy <- matrix(NA,n,nc)
  ii <- numeric(0)

  ciMean <- numeric(0)

  for(i in 1:nn){
    di <- xmat[i,]
    fi <- ymat[i,]
    wi <- which(is.finite(fi) & is.finite(di) & fi > 0)
    if(length(wi) < 5)next
    if(var(fi[wi]) == 0)next

    xxx <- xx[wi]
    ddd <- di[wi]
    fff <- fi[wi]
    di  <- lm(ddd ~ xxx)$residuals
    fi  <- lm(fff ~ xxx)$residuals

    cross <- ccf(di,fi,type='correlation',plot=F)
    cx <- cross$lag
    cy <- cross$acf

    cy <- cy[cx %in% xx]
    cx <- cx[cx %in% xx]
    yy[i,match(cx,xx)] <- cy
    ii <- c(ii,i)
  }

  ci <- apply(yy,2,quantile,c(.5,.025,.975),na.rm=T)
  colnames(ci) <- xx

  nk <- length(ii)   #sample size for good series

  if(BOOTSTRAP){

    nboot <- 2000
    mu <- matrix(NA,nboot,nc)
    for(g in 1:nboot){
      isamp <- sample(ii,nk,replace=T)
      mu[g,] <- apply(yy[isamp,],2,mean,na.rm=T)
    }
  }

  ciMean <- apply(mu,2,quantile,c(.025,.975),na.rm=T)
  colnames(ciMean) <- xx

  if(PLOT){
    par(bty='n')
    plot(xx,ci[1,],type='l',lwd=2,ylim=c(-.6,.6),ylab='Correlation',xlab='Lag',col=2)
    abline(h=0,lwd=2,col='grey')
    abline(v=0,lwd=2,col='grey')

    for(j in 1:3)lines(xx,ci[j,],lty=2)
    for(j in 1:2)lines(xx,ciMean[j,],lty=2,col=2,lwd=2)

    text(xx[1],.5,paste('n = ',nk),pos=4)
  }

  list(lag = xx, ci = ci,  ciMean = ciMean, n = nk)
}


####################################################

points2contour <- function(x,y,q,xlabel,ylabel,main,xp,yp,levs){  

  #creates contours for (x,y) at density q

  xr    <- range(x,na.rm=T)
  yr    <- range(y,na.rm=T)
  xd <- (xr[2] - xr[1])/q
  yd <- (yr[2] - yr[1])/q

  xgrid <- seq(xr[1],xr[2],by=xd)
  ygrid <- seq(yr[1],yr[2],by=yd)

  xf <- cut(x,xgrid)
  yf <- cut(y,ygrid)

  z <- table(xf,yf)

  xmids <- (xgrid - xd/2)[-1]
  ymids <- (ygrid - yd/2)[-1]

  lwdd <- seq(1,length(levs),by=1)

  image(xmids,ymids,z,xlab=xlabel,ylab=ylabel,xlim=c(xp[1],xp[2]),ylim=c(yp[1],yp[2]))
  contour(xmids,ymids,z,add=T,levels=levs,lwd=lwdd)
  title(main)
}

####################################################

histf <- function(vec,minv,maxv)hist(vec,breaks=seq(minv,maxv,by=.02))


####################################################

myrmultinom <- function(size,p){  

  #n multinomial r.v. for a n by ncol(p) matrix of probs
  #each row of p is a probability vector

  n     <- nrow(p)
  J     <- ncol(p)
  y     <- matrix(0,n,J)
  sizej <- rep(size,n)
  sumj  <- rep(0,n)
  dpj   <- rep(1,n)
  pj    <- p
  wj    <- c(1:n)

  for(j in 1:(J-1)){
    a     <- round(pj[wj,1],10)
    y[wj,j] <- rbinom(length(wj),sizej[wj],a)
    sumj  <- sumj + y[,j]
    sizej <- size - sumj
    dpj   <- dpj - p[,j]
    pj    <- matrix(p[,c((j+1):J)]/dpj,nrow(p))
    wj    <- which(sumj < size,arr.ind=T) 
  }

  if(n == 1)y[,J] <- size - sum(y)
  if(n > 1) y[,J] <- size - apply(y,1,sum)
  y

}

####################################################

truncpars <- function(x,lo,hi){        

  #JS Clark
  #fit truncated multivariate normal to posterior x, known lo and hi

  if(is.vector(x))x <- matrix(x,length(x),1)  
  px <- ncol(x)          #dimension of original x

  ww <- which(!is.finite(x),arr.ind=T)
  if(length(ww) > 0){
    ww <- unique(ww[,2])
    wk <- which(!c(1:ncol(x)) %in% ww,arr.ind=T)
    x <- x[,wk]
  }
  muvec <- apply(x,2,mean,na.rm=T)
  sig   <- cov(x,use="complete.obs")
  pk <- ncol(x)
  nn <- nrow(x)

  ngg   <- 1000
  nkeep <- ngg - 300
  nk    <- 0
  mug   <- rep(0,pk)
  cvg   <- rep(0,pk^2)

  for(g in 1:ngg){
    tmp   <- truncmvtnorm(x,muvec,sig,lo,hi)
    muvec <- tmp$mu
    sig   <- tmp$sig

    if(g > nkeep){
      print(muvec)
      nk       <- nk + 1
      mug      <- mug + muvec
      cvg      <- cvg + as.vector(sig)
    }
  }
  mvec <- mug/nk
  covmat <- matrix(cvg/nk,pk,pk)

  if(length(ww) > 0){
    mnew     <- rep(NA,px)
    mnew[wk] <- mvec
    covnew   <- matrix(NA,px,px)
    covnew[wk,wk] <- covmat
    mvec   <- mnew
    covmat <- covnew
  }

  list(mu = mvec, cm = covmat)
}

####################################################

trunclogis <- function(n,lo,hi,bpars,xvars,wp){ 

  #truncated logistic
  # bpars - parameter vector for logit
  # xvars - variables corresponding to bpars, e.g., c(1,x1,x2)
  # wp    - which variable and parameter in the logit model 
  #         (not 1, because bpars[1] is intercept)
  #         (xvars[wp] is not used)

  xlo     <- xvars
  xlo[wp] <- lo
  xhi     <- xvars
  xhi[wp] <- hi
  sl      <- sum(xlo*bpars)
  sh      <- sum(xhi*bpars)
  lflo    <- exp(sl)/(1 + exp(sl))
  lfhi    <- exp(sh)/(1 + exp(sh))

  z <- runif(n,lflo,lfhi)
  (log(z/(1 - z)) - sum(xvars[-wp]*bpars[-wp]))/bpars[wp]
}

####################################################

truncmvtnorm <- function(x,muvec,sig,lo,hi){

  #sample from a truncated normal

  if(is.vector(x))x <- matrix(x,length(x),1)
  if(length(sig) == 1)sig <- matrix(sig,1,1)
  y  <- x*0
  pk <- ncol(x)
  n  <- nrow(x)

  for(j in 1:pk){

    if(j == 1){
      t <- muvec[1]
      w <- sqrt(sig[1,1])
    }
    if(j > 1){

       svec <- c(1:(j-1))
       vmat <- sig[j,svec] %*% solve(sig[svec,svec])
       ymu  <- y[,svec] - matrix(rep(muvec[svec],n),n,(j-1),byrow=T)
       t    <- t(muvec[j] +  vmat %*% t(ymu))
       w    <- as.numeric(sqrt( sig[j,j] - vmat %*% c(sig[svec,j]) ))
    }

    up <- pnorm(x[,j],t,w) - pnorm(lo[j],t,w)
    do <- pnorm(hi[j],t,w) - pnorm(lo[j],t,w)

    add <- w*qnorm(up/do)
    add[!is.finite(add)] <- 0
    y[,j] <- t + add

  }

  muy   <- apply(y,2,mean)
  muvec <- myrmvnorm(1,muy,sig/n)

 #the covariance matrix

  mumat <- matrix(muvec,n,pk,byrow=T)
  sy    <- crossprod(y - mumat) 
   ss   <- solve(sy)
   df   <- n 
  alpha <- myrmvnorm(df,rep(0,pk),ss)
  th    <- solve(crossprod(alpha))

  list(mu = muvec, sig = th)

} 

####################################################

logit <- function(x){log(x/(1-x))}  #logit

####################################################

invlogit <- function(x, log = FALSE){  #inverse logit

  if(log)return(-log(1 + exp(-x)))
  1/(1 + exp(-x))

}
##################################################
invlogt <- function(pc,lss){
  z <- exp(pc[1] + pc[2]*lss)
  z/(1 + z)
}
##############################
invMat <- function(s,NEARPD=F){  #matrix inversion, if NEARPD find closest PD matrix

    testv <- try(chol(s),T)

    if(inherits(testv,'try-error')){
       message('error in invMat')
       if(NEARPD){
         require(Matrix)
         s     <- as.matrix(nearPD(s)$mat)
         testv <- try(chol(s),T)
       }
    }

    chol2inv(testv)
}

################################################
mydmvnorm <- function(x,mean,sigma,log=FALSE){

  #mv normal density

    if (is.vector(x))x <- matrix(x, ncol = length(x))
    if (is.vector(mean))mean <- matrix(mean, ncol = length(x))

    xx   <- sweep(x, 2, mean)

    testv <- try(chol(sigma),T)
    if(inherits(testv,'try-error')){
       message(paste('error in mydmvnorm at g = ',g,sep=''))
       return(mean)
    }

    cov <- chol2inv(testv)
    distval <- rowSums((xx %*% cov) * xx)
    names(distval) <- rownames(xx)

  #  distval <- mahalanobis(xx, mean, sigma)
    logdet   <- sum(log( eigen(sigma, symmetric = TRUE, only.values = TRUE)$values ))
    if(is.na(logdet))return(logdet)
    logretval <- -(ncol(xx) * log(2 * pi) + logdet + distval)/2
    if(log)return(logretval)
    exp(logretval)
}


####################################################

myrmvnorm <- function (n, mu, sigma){

#    ev <- eigen(sigma, sym = TRUE)$values
#    if (!all(ev >= -sqrt(.Machine$double.eps) * abs(ev[1])))
#        warning("sigma is numerically not positive definite")
    sigsvd <- svd(sigma)
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval + mu
}


####################################################

mypmvnorm <- function (lower, upper, mu,sigma){

   corr <- cov2cor(sigma)
   lower <- (lower - mu)/sqrt(diag(sigma))
   upper <- (upper - mu)/sqrt(diag(sigma))
   mean <- rep(0, length(lower))
   RET  <- mvt(lower = lower, upper = upper, df = 0,
                corr = corr, delta = mu, maxpts = 25000, abseps = 0.001,
                releps = 0)
   return(RET$value)
}
####################################################

myqmvnorm <- function (p, interval = c(-10, 10), tail = c("lower.tail", "upper.tail",
    "both.tails"), mean = 0, corr = NULL, sigma = NULL, maxpts = 25000,
    abseps = 0.001, releps = 0, ...)
{
    if (length(p) != 1 || (p <= 0 || p >= 1))
        stop(sQuote("p"), " is not a double between zero and one")
    tail <- match.arg(tail)
    dim <- length(mean)
    if (is.matrix(corr))
        dim <- nrow(corr)
    if (is.matrix(sigma))
        dim <- nrow(sigma)
    lower <- rep(0, dim)
    upper <- rep(0, dim)
    args <- checkmvArgs(lower, upper, mean, corr, sigma)
    dim <- length(args$mean)
    pfct <- function(q) {
        switch(tail, both.tails = {
            low <- rep(-abs(q), dim)
            upp <- rep(abs(q), dim)
        }, upper.tail = {
            low <- rep(q, dim)
            upp <- rep(Inf, dim)
        }, lower.tail = {
            low <- rep(-Inf, dim)
            upp <- rep(q, dim)
        }, )
        pmvnorm(lower = low, upper = upp, mean = args$mean, corr = args$corr,
            sigma = args$sigma, abseps = abseps, maxpts = maxpts,
            releps = releps) - p
    }
    if (tail == "both.tails") {
        interval[1] <- 0
        interval <- abs(interval)
    }
    qroot <- uniroot(pfct, interval = interval, ...)
    names(qroot)[1:2] <- c("quantile", "f.quantile")
    qroot
}


####################################################

fit.tnorm <- function(parvec){  

  #fit parameters for truncated normal
  #vector parvec includes mean, diagonal, & offdiagonals for cov matrix

  mu <- parvec[mi]
  sigmat <- matrix(0,length(mu),length(mu))
  sigmat <- diag(parvec[p2])
  sigmat[si] <- parvec[p3]
  sigmat[cbind(si[,2],si[,1])] <- parvec[p3]

  c1     <- mydmvnorm(x,mu,sigmat,log=T)
  c2     <- log(mypmvnorm(lo,hi,mu,sigmat))
  -sum(c1 - c2)

}

####################################################

smooth.na <- function(x,y){   

  #remove missing values

    wy <- which(!is.finite(y),arr.ind =T)
    if(length(wy) == 0)return(cbind(x,y))
      ynew <- y[-wy]
      xnew <- x[-wy]

    return(cbind(xnew,ynew))
}
####################################################

smooth.ma <- function(y,wt){   

  #moving average filter with weights (w0,w1,...), assumed symmetric

  if(length(wt) > length(y))wt <- wt[1:length(y)]
  nw <- length(wt)
  ny <- length(y)
  w <- c(rev(wt),wt[2:nw])
  ymat <- matrix(NA,ny,length(w))
  kb <- nw
  ke <- ny
  ky <- ny - kb + 1

  kb <- c(nw:-(nw-2)); kb[kb < 1] <- 1
  ke <- c((ny+nw-1): (ny-nw+1)); ke[ke > ny] <- ny
  yb <- rev(kb)
  ye <- rev(ke)

  for(kj in 1:(2*nw-1)){
    ymat[kb[kj]:ke[kj],kj] <- y[yb[kj]:ye[kj]]
  }

  wmat <- matrix(rep(w,ny),ny,length(w),byrow=T)
  wmat[which(is.na(ymat))] <- NA
  wmat <- wmat/apply(wmat,1,sum,na.rm=T)
  newy <- apply(ymat*wmat,1,sum,na.rm=T)
  newy

}

####################################################

distmat <- function(xt,yt,xs,ys){
    xd <- outer(xt,xs,function(xt,xs) (xt - xs)^2)
    yd <- outer(yt,ys,function(yt,ys) (yt - ys)^2)
    t(sqrt(xd + yd)) 
}
####################################################

tnorm <- function(n,lo,hi,mu,sig){   

  #normal truncated lo and hi

  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))

  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig)

  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == Inf]  <- lo[z == Inf]
  z[z == -Inf] <- hi[z == -Inf]
  z
}

nextTimeMat <- function(mat,last=rep(ncol(mat),nrow(mat)),INC=F){

  # mat is n by time matrix, returns same shifted to left

  til <- cbind(c(1:nrow(mat)),last)

  x <- cbind(mat[,-1],mat[,ncol(mat)])
  x[ til ] <- mat[ til ]

  if(INC){
    inc <- x[ cbind(c(1:nrow(mat)),last-1) ] - x[ cbind(c(1:nrow(mat)),last-2) ]
    x[ til ] <- x[ til ] + inc
  }

  x
}

lastTimeMat <- function(mat,first=rep(1,nrow(mat)),INC=F){  

  # mat is n by time matrix, returns same shifted to right
  # first repeats first value at the first time for obs i 

  tif <- cbind(c(1:nrow(mat)),first)
 
  x <- cbind(mat[,1],mat[,-ncol(mat)])
  x[ tif ] <- mat[ tif ]

  if(INC){   #increment first value
    inc <- x[ cbind(c(1:nrow(mat)),first+2) ] - x[ cbind(c(1:nrow(mat)),first+1) ]
    wf  <- which(first > 1)
    ff  <- x[ tif[wf,] ] - inc[wf]
    ff[ff < .1] <- .01
    x[ tif[wf,] ] <- ff
  }
    
  x
}


diffTimeMat <- function(mat,index,first = 0){  # increment matrix from obs by time matrix, pad last value

  md <- mat - lastTimeMat(mat,first) 
  vc <- md[index]
  list(dmat = md, vec = vc)
}


####################################################

tnorm.mvt <- function(avec,muvec,smat,lo=rep(-Inf,length(avec)),hi=rep(Inf,length(avec)),
                      whichSample=c(1:length(avec)),times=1){   

  # truncated multvariate normal
  # muvec is the vector of means
  # smat is the covariance matrix 
  # whichSample indicates which variables to sample

  if(length(lo) == 1)lo <- rep(lo,length(avec))
  if(length(hi) == 1)hi <- rep(hi,length(avec))

  for(j in 1:times){
   for(k in whichSample){

    tmp <- conditionalMVN(avec,muvec,smat,k)
    muk <- tmp$mu
    sgk <- tmp$vr

    if(length(muk) == 0)next

    avec[k] <- tnorm(1,lo[k],hi[k],muk,sqrt(sgk))
   }
  }
  avec
}


priorMVReg <- function(x,y,sigma){

  #Minka (2001) Bayesian Linear Regression

  # V - error covariance (sigma)
  # m - columns in x
  # d - columns in y
  #

  m <- ncol(x)
  d <- ncol(y)

  xx  <- crossprod(x)
  xy  <- crossprod(x,y)
  ixx <- invMat(xx)

  p1    <- invMat(sigma)%*%t(xy)%*%ixx%*%xy
  p2    <- m*d
  alpha <- p2/(sum(diag(p1)) + p2) 

  priorAcov <- kronecker(sigma,ixx/alpha)

  list(alpha = alpha, priorAcov = priorAcov)
}



#######################################################

dmvnormLog <- function(x,mu,sigma){    #multivariate normal 

  #mv normal density

    if (is.vector(x))x <- matrix(x, ncol = length(x))

    distval <- mahalanobis(x, mu, sigma)
    logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    logretval
}


#############################################################

rwish <- function(df,S){

  z  <- matrix(rnorm(df*nrow(S)),df,nrow(S))%*%chol(S)
  crossprod(z)
}

#####################################################
riwish <- function(v,S){

  solve(rwish(v,solve(S)))
}
#######################################

conditionalMVN <- function(x, mu, sigma, cindex){  #cindex is vector index for conditional

  tiny <- min(diag(sigma))*.0001
  nm   <- length(mu)
  if(length(x) != nm)stop('x and mu different length in conditionalMVN')

  x  <- matrix(x,nrow=1)
  mu <- matrix(mu,nrow=1)

  testv <- try(chol(sigma[-cindex,-cindex]),T)
  if(inherits(testv,'try-error')){
      return( list(mu = numeric(0), vr = numeric(0)) )
  }

  sin <- chol2inv(testv)
  p1  <- sigma[-cindex,cindex]%*%sin

  mu1 <- mu[cindex] + p1%*%(x[-cindex] - mu[-cindex])
  vr1 <- sigma[cindex,cindex] - p1%*%sigma[cindex,-cindex]
  if(vr1 < 0)vr1 <- tiny

  list(mu = mu1, vr = vr1)
}
##########################################################
processPars <- function(xgb,xtrue=numeric(0),CPLOT=F,DPLOT=F,
                        sigOnly = F,burnin=1,xlimits = NULL){  

  #xg      - matrix of gibbs chains
  #xtrue   - true values (simulated data)
  #CPLOT   - if T, plot chains
  #DPLOT   - if T, plot density
  #burnin  - analyze chains > burnin
  #xlimits - xlimits for plot
  #sigOnly - plot only parameters that 95% CI does not include 0
  
  if(!is.matrix(xgb))xgb <- matrix(xgb,ncol=1)
  if(is.null(colnames(xgb)))colnames(xgb) <- paste('V',c(1:ncol(xgb)),sep='-')
  
  if(sigOnly){
    wi   <- grep('intercept',colnames(xgb))      #extract covariates for plotting
    btmp <- xgb
    if(length(wi) > 0){
    	btmp <- xgb[,-wi]
      if(length(xtrue) > 0)xtrue <- xtrue[-wi]
    }

    wq   <- apply(btmp,2,quantile,c(.025,.975),na.rm=T)  #extract parameters != 0
    wq   <- which(wq[1,] < 0 & wq[2,] > 0)
    if(length(wq) > 0){
      xgb  <- btmp[,-wq]
      if(length(xtrue) > 0)xtrue <- xtrue[-wq]
    }
   }

  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  if(burnin > 1){
  	     if(burnin > (nrow(xgb) + 100))stop("burnin too large")
  	     xgb <- xgb[-c(1:burnin),]
  }
  if(!is.matrix(xgb))xgb <- as.matrix(xgb)
  nc <- ncol(xgb)
  nf <- round(sqrt(nc),0)

  out <- t(rbind(apply(xgb,2,mean,na.rm=T),apply(xgb,2,quantile,c(.025,.975),na.rm=T)))
  if(!is.null(colnames(xgb)))rownames(out) <- colnames(xgb)
  colnames(out) <- c('mean','0.025','0.975')
  if(length(xtrue) > 0){
    out <- cbind(out,xtrue)
    colnames(out) <- c('mean','0.025','0.975','true value')
  }

  armat <- matrix(0,nc,10)  #for AR model
  
 # for(j in 1:nc)armat[j,] <- ar(xgb[,j],aic=F,order.max = 10)$ar[1:10]

 # if(!is.null(colnames(xgb)))rownames(armat) <- colnames(xgb)
#  colnames(armat) <- paste('AR',c(1:10),sep='-')

  if(CPLOT | DPLOT)par(mfrow=c((nf+1),nf))
  if(CPLOT & DPLOT)par(mfrow=c((nf+1),nc))

  if(CPLOT){
      for(j in 1:nc){
       plot(xgb[,j],type='l')
       abline(h=out[j,],lty=2)
       if(length(xtrue) > 0)abline(h=xtrue[j],col='red')
       title(colnames(xgb)[j])
     }
  }
  xlims <- xlimits
  if(DPLOT){
      for(j in 1:nc){
        xj <- density(xgb[,j])
        if(is.null(xlimits))xlims <- range(xj$x)
        plot(xj$x,xj$y,type='l',xlim=xlims)
        abline(v=out[j,],lty=2)
        if(length(xtrue) > 0)abline(v=xtrue[j],col='red')
        title(colnames(xgb)[j])
     }
  }
  list(summary = signif(out,4)
)

}

#####################
distmat <- function(xt,yt,xs,ys){
    xd <- outer(xt,xs,function(xt,xs) (xt - xs)^2)
    yd <- outer(yt,ys,function(yt,ys) (yt - ys)^2)
    t(sqrt(xd + yd)) 
}
#################################
parSummary <- function(x){

    if(!is.matrix(x)){
      xb <- c(mean(x),sd(x),quantile(x,c(.025,.975)))
    }
    if(is.matrix(x)){
      xb <- cbind(apply(x,2,mean,na.rm=T),apply(x,2,sd,na.rm=T),
            t(apply(x,2,quantile,c(.025,.975),na.rm=T)))
    }
    xb
}
#############################################

y2zMVlogit <- function(y){     #multivar logit to fractions

  if(!is.matrix(y)){
     z <- inv.logit(y)
     return(cbind(z,1 - z))
  }

  zs   <- apply(exp(y),1,sum)
  z1   <- 1/(1 + zs)
  zm   <- exp(y)/ (1 + zs)
  cbind(zm,z1)
}

z2yMVlogit <- function(z){     #fractions to multivar logit

  r <- ncol(z)
  if(r == 2){
    ss <- z[,1]/apply(z,1,sum)
    return(log(ss/(1 - ss) ))
  }

  log(z[,-r]/(1 - apply(z[,-r],1,sum)))
  
}
###########################################################

dmvnormZeroMean <- function(x,sigma){          #MVN density for mean 0

    testv <- try(chol(sigma),T)
    if(inherits(testv,'try-error')){
 #      message(paste('error in mydmvnorm at g = ',g,sep=''))
       tiny <- min(abs(x))/100 + 1e-5
       sigma <- sigma + diag(diag(sigma + tiny))
       testv <- try(chol(sigma),T)
    }

    covx <- chol2inv(testv)
    distval <- rowSums((x %*% covx) * x)
  
    ev <- eigen(sigma, only.values = TRUE)$values 
    if(min(ev) < 0)ev <- nearPD(sigma)$eigenvalues

    logdet   <- sum(log( ev ))
    -(ncol(x) * log(2 * pi) + logdet + distval)/2
}
########################################

getChainMeanVar <- function(chains){  #mean vector, covariance matrix from chains matrix

  chains <- row2Mat(chains)
  nn     <- nrow(chains)
  if(nn == 1){
    mub    <- chains
    vrb    <- diag(1,ncol(chains))
  }
  if(nn > 1){
     mub   <- apply(chains,2,mean)
     vrb   <- cov(chains)
  }
  list(mu = mub, vr = vrb)
}
#####################################
getPost <- function(b,bchain){     #posterior estimate for marginal likelihood

  #gaussian
  # method of Chib 1995, JASA

 # b <- row2mat(b)

  tmp <- getChainMeanVar(bchain)
  mub <- tmp$mu
  vrb <- tmp$vr
  
  mu  <- matrix(mub,nrow(bchain),ncol(bchain),byrow=T) - bchain
 
  log( mean(exp( dmvnormZeroMean(mu,vrb) )) )

}

#######################################3
biVarMoments <- function(x1,x2,wt,PLOT = F, LINES=F, color=1){  #x1, x2 - vectors for variables, wt is weight

  if(length(wt) == 1)wt <- rep(wt,length(x1))

  ww <- which(is.finite(x1) & is.finite(x2))

  x1 <- x1[ww]
  x2 <- x2[ww]
  wt <- wt[ww]

  w1 <- x1*wt
  w2 <- x2*wt
  m1 <- sum(w1)/sum(wt)
  m2 <- sum(w2)/sum(wt)

  v1  <- sum(wt*x1^2)/sum(wt) - m1^2
  v2  <- sum(wt*x2^2)/sum(wt) - m2^2
  c   <- sum(wt*(x1 - m1)*(x2 - m2))/sum(wt)

  if(PLOT | LINES){
    tmp <- list(loc = c(m1,m2), cov = matrix(c(v1,c,c,v2),2,2), d2 = qchisq(.95,1) )
    tmp <- predict.ellipsoid(tmp)
    if(PLOT)plot(tmp[,1],tmp[,2],type='l',col=color)
    if(LINES)lines(tmp[,1],tmp[,2],type='l',col=color)
 #   tmp <- list(loc = c(m1,m2), cov = matrix(c(v1,c,c,v2),2,2), d2 = qchisq(.99,1) )
 #   tmp <- predict.ellipsoid(tmp)
 #   lines(tmp[,1],tmp[,2],type='l',col=color)
  }

    list(mu = c(m1,m2), var = matrix(c(v1,c,c,v2),2,2) )
}
#############################################################

rtrunc <- function(n,lo,hi,p1,p2,FN){    #truncated using inv dist sampling
   #truncated 2 parameter distribution
	
  if(FN == 'invGamma'){
  	z1 <- 1 - pgamma(1/lo,p1,p2)
  	z2 <- 1 - pgamma(1/hi,p1,p2)
   z  <- 1 - runif(n,z1,z2)
   return(1/qgamma(z,p1,p2))
  }

  pf <- paste('p',FN,sep='')
  qf <- paste('q',FN,sep='')
  pf1 <- match.fun(pf)
  qf1 <- match.fun(qf)
  
  z1 <- pf1(lo,p1,p2)
  z2 <- pf1(hi,p1,p2)
  z  <- runif(n,z1,z2)
  
  qf1(z,p1,p2)
}

acceptMH <- function(p0,p1,x0,x1,BLOCK=F){   #accept for M, M-H
	# if BLOCK, then accept as a block,
	# otherwise, accept individually

  nz          <- length(x0)  #no. to accept
  if(BLOCK)nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a,arr.ind=T)
  
  if(BLOCK & length(keep) > 0)x0 <- x1
  if(!BLOCK)                  x0[keep] <- x1[keep]           
  ac <- length(keep)        

  list(x = x0, accept = ac)
}

predVsObs <- function(o,p){ 
	
  #o  - length n vector of obs or true values
  #p - ng by n matrix of estimates
  
  n <- length(o)
  y <- apply(p,2,quantile,c(.5,.025,.975))

  plot(o,y[1,])
  for(j in 1:n)lines(c(o[j],o[j]),y[2:3,j])
  abline(0,1,lty=2)
  y
}
#######################################################
plotObsPred <- function(obs,yMean,ySE=NULL,nbin=NULL,nPerBin=NULL,log=F,ylimit=NULL,
                        xlabel='Observed',ylabel='Predicted'){

   ptcol <- 'black'
   ptcol <- 'wheat'
   if(!is.null(nbin))ptcol <- 'grey'

   if(is.null(ylimit)){
     if(!log)plot(obs,yMean,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel)
     if(log) plot(obs,yMean,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,log='xy')
   }
   if(!is.null(ylimit)){
     if(!log)plot(obs,yMean,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,ylim=ylimit)
     if(log) plot(obs,yMean,col=ptcol,cex=.3,xlab=xlabel,ylab=ylabel,log='xy',ylim=ylimit)
   }
   if(!is.null(ySE)){
     ylo <- yMean - 1.96*ySE
     yhi <- yMean + 1.96*ySE
     for(i in 1:length(obs))lines(c(obs[i],obs[i]),c(ylo[i],yhi[i]),col='grey',lwd=2)
   }

   if(!is.null(nbin) | !is.null(nPerBin)){

     if(is.null(nbin))nbin <- 20
     bins <- seq(min(obs,na.rm=T),max(obs,na.rm=T),length=nbin)

     if(!is.null(nPerBin)){
         nbb <- nPerBin/length(obs)
         nbb <- seq(0,1,by=nbb)
         if(max(nbb) < 1)nbb <- c(nbb,1)
         bins <- quantile(obs,nbb,na.rm=T)
         nbin <- length(bins)
     }

     wide <- diff(bins)/2
     for(k in 1:(nbin-1)){
        q  <- quantile(yMean[obs >= bins[k] & obs <= bins[k+1]],c(.5,.025,.158,.841,.975),na.rm=T)
        xx <- mean(bins[k:(k+1)])
        lines(c(xx,xx),q[c(2,5)],lwd=3)
        lines(c(xx-.6*wide[k],xx+.6*wide[k]),q[c(1,1)],lwd=3)
        rect(xx-.4*wide[k],q[3],xx+.4*wide[k],q[4],col='grey')
      }
   }

 #  points(obs,yMean,col=ptcol,cex=.3)
}
#######################################################
deviance <- function(y,x,b,s=0,LIKE){
	    
    if(LIKE == 'norm')  dv <- dnorm(y,x%*%b,sqrt(s),log=T) 
    if(LIKE == 'pois')  dv <- dpois(y,exp(x%*%b),log=T)
    if(LIKE == 'binom') dv <- dbinom(y,1,invlogit(x%*%b),log=T)
    if(LIKE == 'mvnorm')dv <- dmvnorm(y,x%*%b,s,log=T)
    if(LIKE == 'multinom')dv <- multinomLike(y,x,b)
    -2*dv
}
#####################################33
dinvGamma <- function(x,a,b,log=FALSE){
	
	p <- a*log(b) - lgamma(a) - (a+1)*log(x) - b/x
	if(log)return(p)
	exp(p)
}
###########################################
dbivarNormFromCols <- function(y,mu,S){    #bivariate norm, covariances for each y in columns, log density

  #y  - n by 2
  #mu - n by 2 or scalar 0
  #S  - n by 4: S11, S12, S21, S22

  n <- nrow(y)
  if(length(mu) == 0)mu <- y*0

  y <- y - mu

  invS <- invertcol2(S[,1],S[,2],S[,3],S[,4])
  z    <- y[,1]^2*(invS[,1] + invS[,2])+ y[,2]^2*(invS[,3] + invS[,4])

  ldet <- log(S[,1]*S[,4] - S[,2]*S[,3])

  -(n*log(2*pi) + ldet + z)/2
}

#############################################################

invertcol2 <- function(S1,S12,S21,S2){   
    # inverse of matrix supplied and returned as vectors: var1, cov12, cov12, var2

  dt  <- S1*S2 - S12*S21
  I1  <- S2
  I2  <- S1
  I12 <- -S12
  I21 <- -S21

  cbind(I1,I12,I21,I2)/dt
}
############################################################

initialStatesSS <- function(likelihood,y,priorObsEr,priorObsErWt,bias=1){

  library(stats)
  if(!is.matrix(y))y <- matrix(y,ncol=1)
  r <- ncol(y)

  n    <- length(y)
  time <- c(1:n)
  wm   <- which(is.na(y))
  notMiss <- c(1:n)
  
  x <- y*bias

  if(length(wm) > 0){
   notMiss <- notMiss[-wm]
   x[wm]   <- predict(smooth.spline(time[-wm],y[-wm]),time)$y[wm]
  }

  propSd <- sd(diff(as.vector(y)),na.rm=T)
  
  if(likelihood == 'dpois')x <- log(y + .1)
  
  s1s2 <- c(0,0)
  if(likelihood == 'dnorm'){
   	 s1 <- priorObsErWt
   	 s2 <- priorObsEr*(s1 - 1)
   	 s1s2 <- c(s1,s2)
  }

  list(x = as.vector(x), propSd = propSd,notMiss = notMiss, miss = wm, s1s2 = s1s2)
}


############################################################
condProb <- function(mu,V,y){     #bivariate conditional on y for 2nd variables

  cm <- mu[1] + V[1,2]/V[2,2]*(y - mu[2])
  cV <- V[1,1] - V[1,2]/V[2,2]*V[2,1]

  list(muY = cm, cV = cV)
}
###########################################

multiLogitStates <- function(b,discStates,contStates,maxD){    # continuous based on discrete states

      tiny <- 1e-20
      nn   <- length(discStates)
 
      discStates[is.na(discStates)] <- 0
      prob <- discStates*0 + 1

      if(maxD == 1)return(log(prob))

      wk <- which(discStates == 1 & is.finite(contStates),arr.ind=T)

      prob[wk] <- invlogt(b[1,],contStates[wk])

      if(maxD > 2){
        for(k in 2:(maxD-1)){
         wk <- which(discStates == k & is.finite(contStates),arr.ind=T)
         prob[wk] <- invlogt(b[k,],contStates[wk]) - 
                                  invlogt(b[(k-1),],contStates[wk])
        }
      }
      prob[discStates == maxD] <- 1 - invlogt(b[(maxD-1),],contStates[discStates == maxD])

    prob[prob < tiny] <- tiny
    log(prob)
}

############################################
breaks2pars <- function(cmat,breaks){        #break points and coefficents to ordinal multinomial logit pars

  nd <- length(breaks)
  c0 <- rep(0,nrow(cmat))

  for(k in 1:nd){

    qk <- 0
    if(k > 1){
      for(j in 1:(k-1)){
       ej <- exp(c0[j] + cmat[j,2]*breaks[k])
       qk <- qk + ej/(1 + ej)
      } 
    }
    D     <- -log(1/(.5 + qk) - 1)
    cc0   <- D - cmat[k,2]*breaks[k]
    if(k > 1){
      if(cc0 < cmat[k-1,1]){
        cc0 <- cmat[k-1,1] + .5
        cmat[k,2] <- (D - cc0)/breaks[k]
      }
    }
    cmat[k,1] <- cc0
  }
  cmat
}
#################################################
proposeBreak <- function(cmat,br,brLims){      #propose new breakpoints given limits brLims

  c1 <- cmat[,1]
  c2 <- cmat[,2]

  nd <- nrow(cmat)
  brNew <- tnorm(nd,brLims[1,],brLims[2,],br,.1)

  if(nd == 1){
     c1New <- tnorm(1,1,100,cmat[1],.1)
     c2New <- (log(.5) + log(1 + exp(c1New + cmat[2]*brNew)) - c1New)/brNew
     cnew <- matrix(c(c1New,c2New),1,2)
  }

  if(nd > 1)cnew <- getCmat(cmat,brNew)

  list(cmat = cmat,br = brNew)
}
##########################################
pars2p <- function(cmat,h){               #multinomial logit pars and scale h to Pr for multinom logit

  nd <- nrow(cmat)
  nh <- length(h)

  c1 <- matrix(cmat[,1],nh,nd,byrow=T)
  c2 <- matrix(cmat[,2],nh,nd,byrow=T)
  hh <- matrix(h,nh,nd)
  eh <- exp(c1 + c2*hh)

  theta <- matrix(0,nh,nd+1)
  sumt  <- rep(0,nh)

  for(k in 1:nd){

     tk   <- eh[,k]/(1 + eh[,k]) - sumt
     theta[,k] <- tk 
     sumt <- sumt + tk
  }
  theta[,nd+1] <- 1 - apply(theta,1,sum)

  theta
}


######################################################
plotLogit <- function(lims,breaks,h,cmat){  #plot multinomial ordinal logit

  tmp  <- pars2p(cmat,h)

  plot(h,tmp[,1],type='l',lwd=2,xlab='latent health scale',ylab='Probabilty')

  for(j in 1:ncol(lims)){
 #   polygon(c(lims[,j],rev(lims[,j])),c(0,0,1,1),border=NA,col='grey')
    lines(h,tmp[,j],col=j,lwd=2)
    l1   <- 0
    if(j > 1)l1 <- breaks[j-1]
    midx <- (l1 + breaks[j])/2
    text(midx,.8,j,col=j,cex=1.4)
  }
  lines(h,tmp[,j+1],col=j+1,lwd=2)
  midx <- ( breaks[j] + 100 )/2
  text(midx,.8,j+1,col=j+1,cex=1.4)
  abline(h=.5,lty=2)
}   	
###########################

plotSetup <- function(xtic,ytic,xvals = xtic, yvals = ytic, xlabel=' ',ylabel=' ',endFactor=c(.05,.05),
                      fc = 'azure2',lc='white'){

  xr <- range(xtic)
  yr <- range(ytic)

  xlimit <- c(xr[1] - diff(xr)*endFactor[1],xr[2] + diff(xr)*endFactor[1])
  ylimit <- c(yr[1] - diff(yr)*endFactor[2],yr[2] + diff(yr)*endFactor[2])

  plot(-100,-100,xlim=xlimit,ylim=ylimit,ylab=ylabel,cex=1.2,xlab=xlabel,xaxt='n',yaxt='n')

  rect(xlimit[1],ylimit[1],xlimit[2],ylimit[2],col=fc,border=NA)

  axis(1,at=xtic,labels=xvals,cex.lab=1.2)
  axis(2,at=ytic,labels=yvals,cex.lab=1.2)

  abline(v=xtic,lwd=3,col=lc)
  abline(h=ytic,lwd=3,col=lc)
}



getQuant <- function(data,dim,q){

  quantile(data,dim,q,na.rm=T)
}


####################################################

myBoxPlot <- function(mu,...,ORDER = F,xvalues=c(1:length(mu)),boxcol='black',
                      xtic=c(1,length(mu)),ytic=NULL,plabels=NULL,shadelo=NULL,shadehi=NULL,add=F){

  #mu has median
  #... has lo to hi as paired columns, inner to outer
  #shadehi - shade above this value
  #shadelo - shade below this value

  n  <- length(mu)
  xi <- xvalues
  or <- c(1:n)
  if(ORDER)or <- order(mu,decreasing=T)

  x  <- list(...)
  nq <- length(x)

  if(is.null(ytic))ytic <- signif(range(mu),1)

  if(!add){
    plotSetup(xtic,ytic)
  }

  if(!is.null(shadehi))rect(xtic[1],shadehi,max(xtic),max(ytic),col='bisque1',border=NA)
  if(!is.null(shadelo))rect(xtic[1],min(ytic),max(xtic),shadelo,col='bisque1',border=NA)

  points(xi,mu[or],pch=3,col=boxcol)

 # ji <- c(1:nq)/diff(range(xi))*3

  ji <- .2*diff(xi)
  ji <- c(ji[1],ji)

  jk <- seq(.2,.8,length=nq)

  for(i in 1:n){
    maxi <- ytic[1]
    for(k in 1:nq){
      yi <- x[[k]][or[i],]
      if(is.na(yi[1]))next
      
      ki <- jk[nq - k + 1]*ji[i]

      x1 <- xi[i] - ki
      x2 <- xi[i] + ki

      rect(x1,yi[1],x2,yi[2],lwd=2,col=boxcol,border=boxcol)
      if(k == nq)lines(c(xi[i],xi[i]),yi,lwd=2,col=boxcol)
      if(max(yi) > maxi)maxi <- max(yi)
    }
    if(maxi > max(ytic))maxi <- max(ytic) - .1*diff(range(ytic))
    if(!is.null(plabels))text(i,maxi,plabels[or[i]],srt=90,pos=4,col=boxcol)
  }

}

trim <- function(xx,p){     #p - lower,upper %tile

  xq <- quantile(xx,p,na.rm=T)
  
  xx[xx < xq[1]] <- xq[1]
  xx[xx > xq[2]] <- xq[2]

  xx
}

#################################################
weightAve <- function(x,wt){  #x vector for variable, wt is weight

  w1 <- x*wt
  m1 <- sum(w1)/sum(wt)

  v1  <- sum(wt*x^2)/sum(wt) - m1^2

  c(m1,v1)
}

cov2cor <- function(covmat){  #covariance matrix to correlation matrix

  d    <- nrow(covmat)
  di   <- diag(covmat)
  covmat/ sqrt(matrix(di,d,d)) / sqrt(matrix(di,d,d,byrow=T))
}
#####################################################
corPlot <- function(cmat,slim=NULL,diag0=F,textSize=1){  #correlation or covariance matrix

  if(diag0)diag(cmat) <- 0

  cmat[lower.tri(cmat)] <- 0

  d <- nrow(cmat)
  xtext <- rep(c(1,100),d/2)
  if(length(xtext) < d)xtext <- c(xtext,1)

  if(d < 20)xtext <- xtext*0 + 1

  xtext <- xtext*0 + 1

  ncol <- 100
  colseq <- mapColors(ncol)


  if(is.null(slim))slim = range(cmat)
  slim  <- signif(slim,1)
  scale <- seq(slim[1],slim[2],length.out=ncol)

  ww   <- as.matrix(expand.grid(c(1:d),c(1:d)))

  ww  <- ww[ww[,1] <= ww[,2],]
  ww  <- ww[order(ww[,1]),]

  icol <- findInterval(cmat[ww],scale,all.inside=T)
  coli <- colseq[icol]

  symbols(ww[,1],ww[,2]+1,squares=rep(.9,nrow(ww)),xlim=c(0,d+2),ylim=c(0,d+2),
              fg=coli,bg=coli,inches=F,axes=F,xlab=' ',ylab=' ')

  for(kk in 1:(d+1)){
        ks <- kk - .5
        if(kk <= d)lines(c(ks,ks),c(ks+1,d+1.5),col='grey')
        if(kk > 1) lines(c(.5,ks-1),c(ks,ks),col='grey')
        if(kk <= d)lines(c(ks,ks+1),c(ks+1,ks),col='grey')
  }
  lines(c(.5,ks-1),c(ks+1,ks+1),col='grey')
  text(c(1:d)+.1*xtext,c(1:d)+.4,colnames(cmat),pos=4,srt=-45,cex=textSize)

  colorLegend(c(d+1,d+2),c(0,d/2.2),ytick=c(slim[1],0,slim[2]),scale,colseq,labside='left')
}
##################################################

zeroOneScale <- function(x) (x - min(x,na.rm=T))/(max(x,na.rm=T) - min(x,na.rm=T))

######################################################

mapContoursUSA <- function(xx,yy,z,zlevs=NULL,yname=NULL){ #must have lon,lat in main program

  map('state',interior=F,col='grey',xlim=c(-97,-65))
  map('state',boundary=F,col='grey',add=T)

  values2contour(xx,yy,z,q=1,zlevs,lwd=1,col='black',add=T)
  if(!is.null(yname))title(yname)
 
}
######################################################

arrowField <- function(xy0,xy1=NULL,directionLong=NULL,angle=20,col=1){

  #xy0, xy1, 2 columns each
  #directionLong in radians, length

  if(is.null(xy1) & is.null(directionLong))stop('either xy1 or directionLong must be given')

  if(!is.null(directionLong)){
     xy1     <- xy0*0
     xy1[,1] <- xy0[,1] + directionLong[,2]*cos(directionLong[,1])
     xy1[,2] <- xy0[,2] + directionLong[,2]*sin(directionLong[,1])
     long    <- sqrt( (xy1[,1] - xy0[,1])^2 + (xy1[,2] - xy0[,2])^2)
  }

  lwide <- 2*zeroOneScale(long)

  arrows(xy1[,1],xy1[,2],xy0[,1],xy0[,2],length=.05*long,angle,col,lwd=lwide)
}

#################################################################3
mapArrowsUSA <- function(xx,yy,direction,long,yname=NULL){

  require(spatial)

  map('state',interior=F,col='grey',xlim=c(-97,-65))
  map('state',boundary=F,col='grey',add=T)
      

    wf <- which(is.finite(direction) & is.finite(long))

    lo <- xx[wf]
    la <- yy[wf]

    arrowField(cbind(lo,la),directionLong=cbind(direction,long)[wf,],angle=20,col=1)
 
  if(!is.null(yname))title(yname)

  map('state',interior=F,add=T,lwd=4,col='white')
  map('state',interior=F,add=T)

 # map('state',boundary=F,col='white',lwd=3,add=T)
  map('state',boundary=F,col='grey',add=T)
 # title(mapname)

}
#######################################3


mapPointsUSA <- function(xx,yy,...,sym='values',cols='black',
                         fill='none',yname=NULL,scale=NULL,KRIG=F){  

  require(spatial)

  map('state',interior=F,col='grey',xlim=c(-97,-65))
  map('state',boundary=F,col='grey',add=T)

  zz  <- list(...)

  if(length(zz) > 1){
    if(length(sym) == 1) sym  <- rep(sym,length(x))
    if(length(cols) == 1)cols <- rep(cols,length(x))
    if(length(fill) == 1)fill <- rep(fill,length(x))
  }

  if('ramp' %in% cols | 'ramp' %in% fill){

    ncol   <- 100
    colseq <- mapColors(ncol)

  }

  for(j in 1:length(zz)){

    xj <- zj <- zz[[j]]
    wf <- which(is.finite(xj))
    xj <- xj[wf]

    lo <- xx[wf]
    la <- yy[wf]

    if(sym[j] == 'ones'){
       w0 <- which(xj == 1)
       xj <- rep(.1,length(w0))
       ll <- la[w0]
       lo <- lo[w0]
    }

    fj <- cols[j]
    bj <- fill[j]
    if(bj == 'none')bj <- NA

    if(cols[j] == 'ramp' | fill[j] == 'ramp'){
      fj <- bj <- colseq[findInterval(xj,seq(min(xj,na.rm=T),max(xj,na.rm=T),length=ncol))]
    }

    if(KRIG){
        tmp <- surf.gls(2,expcov,x=lo,y=la,xj,d=.7)
        ps  <- prmat(tmp,xl=min(lo),xu=max(lo),yl=min(la),yu=max(la),n=50)
        xy  <- expand.grid(ps$x,ps$y)
        fj <- bj <- colseq[findInterval(ps$z,seq(min(ps$z,na.rm=T),max(ps$z,na.rm=T),length=ncol))]
        symbols(xy[,1],xy[,2],squares=rep(.5,nrow(xy)),inches=F,fg=fj,bg=bj,add=T)
        zj <- ps$z
          
    }

    if(!KRIG){
      if(is.null(scale)) xj <- zz[[j]]
      if(!is.null(scale))xj <- rep(scale,length(zz[[j]]))
      symbols(xx,yy,circles=xj,inches=F,fg=fj,bg=bj,add=T)
      zj <- xj
    }

    tscale <- signif(range(zj,na.rm=T),1)
    ttic   <- tscale
    if(tscale[1] < 0 & tscale[2] > 0)ttic   <- c(tscale[1],0,tscale[2])
    tscale <- seq(tscale[1],tscale[2],length=ncol)

    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    colorLegend(c(-70,-69),c(26,35),ytick=ttic,scale=tscale,colseq,labside='left')
  }

  if(!is.null(yname))title(yname)

  map('state',interior=F,add=T,lwd=4,col='white')
  map('state',interior=F,add=T)

 # map('state',boundary=F,col='white',lwd=3,add=T)
  map('state',boundary=F,col='grey',add=T)
 # title(mapname)

   

}
############################################3
points2angle <- function(xy0,xy){             #(x,y) columns in xy from reference (x,y) in xy0, returns angle

  h  <- distmat(xy0[1],xy0[2],xy[,1],xy[,2])
  dx <- xy[,1] - xy0[1]
  dy <- xy[,2] - xy0[2]

  list(theta = atan2(dy,dx), dist = h)
}
#################################################
localSlope <- function(xk,yk,zk,direction=NULL){ #if not null, direction is theta

  x0 <- mean(xk)
  y0 <- mean(yk)
  z0 <- mean(zk)

  xk <- xk - x0
  yk <- yk - y0
  zk <- zk - z0

  u1 <- u2 <- 1

  X <- cbind(xk,yk)

  if(length(xk) == 2){
    fx <- diff(zk)/diff(xk)
    fy <- diff(zk)/diff(yk)
  }

  if(length(xk) > 3){
    b <- invMat(crossprod(X),NEARPD=T)%*%crossprod(X,zk)
    fx <- b[1,]
    fy <- b[2,]
  }

  if(is.null(direction)){     #maximum derivative
    grade <- sqrt(fx^2 + fy^2)
    theta <- atan2(fy,fx)
  }

  if(!is.null(direction)){    #derivative in direction theta
     theta <- direction
     u1    <- cos(direction)
     u2    <- sin(direction)
     grade <- fx*u1 + fy*u2
  }

  list(theta = theta, grade = grade)
}


mapColors <- function(ncol){

    colF   <- colorRampPalette(c('darkblue','blue','lightblue','green','lightgreen','yellow','orange','red','brown'))
    colF(ncol)
}




#############################################

climGradient <- function(xx,yy,z1,z2,dthreshold=1,mapname=NULL,climGrad=F,climDir=F,vegDir=F,ABS=T,DYDX=T,DTHETA=F){  

# directional gradient in z2 relative to maximum gradient in z1
# climGrad - map of climate gradient derivative
# climDir  - map of climate direction
# vegDir   - vegetation, rather than climate direction
# abs value- when there are multiple species

  require(spatial)

  ncol <- 100
  cols <- mapColors(ncol)

  if(!is.matrix(z2))z2 <- as.matrix(z2)
  r   <- ncol(z2)

  out <- matrix(NA,length(xx),5)
  colnames(out) <- c('theta','grad','thetaV','dtheta','dydx')

  distance <- distmat(xx,yy,xx,yy)

  for(i in 1:length(xx)){

    wi <- unique(which(distance[i,] < dthreshold))
    if(length(wi) < 4)wi <- order(distance[i,])[1:4]
    tmp <- localSlope(xx[wi],yy[wi],z1[wi])
    theta <- tmp$theta
    dx    <- tmp$grade

    wz <- which(z2[wi,] == 0)  #zeros do not occur at this location
    cdir <- theta

    gvec <- vslp <- rep(NA,r)

    for(j in 1:r){
      zj  <- z2[wi,j]
   #   wj  <- which(zj != 0)
   #   if(length(wj) < 4)next
      if(vegDir)cdir <- NULL
   #   tmp <- localSlope(xx[wi[wj]],yy[wi[wj]],zj[wj],direction=cdir)

      vslp[j] <- localSlope(xx[wi],yy[wi],zj)$theta
      tmp  <- localSlope(xx[wi],yy[wi],zj,direction=cdir)
      gvec[j] <- tmp$grade
    }

    dydx <- gvec/dx
    if(ABS)dydx <- abs(dydx)
    vslp    <- mean(vslp,na.rm=T)
    dtheta  <- atan2(sin(theta - vslp),cos(theta - vslp))
    out[i,] <- c(theta,dx,vslp,dtheta,mean(dydx,na.rm=T))
  }

  qGrad <- c(.01,.99)
  qdydx <- c(.01,.99)

 # tmap <- cos( out[,'theta'] )
  tmap <- out[,'theta']
  
  gmap <- zeroOneScale( trim( out[,'grad'],qGrad) )
  
  dtmp <- trim( out[,'dydx'],qdydx)
 # dtmp <- log10(dtmp)
  dtmp[dtmp == -Inf] <- min(dtmp[dtmp != -Inf],na.rm=T)

  dmap <- zeroOneScale( dtmp )
  dmap[out[,'dydx'] == 0] <- min(dmap,na.rm=T)

  mmap <- zeroOneScale( abs(out[,'dtheta']) )


 # dmap <- zeroOneScale( trim( out[,'dydx'],qdydx) )

  mscale <- signif( c(0,pi),2)
  tscale <- signif(c(-pi,pi),2)
  gscale <- signif(quantile(out[,'grad'],qGrad,na.rm=T),1)
  dscale <- signif(quantile(out[,'dydx'],qdydx,na.rm=T),1)

  ttic   <- c(tscale[1],0,tscale[2])
  gtic   <- gscale
  dtic   <- dscale
  mtic   <- mscale

  if(gtic[1] < 0 & gtic[2] > 0)gtic <- c(gtic[1],0,gtic[2])
  if(dtic[1] < 0 & dtic[2] > 0)dtic <- c(dtic[1],0,dtic[2])

  tscale <- seq(tscale[1],tscale[2],length=ncol)
  gscale <- seq(gscale[1],gscale[2],length=ncol)
  dscale <- seq(dscale[1],dscale[2],length=ncol)
  mscale <- seq(mscale[1],mscale[2],length=ncol)

  if(climDir){
    mapPointsUSA(xx,yy,tmap,sym='values',cols='ramp',
                         fill='none',yname=paste(mapname,'theta',sep=' '),KRIG=T)
    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)

    colorLegend(c(-70,-69),c(26,35),ytick=ttic,scale=tscale,cols,labside='left')
  }
  if(climGrad){
    mapPointsUSA(xx,yy,gmap,sym='values',cols='ramp',
                         fill='none',yname=paste(mapname,'Gradient',sep=' '),KRIG=T)
    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    colorLegend(c(-70,-69),c(26,35),ytick=gtic,scale=gscale,cols,labside='left')
  }

  if(DYDX){
    mapnew <- paste(mapname,'dxdy',sep=' ')
 
    mapPointsUSA(xx,yy,dmap,sym='values',cols='ramp',
                         fill='none',yname=mapnew,KRIG=T)
    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    colorLegend(c(-70,-69),c(26,35),ytick=dtic,scale=dscale,cols,labside='left')
  }
  if(DTHETA){
    mapnew <- paste(mapname,'del theta',sep=' ')
 
    mapPointsUSA(xx,yy,-(mmap-1),sym='values',cols='ramp',
                         fill='none',yname=mapnew,KRIG=T)
    polygon(c(-75,-65,-65,-75),c(25,25,36,36),col='white',border=NA)
    colorLegend(c(-70,-69),c(26,35),ytick=rev(-mtic),scale=mscale,cols,labside='left')
  }

  invisible(out)
}

################################################3333

associationPlot <- function(z,groupIndex){       

  #z - n rows by c columns give abundance for each species
  #groupIndex - length n

  allnames <- sort(unique(groupIndex))
  ns       <- length(allnames)
  nf       <- ncol(z)

  ncol <- 50
  scale <- seq(-1,1,length.out=ncol)
  
  freqY <- numeric(0)

  nr <- round(ns/2,0) + 1
  par(mfrow=c(nr,2),bty='n',mar=c(1,1,3,1))

  for(i in 1:ns){

      yi <- z[groupIndex == allnames[i],]

      yc <- cor(yi)

      corPlot(yc,slim=c(-1,1),diag0=T)
      
      title(hostNames[i])
 
  }

  dev.print(device=postscript,file='associationPlot.ps',width=6,horizontal=F)
}
############################################################
colorLegend <- function(x,y,ytick=NULL,scale,cols,labside='right'){  #x = (x1,x2), y = (y1,y2)

#  rect(x[1],y[1],x[2],y[2])

  nn <- length(scale)
  ys <- seq(y[1],y[2],length=nn)

  for(j in 1:(length(scale)-1)){
    rect(x[1],ys[j],x[2],ys[j+1],col=cols[j],border=NA)
  }
  if(!is.null(ytick)){

    dx <- diff(x)
    yy <- ys[findInterval(ytick,scale)]
    lx <- c(x[2],x[2]+dx/3)
    tx <- x[2]+dx/2.5
    tp <- 4
    if(labside == 'left'){
         lx <- c(x[1]-dx/3,x[1])
         tx <- x[1]-dx/2.5
         tp <- 2
    }
    for(j in 1:length(ytick)){
       lines(lx,c(yy[j],yy[j]))
       text(tx,yy[j],ytick[j],pos=tp)
    }
  }
}

#################################################################
gibbsLoop <- function(LIKE,ng,x,y,b,startpred=0,sigma = 0,z = numeric(0)){

  kx     <- length(b)
  
  bgibbs <- matrix(NA,ng,kx)
  colnames(bgibbs) <- paste('b',c(1:kx),sep='-')
  sgibbs <- matrix(NA,ng,max(1,length(sigma)))
  
  r <- 1                       #responses
  if(is.matrix(y))r <- ncol(y)
  
  if(sigma[1] == 0)sgibbs <- numeric(0) #variances
  
  sg <- sigma
  
  if(is.matrix(sigma)){                  #Wishart prior
    sg        <- prior.W
    colnames(sgibbs) <- outer(rownames(sigma),rownames(sigma),paste,sep='_') 
    colnames(bgibbs) <- as.vector(outer(colnames(x),colnames(b),paste,sep='_'))
  }
  
  pBVar <- solve(crossprod(x))

  pred <- pred2 <- rep(0,nrow(x)*r)

  bg <- b
  yg <- y
  
  y1 <- yg
  
  if(LIKE == 'pois') bg <- pBVar%*%crossprod(x,log(y + .1))
  if(LIKE == 'binom'){
  	yl <- y
  	yl[yl == 0] <- .1
  	yl[yl == 1] <- .9
  	bg <- pBVar%*%crossprod(x,logit(yl))
  }
  
  if(LIKE == "mvnorm-multinom"){
  	y1   <- prob2Logit(yg)
  }
  if(LIKE == 'multinom')pBVar <- diag(.01,k*(r-1))
  
  dev <- 0
  np  <- 0

  for(g in 1:ng){

print(g)
    bg <- bUpdateGibbs(x,y1,bg,LIKE,sigma=sg,pBVar)
    
    if(sigma[1] > 0 & length(sigma) == 1){
    	 sg <- sigmaUpdate(x,y1,bg)
    }
    
    if(length(grep('mvnorm',LIKE)) > 0){
    	 sinv <- rwish(prior.WDF + nrow(y1),y1 - x%*%bg)
    	 sg   <- solve(sinv)
    }
    
    if(LIKE == "mvnorm-multinom"){
    	y1 <- ysampMvnormMultinom(x,y1,z,bg,sg)$y
    	yg <- logit2Prob(y1)
    }
    
    dev <- dev + sum(deviance(y1,x,bg,sg,LIKE))

    if(g > startpred){
      py    <- as.vector(simY(x,bg,LIKE,r,size,sg))
      pred  <- pred + py
      pred2 <- pred2 + py^2
      np    <- np + 1
    }

    bgibbs[g,] <- bg
    if(length(sgibbs) > 0)sgibbs[g,] <- sg
    
    if(g %in% c(100,200,500,1000)){
    	pBVar <- .1*cov(bgibbs[20:g,])
    }
  }

  ymean <- pred/np
  yse   <- sqrt(pred2/np - ymean^2)
  
  if(r > 1){
  	ymean <- matrix(ymean,n,r)
  	yse   <- matrix(yse,n,r)
  }
  bmean <- apply(bgibbs,2,mean)
  bmean <- matrix(bmean,nrow(bg),ncol(bg))
  smean <- numeric(0)
  if(length(sg) > 1) smean <- matrix(apply(sgibbs,2,mean),nrow(sg),ncol(sg))
  if(length(sg) == 1)smean <- mean(sgibbs)
  
  meanDev <- dev/ng
  
  pd  <- meanDev - sum(deviance(y1,x,bmean,smean,LIKE))
  dic <- 2*pd + meanDev

  list(bgibbs = bgibbs,sgibbs = sgibbs, ymean = ymean, yse = yse, dic = dic)
}
####################################################
simY <- function(x,b,LIKE,r = 1,size=rep(1,nrow(x)),sigma = 0,Effort = 1){     #simulate response

  u <- x%*%b

  if(LIKE == 'norm')return( rnorm(length(u),u,sqrt(sigma)) )
  if(LIKE == 'pois'){
  	 u <- exp(u)*Effort
  	 return( rpois(nrow(x),u) )
  } 
  if(LIKE == 'binom'){
  	 u <- inv.logit(u)
  	 return( rbinom(nrow(x),1,u) )
  }
  if(LIKE == 'multinom'){
     zs <- apply(exp(u),1,sum)
     z1 <- 1/(1 + zs)
     zm <- exp(u)/ (1 + zs)
     u  <- cbind(zm,z1)
     return( myrmultinom(size,u) )
  }
  if(LIKE == 'mvnorm'){
    u <- myrmvnorm(n,u,sigma)
    return( u )
  }
  if(LIKE == 'mvnorm-multinom'){
    u <- myrmvnorm(n,u,sigma)
    zs <- apply(exp(u),1,sum)
    z1 <- 1/(1 + zs)
    zm <- exp(u)/ (1 + zs)
    u  <- cbind(zm,z1)
    return( myrmultinom(size,u) )
  }
  numeric(0)
}
##########################################################

sensIntercept <- function(bgibbs,xnames,ynames){     

  #sensitivity coeffs for multinomial regression
  #bgibbs - multinomial logit par chains
  #xnames - length-k vector of names for input x 
  #ynames - length-(r-1) vector of names for response y

  r <- length(ynames)
  k <- length(xnames)

  fullnames <- as.vector(outer(xnames,ynames[-r],paste,sep='-')) #columns of bgibbs

  sbgibbs <- numeric(0)

  speccol <- matrix(unlist(strsplit(fullnames,'-')),ncol=2,byrow=T)[,2]

  for(h in 1:(r-1)){
    wh1 <- which(speccol == ynames[h])
    wh2 <- wh1[ !wh1 %in% c(grep('health',fullnames),grep('gap',fullnames)) ]
    wh0 <- wh1[ !wh1 %in% wh2 ]
    fmean  <- apply(bgibbs[,wh2],1,mean)
    int    <- sweep(bgibbs[,wh2],1,fmean) 
    int    <- cbind(fmean,bgibbs[,wh0],int)
    colnames(int)[1] <- paste('int',ynames[h],sep='-')
    sbgibbs <- cbind(sbgibbs,int)
  }
  sbgibbs
}
#######################################################
multinomLike <- function(y,x,b){  #log likelihood multinomial logit
	
  tiny <- 1e-20
  huge <- 1 - tiny

     z <- x%*%b
     zs    <- apply(exp(z),1,sum)
     z1    <- 1/(1 + zs)
     zm    <- exp(z)/ (1 + zs)
     z2  <- cbind(zm,z1)
     z2[z2 < tiny] <- tiny
     z2[z2 > huge] <- huge
     y*log(z2)
}
####################################################
bmultiProp <- function(b = matrix(0,k,r-1),pBVar=diag(.1,k*(r-1))){  
	
    bvec <- as.vector(b)
    if(length(bvec) < 100) cvec  <- myrmvnorm(1,t(bvec),pBVar)
    if(length(bvec) >= 100)cvec <- rnorm(length(bvec),t(bvec),sqrt(diag(pBVar)) )
    c    <- matrix(cvec,nrow(b),ncol(b))
    
    if('lob' %in% ls()){                         #if lob and hib available, use tnorm.mvt
    	cvec <- tnorm.mvt(t(bvec),pBVar,lob,hib)
      c    <- matrix(cvec,nrow(b),ncol(b))
    }

  list(c = c, cvec = cvec)
}

####################################################
bUpdateMNom <- function(x,y,b,pBVar){

  bvec <- as.vector(b)

  tmp  <- bmultiProp(b,pBVar)
  c    <- tmp$c
  cvec <- tmp$cvec
  
  pnow <- multinomLike(y,x,b)
  pnew <- multinomLike(y,x,c)

  pnow <- sum(pnow) + mydmvnorm(bvec,priorB,priorVB,log=T)
  pnew <- sum(pnew) + mydmvnorm(cvec,priorB,priorVB,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)b <- c
  b
}
###########################################33
bUpdateGibbs <- function(x,y,b,LIKE,sigma = 0,pBVar=0){

  if(LIKE == 'norm')    return( bUpdateNorm(x,y,sigma) )
  if(LIKE == 'multinom')return( bUpdateMNom(x,y,b,pBVar) )
  if(LIKE %in% c('mvnorm','mvnorm-multinom'))return( bUpdateMVNorm(x,y,b,sigma) )

  b <- matrix(b,length(b),1)
  c <- t(myrmvnorm(1,t(b),pBVar))

  znow <- x%*%b
  znew <- x%*%c

  if(LIKE == 'pois'){
     pnow <- dpois(y,exp(znow),log=T)
     pnew <- dpois(y,exp(znew),log=T)
  }
  if(LIKE == 'binom'){
     pnow <- dbinom(y,1,inv.logit(znow),log=T)
     pnew <- dbinom(y,1,inv.logit(znew),log=T)
  }

  pnow <- sum(pnow) + mydmvnorm(t(b),priorB,priorVB,log=T)
  pnew <- sum(pnew) + mydmvnorm(t(c),priorB,priorVB,log=T)

  a <- exp(pnew - pnow)
  z <- runif(1,0,1)
  if(z < a)b <- c
  b
}
####################################################
sigmaUpdate <- function(x,mu,s1,s2){

  u1 <- s1 + length(y)/2
  u2 <- s2 + .5*sum( (y - mu)^2 )
  1/rgamma(1,u1,u2)
}

######################################
wishsamp <- function(x,yy,b){   #sample from Inv Wishart

   r    <- ncol(b)
   scp  <- crossprod((yy - x %*% b))
   vmat <- solve(scp + prior.W*prior.WDF)
   v2   <- ns + prior.WDF
   stmp <- myrmvnorm(v2,matrix(0,v2,r),vmat)
   crossprod(stmp)
}
###################################
bUpdateMVNorm <- function(x,yy,b,sigma,alpha=0,lo=NULL,hi=NULL,X=NULL){  # update b's for mvnorm

  r      <- ncol(yy)
  k      <- ncol(x)

  if(is.null(X))X <- crossprod(x)
  
  bigv   <- invMat(X)      #multivariate scale invariate prior (Minka)
  smallv <- crossprod(x,yy)
                 
  mu     <- bigv%*%smallv/(alpha + 1)
  vaa    <- kronecker(sigma,bigv/(alpha + 1))   
  if(is.null(lo)) bg <- matrix( rmvnorm(1,as.vector(mu),sigma=vaa) ,k,r,byrow=F)
  if(!is.null(lo))bg <- matrix(tnorm.mvt(as.vector(b),as.vector(mu),vaa,lo,hi,times=1),k,r)
  
  colnames(bg) <- colnames(yy)
  rownames(bg) <- colnames(x)
  bg
}

#################################



diam2mass <- function(d,alloB,alloL){  #allo has c(int,slope) on log10 scale

  b <- 10^alloB[,1] *d^alloB[,2]
  l <- 10^alloL[,1] *d^alloL[,2]
  l[is.finite(l) & l > b] <- b[is.finite(l) & l > b]*.9
  cbind(b,l)

}

stemMass2diam <- function(stem,leaf,alloB){

  10^( (log10(stem + leaf) - alloB[,1])/alloB[,2] )
} 

mass2diam <- function(mass,allo){

  10^((log10(mass) - allo[,1])/allo[,2])
}


#################################
ysampMvnormMultinom <- function(x,y,z=NULL,b,sigma){  

  #sample y on MVlogit scale
  #z - counts, same dimensions as y

  r  <- ncol(y)
  k  <- nrow(b)
  n  <- nrow(y)
  
  propy <- matrix(rnorm(length(y),y,.01),n,r,byrow=F)
  
  pnow <- rep(0,n)
  pnew <- rep(0,n)
  
  for(i in 1:n){
    pnow[i] <- mydmvnorm(y[i,],(x %*% b)[i,], sigma,log=T)
    pnew[i] <- mydmvnorm(propy[i,],(x %*% b)[i,], sigma,log=T)
   }

  zs    <- apply(exp(y),1,sum)
  z1    <- 1/(1 + zs)
  zm    <- exp(y)/ (1 + zs)
  znow  <- cbind(zm,z1)

  zs    <- apply(exp(propy),1,sum)
  z1    <- 1/(1 + zs)
  zm    <- exp(propy)/ (1 + zs)
  znew  <- cbind(zm,z1)

  if(!is.null(z)){
    pnow <- pnow + z*log(znow)
    pnew <- pnew + z*log(znew)
  }

  pnow <- apply(pnow,1,sum)
  pnew <- apply(pnew,1,sum)

  a  <- exp(pnew - pnow)
  zz <- runif(length(a),0,1)
  y[zz < a,] <- propy[zz < a,]
  accept <- length(zz[zz < a])

  list(y = y, a = accept)
}

####################################################
getSens <- function(xv,bchain){

  #names(xv) are repeated for each class in bgibbs
  #ncol(bgibbs) = (r - 1)*length(xv)

  kk <- length(xv)

  nsim <- 2000
  sens <- matrix(NA,nsim,ncol(bchain))
  colnames(sens) <- colnames(bchain)

  for(j in 1:nsim){

    gj  <- sample(ng,1)
    bgj <- matrix(bchain[gj,],kk,r-1)
    sens[j,] <- as.vector( multinomSens(bgj,xvals=xv) )
  }
  sens
}

###############################################
plotSens <- function(sens,ytic,xnames,ynames,SIGONLY=F){  #plot multinomial sensitivities

# sens - output from getSens

  colF <- colorRampPalette(c('darkblue','blue','green','yellow','orange','red'))

  if('int' %in% xnames)xnames <- xnames[xnames != 'int']

  r <- length(ynames)
  k <- length(xnames)

  fullnames <- as.vector(outer(xnames,ynames[-r],paste,sep='-')) #columns of bgibbs

  tmp <- processPars(sens)$summary 

  par(mfrow=c(1,1),bty='n')

  xtic <- c(0:k)+.5

  plotSetup(xtic,ytic,xvals=rep(' ',length(xtic)),ylabel='Sensitivity')

  text(c(1:k),ytic[1],xnames,srt=90,cex=1.2)

  specindex <- matrix(unlist(strsplit(fullnames,'-')),ncol=2,byrow=T)[,2]

  jseq <- seq(-1,1,length=r)*.3
  cols <- colF(r-1)

  sigTable <- matrix(0,r-1,k)
  rownames(sigTable) <- ynames[-r]
  colnames(sigTable) <- xnames

  wwj <- matrix( unlist(strsplit(rownames(tmp),'-')),ncol=2,byrow=T)[,2]

  for(j in 1:(r-1)){

	tj <- tmp[wwj == ynames[j],]
        w0 <- grep('int',rownames(tj))
        if(length(w0) > 0)tj <- tj[-w0,]
	if(length(tj) == 0)next
        neg <- which(tj[,2] < 0 & tj[,3] < 0)
        pos <- which(tj[,2] > 0 & tj[,3] > 0)
        sigTable[j,neg] <- -1
        sigTable[j,pos] <- 1
	for(jk in 1:nrow(tj)){
                if(SIGONLY & !jk %in% neg & !jk %in% pos)next
		lines( c(jk,jk)+jseq[j],tj[jk,2:3],lwd=4,col=cols[j])
		points( jk+jseq[j],tj[jk,1],col=cols[j],pch=3)
	}
  }	


  pathScore <- rowSums(sigTable[,-c(1:3)])
  sigTable <- cbind(sigTable,pathScore)

  hostScore <- colSums(sigTable)
  sigTable  <- rbind(sigTable,hostScore)

  print(sigTable)

  legend('topleft',ynames[-r],text.col=cols[-r],cex=.9,bty='n',ncol=4)
  dev.copy2pdf(file='multinomModel.pdf')

  invisible(sigTable)
}

#######################################################
multinomSens <- function(bb,x=NULL,xvals=NULL){  #sensitivity coefficients multinomial logit
	
  tiny <- 1e-20
  huge <- 1 - tiny
  kk   <- nrow(bb)
  rr   <- ncol(bb)

  if(is.null(xvals))xvals <- matrix(apply(x,2,mean),1)

  z     <- xvals%*%bb
  zs    <- apply(exp(z),1,sum)
  zm    <- exp(z)/ (1 + zs)
  theta <- matrix(zm,kk,rr,byrow=T)
  bb*theta*(1 - theta)
}
####################################################

prob2Logit <- function(y){     #fractions to multivar logit      
  
  log(y[,-r]/(1 - apply(y[,-r],1,sum)))
  
}

#######################################
plotPosteriorOrder <- function(xx,yy,x1 = NULL, y1 = NULL, 
                               sep=max(yy)/nrow(yy)/2,xtic,ymax=NULL,xlab=' ',ALTLABEL=F,textSize=1){    

  #  y - each row is a density
  #  x - each row is the value for the density in y
  #  x1, y1 - second densities to be plotted with x,y
  #  names are rownames for y
  #  sep - vertical separation between densities
  
  nc      <- nrow(xx)
  wsord   <- apply(yy,1,which.max)
  xrw     <- range(xx)
  yrw     <- round(sep*11 + 1.2*max(yy),-2)
  word    <- order( xx[cbind(1:nrow(xx),wsord)] )
  specOrd <- rownames(yy)[word]

  if(!is.null(x1)){

    nc      <- nc + nrow(x1)
    wnord   <- apply(y1,1,which.max)
    xrw     <- range(c(xx,x1))
    yrw     <- round(sep*11 + 1.2*max(rbind(yy,y1)),-2)
    word    <- order( c(xx[cbind(1:nrow(xx),wsord)],x1[cbind(1:nrow(x1),wnord)]) )
    specOrd <- c(rownames(yy),rownames(y1))[word]

  }

  postext <- 4

  if(is.null(ymax))ymax <- yrw

  specOrd <- rev(unique(specOrd))
  jj      <- rev((c(1:length(specOrd)) - 1)*sep)
  if(max(jj) < ymax)ymax <- max(jj) + sep*2
  yvalw   <- seq(0,ymax,length=3)

  colS     <-  colorRampPalette( c("darkblue","orange","red","brown","black") )
  colSpecs <- colS(nc)

  plotSetup(xtic,yvalw,xlabel=xlab,yvals=rep(' ',length(yvalw)),ylabel='Density')

  for(j in 1:length(specOrd)){

    w0 <- which(rownames(yy) == specOrd[j])

    dx <- diff(xx[w0,])[1]
    cj <- cumsum(dx*yy[w0,])
    rj <- xx[w0,findInterval(c(.025,.975),cj)]

    yj <- jj[j] + yy[w0,]
    yj[yj > ymax] <- ymax

    if(rj[1] > 0)lines(c(xtic[1],max(xtic)),c(jj[j],jj[j]),lwd=6,col='white')
    if(rj[2] < 0)lines(c(xtic[1],max(xtic)),c(jj[j],jj[j]),lwd=6,col='wheat')
  }

  abline(v=0,col='grey',lwd=3)

  for(j in 1:length(specOrd)){

    w0 <- which(rownames(yy) == specOrd[j])
    w1 <- numeric(0)
    if(!is.null(y1))w1 <- which(rownames(y1) == specOrd[j])

    if(length(c(w0,w1)) == 0)next
  
    xl <- numeric(0)
    if(length(w0) > 0){
       xl <- range(xx[w0,])
  #     if(diff(xl) < .1){
         yj <- jj[j] + yy[w0,]
         yj[yj > ymax] <- ymax
         lines(xx[w0,],yj,lwd=2,col=colSpecs[j])
         polygon(c(xx[w0,],xx[w0,1]),c(yj,yj[1]),border=colSpecs[j])
             polygon(c(xx[w0,],xx[w0,1]),c(yj,yj[1]),border=colSpecs[j],col=colSpecs[j])
  #     }
    }
    if(length(w1) > 0){
       x2 <- range(x1[w1,])
       if(diff(x2) < .1){
         yj <- jj[j] + y1[w1,]
         yj[yj > ymax] <- ymax
         lines(x1[w1,],yj,lwd=2,col=colSpecs[j])
         polygon(c(x1[w1,],x1[w1,1]),c(yj,yj[1]),border=colSpecs[j])
         xl <- range(c(xl,x2))
       }
    }
  #  lines(xl,c(jj[j],jj[j]),col=colSpecs[j],lty=2,lwd=2)

    yt <- 1.2*max(yy[w0,])
    xt <- xx[w0,which.max(yy[w0,]) - 12]

    xxl <- max(xl)

    if(!ALTLABEL){
      postext <- 4
      xxl <- max(xl)
    }

    if(xxl > max(xtic) & !ALTLABEL){
        xxl <- min(xl)
        postext <- 2
    }

   if(postext == 4)xxl <- max(xl)
   if(postext == 2)xxl <- min(xl)

   text(xxl,jj[j],specOrd[j],col=colSpecs[j],pos=postext,cex=textSize)

   if(ALTLABEL & postext == 4){
     postext <- 2
   } else {postext <- 4}

  }

  scc <- signif(sep,1)
  ys2 <- yvalw[length(yvalw)]
  ys1 <- ys2 - scc
  arrows(xtic[1],ys1,xtic[1],ys2,angle=90,length=.06,lwd=2,code=3)
  text(xtic[1],ys2-.5*scc,scc,pos=4,cex=textSize)

}

getScoreNorm <- function(x,mu,xvar){  #Gneiting and Raftery's proper scoring rule

  #outcome x, prediction mean variance (mu, xvar)

  - ( (x - mu)^2)/xvar - log(xvar)

}

getScoreBinary <- function(x,p){    #binary outcome x (0,1), probability p

  (log(p))^x*(log(1 - p))^(1 - x)

}


