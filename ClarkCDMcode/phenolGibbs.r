

source('bbFunctions.r')

#library(Matrix)

#library(mvtnorm)
#library(tmvtnorm)
#library(gmm)
#library(utils)
#library(stats4)

threshMu <- 0
threshSd <- 1

colVector <- c("brown","tan4","green","darkgreen")

colF      <-  colorRampPalette(colVector)
colVals   <- colF(100)
colStates <- colF(nclass)
stateSeq  <- colVals


factors <- c('seed','site','gap')

jspec <- which(specList == gibbsSpec)
pData <- pDataAll[[jspec]]
specData <- sDataAll[[jspec]]

if(length(pData) == 0)next



dt  <- 2   #day intervals to model
nyr <- length(yrVec)

tmp <- timeSetup()
#germ  <- tmp$germ        #by jdAll date
#died  <- tmp$died
pData <- tmp$pData
specData <- tmp$specData
jdYr     <- tmp$jdYr
jdAll    <- tmp$jdAll
sampleDates <- tmp$sampleDates
sampleIndex <- tmp$sampleIndex
jdIndex    <- tmp$jdIndex
yrIndex    <- tmp$yrIndex
dtvec      <- c(0,diff(jdAll))      #dt from t to t+1
yrBreak    <- which(dtvec > dt)       #yr transitions
#yrBreak   <- c(yrBreak,yrVec[yrVec > max(yrBreak)]+1)
germIndex   <- tmp$germTime
germTime    <- germIndex
diedIndex   <- tmp$diedTime
diedTime   <- diedIndex
germYr     <- tmp$germYr
diedYr     <- tmp$diedYr
treeYr     <- sum(diedYr - germYr)

nt <- length(jdAll)

if(length(pData) == 0)next
if(nrow(pData) < 10)next


siteCode <- seedCode <- c('DF','HF')  #factors have HF = 1


#print(specList[jspec])

n   <- nrow(pData)
nd  <- ncol(pData)
nyr <- max(yrIndex)
years      <- c(2009:2020)[1:nyr]

nclass <- max(pData,na.rm=T)
times <- c(1:nt)
germIndex <- cbind(c(1:n),germIndex+1)
diedIndex <- cbind(c(1:n),diedIndex+1)

tmat <- matrix(times,n,nt,byrow=T)

if('seed' %in% xnames){    #south 0, north 1
  tmp <- seedOrigin()
  seed <- tmp$seed
  if(sum(seed) < 10 | sum(1 - seed) < 10){
    wg <- grep('seed',xnames)
    xnames <- xnames[-wg]
  }
}
gap <- specData[,'gap']
  

#print(xnames)

#priors


xsummer <- rep(0,length(xnames))
xsummer[xnames == 'AT'] <- 10
xsummer[xnames == 'SM'] <- -.3


tmp <- initialStates()
   sstate <- tmp$sstate
   hInit  <- tmp$h
   dataMat <- tmp$dataMat
   hIndex <- tmp$hIndex
   xIndex <- tmp$xIndex
   fix1   <- tmp$fix1
   fix6   <- tmp$fix6


pcor <- matrix(NA,n,nt)
pcor[hIndex] <- 1
pcor[,nt] <- 1
hInit  <- hInit*pcor
sstate <- sstate*pcor

h <- hInit
dataIndex <- which(is.finite(dataMat),arr.ind=T)
	 
  #data matrix and initialize########################################
	 

  cham <- matrix(unlist(strsplit(as.character(specData[,'Chamber']),'_')),
                  ncol=2,byrow=T)[,1]
  site <- specData[,'site']

  nx <- length(xnames)
  
  tmp <- makeXmat()
   xmat   <- tmp$xmat
   x3d    <- tmp$x3d
   meanX  <- tmp$meanX
   rangeX <- tmp$rangeX
   xnames <- colnames(xmat)
   nx     <- length(xnames)

   pathName <- 'model'
   for(k in 2:nx)pathName <- paste(pathName,xnames[k],sep='-')
  tmp <- inProcEnvData(PLOT=T,warmTemp=3,chillTemp=3)

#  dev.print(device=postscript,file=paste(gibbsSpec,'Env.ps',sep='_'),
#          width=6,horizontal=F)

  if(!pathName %in% list.files())dir.create(pathName)
  filename <- paste(pathName,'/',gibbsSpec,sep='')

  if(paste(filename,'.RData',sep='') %in% list.files())next

ddmat      <- getDD(0)
firstDates <- date2Phen()           #first obs date and DD by stage
if(length(firstDates) == 0)next     #none yet at BB


obsPrior <- matrix(1e-10,nclass,nclass)
diag(obsPrior) <- n*nt
obsPrior[row(obsPrior) == col(obsPrior)-1] <- 0.1
obsPrior[row(obsPrior) == col(obsPrior)+1] <- 0.1
obsPrior[row(obsPrior) == col(obsPrior)-2] <- .000001
obsPrior[row(obsPrior) == col(obsPrior)+2] <- .000001

obsErr <- obsPrior
obsErr <- obsErr/matrix(apply(obsErr,2,sum),nclass,nclass,byrow=T)

lo <-  c(5,20,30,55,95)
hi <- c(10,30,55,80,99)
hseq <- seq(0,100,by=1)
maxSlope <- -.20           #slope must be steeper than this

breakLims <- rbind(lo,hi)
breaks   <- apply(breakLims,2,mean)
priorBQ  <- makePriorB(c(1,10),c(-2,-.6),breaks,hseq)
vb       <- rep(10,10)
vb[c(1,5,6,10)] <- 1/n/nt
priorVBQ <- diag(vb)
plotLogit(breakLims,breaks,hseq,priorBQ)

priorB   <- matrix(0,length(xnames),1)
priorB[xnames == 'AT',1] <- .05
priorB[xnames == 'SM',1] <- .5
priorB[xnames == 'h2',1] <- -2/100

priorVB  <- rep(10,length(xnames))
#priorVB[xnames == 'h'] <- 1/n/n/nd/nd
priorVB <- diag(priorVB)

priorIVB <- solve(priorVB)
loBA      <- priorB*0 - 100
loBA[xnames == 'AT'] <- 0
#loBA[xnames == 'h']  <- 0
#loBA[xnames == 'gap'] <- 0

hiBA      <- priorB*0 + 100
#hiBA[xnames == 'h2'] <- 0


maxS <- 2    # dh/dt approx 1 degree per day
mus <- .1
s1  <- n
s2  <- mus*(s1 - 1)

rangeTemp <- matrix(NA,n,nt)
rangeTemp[xIndex] <- xmat[,'AT']
rangeTemp <- apply(rangeTemp,2,quantile,c(.5,.025,.975),na.rm=T)
plot(jdAll,rangeTemp[1,],type='l')
abline(v=jdAll[yrBreak],lty=2)


bq     <- priorBQ
breakg <- breaks

sigma <- 1

missX   <- which(is.na(xmat),arr.ind=T)
missRow <- sort(unique(missX[,1]))


z <- makeY(h[xIndex],h[hIndex],.1)
wz <- which(is.finite(z))

bprop <- solve(crossprod(xmat[wz,]))
bg    <- bprop%*%crossprod(xmat[wz,],z[wz])
bg[bg < loBA] <- loBA[bg < loBA] + .01
bg[bg > hiBA] <- hiBA[bg > hiBA] - .01

hprop <- .01
hpropSD <- .1

t11     <- table(hIndex[,1])
iprint <- as.numeric(names(t11)[which.max(t11)])

#iprint <- 221

#iprint <- 1

print1 <- yrBreak[findInterval(germIndex,yrBreak)[iprint] + 1] 
print2 <- c(yrBreak,nt)[findInterval(germIndex,yrBreak)[iprint] + 2] 
#print2 <- nt
wii    <- c(print1:print2)


#plotData()
#dev.print(device=postscript,file=paste(filename,'Data.ps',sep='_'),
#          width=6,horizontal=F)

###########################################################################

bqgibbs <- matrix(NA,ng,length(bq))
colnames(bqgibbs) <- as.vector( outer( 1:(nclass-1),1:2,paste,sep=',') )
bggibbs <- matrix(NA,ng,length(bg)); colnames(bggibbs) <- xnames
sgibbs  <- rep(0,ng)

hg <- hg2 <- h*0
ogibbs <- matrix(NA,ng,length(obsErr))

predObs <- h*0

bcount <- 0
acount <- 0
pv1 <- log(.01)
pv2 <- log(.0005)

whichxh <- which(xnames == 'h')
xmatTemp <- matrix(0,n,nt)

counts <- counth <- stateSum <- h*0
sumdev <- 0

for(g in 1:ng){

print(g)

  for(j in 1:5){
    tmp <- statePars(bq,sstate[hIndex],h[hIndex])
    bq     <- tmp$bgg
    breakg <- tmp$breakg
    aa     <- tmp$accept
    if(aa == 1)acount <- 0
    if(aa == 0)acount <- acount + 1
  }
 
  tmp <- updateB()
  bg    <- tmp$bg
  sigma <- tmp$sigma

  tmp <- updateStates() 
  h      <- tmp$h        #latent continuous state
  sstate <- tmp$s        #latent discrete state
  pObs   <- tmp$predObs
  counts <- counts + tmp$counts
  counth <- counth + tmp$counth
 # h <- inc2H(h,hIndex,xIndex)
  devg <- tmp$dev

  if('h' %in% xnames){
    xmat[,'h'] <- h[xIndex]
    xmat[,'h2'] <- h[xIndex]^2
    xmatTemp[xIndex] <- xmat[,'h']
    x3d[,,'h'] <- xmatTemp
    xmatTemp[xIndex] <- xmat[,'h2']
    x3d[,,'h2'] <- xmatTemp
  }

  obsErr <- updateObsError()
#print('obs')
  bqgibbs[g,] <- bq
  bggibbs[g,] <- bg
  sgibbs[g]   <- sigma

  if(g >= burnin){
    hg  <- hg + h
    hg2 <- hg2 + h^2
    predObs <- predObs + pObs
    stateSum <- stateSum + sstate
    sumdev <- sumdev + devg
  }

  ogibbs[g,] <- obsErr
  
#  print(sigma)
#  print(bg)
#  print(bq)
  print( round(rbind(h[iprint,wii],sstate[iprint,wii],dataMat[iprint,wii])),1)   
#   print(rbind(sstate[iprint,wii],dataMat[iprint,wii]))

}

ng <- g-1
ntot <- ng - burnin

predObs <- predObs*pcor
stateMu <- stateSum*pcor

bqgibbs <- bqgibbs[burnin:ng,]
bggibbs <- bggibbs[burnin:ng,]
sgibbs  <- sgibbs[burnin:ng]
ogibbs  <- ogibbs[burnin:ng,]

bqout <- processPars(bqgibbs,CPLOT=T)$summary

bgout <- processPars(bggibbs,xtrue=rep(0,ncol(bggibbs)),CPLOT=T)$summary

#bgout <- processPars(bggibbs,xtrue=rep(0,ncol(bggibbs)),DPLOT=T)$summary

sout <- processPars(cbind(sgibbs,ogibbs[,1:5]),CPLOT=T)$summary

sout <- processPars(cbind(sgibbs,ogibbs[,1:5]),DPLOT=T)$summary

outTable <- parTable()$parOut

predvals <- predObs/ntot
stateMu  <- stateSum/ntot

oerror <- matrix(apply(ogibbs,2,mean),nclass,nclass)

#DD calculations

hmean <- hg/ntot
hse   <- sqrt(hg2/ntot - hmean^2)
hmean[hmean > 100] <- 100

smean <- stateSum/ntot




save(list = ls(all=TRUE), file = paste(filename,'.RData',sep=''))



