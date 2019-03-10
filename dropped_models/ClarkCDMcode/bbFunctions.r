
library(cluster)

envCol2treatment <- function(siteName,mat,cvar){  #column names from env files to treatment

  gcol <- grep(cvar,colnames(mat))
  dcol <- unlist(strsplit(colnames(mat)[gcol],cvar))
  
  if(siteName == 'DF'){
          amb   <- c('G01','G04','G07','S03','S04','S08','G1','G4','G7','S3','S4','S8')
          plus3 <- c('G03','G05','G09','S02','S05','S09','G3','G5','G9','S2','S5','S9')
          plus5 <- c('G02','G06','G08','S01','S06','S07','G2','G6','G8','S1','S6','S7')
          cont  <- c('G10','G11','G12','S10','S11','S12')
  }
  if(siteName == 'HF'){
          amb   <- c('G02','G04','G07','S02','S05','S07')
          plus3 <- c('G01','G06','G08','S03','S06','S08')
          plus5 <- c('G03','G05','G09','S01','S04','S09')
          cont  <- c('G10','G11','G12','S10','S11','S12')
  }

    zero  <- which(dcol %in% amb | dcol %in% cont)
    three <- which(dcol %in% plus3)
    five  <- which(dcol %in% plus5)
    gaps  <- grep('G',dcol)
    shade <- grep('S',dcol)

  list(zero = gcol[zero], three = gcol[three], five = gcol[five], gaps = gcol[gaps], shade = gcol[shade])
}


getHotCold <- function(threshold,cDates,wDates,jd,mat){

    chillDays <- warmDays <- chillPlot <- warmPlot <- jdw <- jdc <- numeric(0)
    wd  <- which(jd > cDates[1] & jd < cDates[2])

    if(length(wd) == 0){

      return(list(chillDays = chillDays, warmDays = warmDays, chillPlot = chillPlot, 
             warmPlot = warmPlot, jdwarm = jdw, jdcold = jdc))
    }

    jdc <- jd[wd]
    ct  <- threshold - apply(mat[wd,],1,mean,na.rm=T)
    ct[ct < 0] <- 0
    chillDays <-  -cumsum(ct)

    ct <- threshold - mat[wd,]
    ct[ct < 0] <- 0
    chillPlot <- apply(ct,2,sum,na.rm=T)

    wd  <- which(jd > wDates[1] & jd < wDates[2])
    jdw <- jd[wd]
    wt  <- apply(mat[wd,],1,mean,na.rm=T) - threshold
    wt[wt < 0] <- 0
    warmDays <- cumsum(wt)

    wt <- mat[wd,] - threshold
    wt[wt < 0] <- 0
    warmPlot <- apply(wt,2,sum,na.rm=T)

    list(chillDays = chillDays, warmDays = warmDays, chillPlot = chillPlot, warmPlot = warmPlot, jdwarm = jdw, jdcold = jdc)
}


inProcEnvData <- function(PLOT=F,warmTemp=3,chillTemp=3){

  path    <- 'processedData/'
  
  DF   <- read.csv(paste(path,"daily-Duke-AT.csv",sep=''),header=T)
  DFSM <- read.csv(paste(path,"daily-Duke-SM.csv",sep=''),header=T)
  HF   <- read.csv(paste(path,"daily-Harvard-AT.csv",sep=''),header=T)
  HFSM <- read.csv(paste(path,"daily-Duke-SM.csv",sep=''),header=T)
  DFST <- read.csv(paste(path,"daily-Duke-ST.csv",sep=''),header=T)
  HFST <- read.csv(paste(path,"daily-Harvard-ST.csv",sep=''),header=T)

  tmp <- envCol2treatment('DF',DF,'AT')
  ambDF <- tmp$zero
  p3DF  <- tmp$three
  p5DF  <- tmp$five
  gapDF <- tmp$gap
  shaDF <- tmp$shade

  tmp <- envCol2treatment('HF',HF,'AT')
  ambHF <- tmp$zero
  p3HF  <- tmp$three
  p5HF  <- tmp$five
  gapHF <- tmp$gap
  shaHF <- tmp$shade

  tmp <- envCol2treatment('DF',DFST,'ST')
  ambDFST <- tmp$zero
  p3DFST  <- tmp$three
  p5DFST  <- tmp$five
  gapDFST <- tmp$gap
  shaDFST <- tmp$shade

  tmp <- envCol2treatment('HF',HFST,'ST')
  ambHFST <- tmp$zero
  p3HFST  <- tmp$three
  p5HFST  <- tmp$five
  gapHFST <- tmp$gap
  shaHFST <- tmp$shade

  y2plot <- sort(unique(c(DF[,'yeary'],HF[,'yeary'])))

  if(PLOT){
  
  colVec <- c('black','grey55','blue2','slateblue2')

  par(mfrow=c(3,1),bty='n',mar=c(1,4,1,2),cex.axis=1.4,cex.lab=1.4)

  ylimit <- c(-15,35)
  xlimit <- c(0,max(jdAll) + 120)

  plot(c(-100,-100),c(0,0),xlim=xlimit,ylim=ylimit,xaxt='n',xlab=' ',ylab='degrees C')
  axis(side=1,at=moVec,labels=F)

  rect(xlimit[1],ylimit[1],xlimit[2],ylimit[2],col='azure2',border=NA)

  abline(h=0,lwd=3,col='white')
  abline(v=yrVec,lwd=3,col='white')

    jdd <- DF[,'ydate'] 
    tda  <- apply(DF[,ambDF],1,mean,na.rm=T)
    lines(jdd,tda)

    tdp  <- apply(DF[,c(p3DF,p5DF)],1,mean,na.rm=T) - tda
    lines(jdd,tdp,lwd=2,col=colVec[2])

    jdh <- HF[,'ydate'] 
    tha  <- apply(HF[,ambHF],1,mean,na.rm=T)
    lines(jdh,tha,col=colVec[3])

    thp  <- apply(HF[,c(p3HF,p5HF)],1,mean,na.rm=T) - tha
    lines(jdh,thp,col=colVec[4],lwd=2)

    text(.15*max(jdAll),28,'a) Air Temperature',cex=1.3)


  plot(c(-100,-100),c(0,0),xlim=xlimit,ylim=ylimit,xaxt='n',xlab=' ',ylab='degrees C')
  axis(side=1,at=moVec,labels=F)

  rect(xlimit[1],ylimit[1],xlimit[2],ylimit[2],col='azure2',border=NA)

  abline(h=0,lwd=3,col='white')
  abline(v=yrVec,lwd=3,col='white')

    jdd <- DFST[,'ydate'] 
    tda  <- apply(DFST[,ambDF],1,mean,na.rm=T)
    lines(jdd,tda)

    tdp  <- apply(DFST[,c(p3DF,p5DF)],1,mean,na.rm=T) - tda
    lines(jdd,tdp,col=colVec[2],lwd=2)

    jdh <- HFST[,'ydate'] 
    tha  <- apply(HFST[,ambHF],1,mean,na.rm=T)
    lines(jdh,tha,col=colVec[3])

    thp  <- apply(HFST[,c(p3HF,p5HF)],1,mean,na.rm=T) - tha
    lines(jdh,thp,,col=colVec[4],lwd=2)

    text(.15*max(jdAll),28,'b) Soil Temperature',cex=1.3)



  ylimit <- c(-8000,7000)

  par(mar=c(4,4,1,2))

  plot(c(-1000,-1000),c(0,0),xlim=xlimit,ylim=ylimit,xaxt='n',xlab='Year',ylab='DD',yaxt='n')


  rect(xlimit[1],ylimit[1],xlimit[2],ylimit[2],col='azure2',border=NA)

  axis(side=1,at=c(0,yrVec[-length(yrVec)])+365/2,labels=y2plot)
  axis(side=1,at=moVec,labels=F)


  axis(side=2,at=seq(0,5000,by=1000),labels=c(' ',' ',' ',' ',' ',5000))
  axis(side=2,at=seq(-5000,0,by=1000),labels=c(500,' ',' ',' ',' ',0),line=.5)

  abline(h=0,lwd=3,col='white')
  abline(v=yrVec,lwd=3,col='white')
  }

  startChill <- moVec[names(moVec) == 'oct']
  endChill   <- moVec[names(moVec) == 'feb']
  startChill <- startChill
  endChill   <- c(endChill[-1],endChill[length(endChill)] + 300)

  startWarm <- c(1,yrVec)
  endWarm   <- c(yrVec-70,yrVec[length(yrVec)] + 365 - 70)

  chillDates <- cbind(startChill,endChill)
  warmDates  <- cbind(startWarm,endWarm)  


  coldDFa <- coldHFa <- coldDFs <- coldHFs <-hotDF <- hotHF <- numeric(0)


  for(y in 1:length(y2plot)){

    tmp <- getHotCold(warmTemp,chillDates[y,],warmDates[y,],DF[,'ydate'],DF[,-c(1:3)])
      coldDFa <- rowBind(coldDFa,tmp$chillPlot,y2plot[y]+1)
      hotDF   <- rowBind(hotDF,tmp$warmPlot,y2plot[y])

    tmp <- getHotCold(warmTemp,chillDates[y,],warmDates[y,],DFST[,'ydate'],DFST[,-c(1:3)])
      coldDFs <- rowBind(coldDFs,tmp$chillPlot,y2plot[y]+1)

    tmp <- getHotCold(warmTemp,chillDates[y,],warmDates[y,],HF[,'ydate'],HF[,-c(1:3)])
      coldHFa <- rowBind(coldHFa,tmp$chillPlot,y2plot[y]+1)
      hotHF   <- rowBind(hotHF,tmp$warmPlot,y2plot[y])

   tmp <- getHotCold(warmTemp,chillDates[y,],warmDates[y,],HFST[,'ydate'],HFST[,-c(1:3)])
      coldHFs <- rowBind(coldHFs,tmp$chillPlot,y2plot[y]+1)


   tmp <- getHotCold(warmTemp,chillDates[y,],warmDates[y,],DF[,'ydate'],DF[,ambDF])
      if(PLOT)lines(tmp$jdcold,tmp$chillDays*10,lwd=2)
      if(PLOT)lines(tmp$jdwarm,tmp$warmDays,lwd=2)
      cold <- min(tmp$chillDays)

   tmp <- getHotCold(warmTemp,chillDates[y,],warmDates[y,],DF[,'ydate'],DF[,c(p3DF,p5DF)])
      if(PLOT)lines(tmp$jdcold,tmp$chillDays*10,lwd=2,col=colVec[2])
      if(PLOT)lines(tmp$jdwarm,tmp$warmDays,lwd=2,col=colVec[2])
      cold <- c(cold,min(tmp$chillDays))

   tmp <- getHotCold(warmTemp,chillDates[y,],warmDates[y,],HF[,'ydate'],HF[,ambHF])
      if(PLOT)lines(tmp$jdcold,tmp$chillDays*10,lwd=2,col=colVec[3])
      if(PLOT)lines(tmp$jdwarm,tmp$warmDays,lwd=2,col=colVec[3])
      cold <- c(cold,min(tmp$chillDays))

   tmp <- getHotCold(warmTemp,chillDates[y,],warmDates[y,],HF[,'ydate'],HF[,c(p3HF,p5HF)])
      if(PLOT)lines(tmp$jdcold,tmp$chillDays*10,lwd=2,col=colVec[4],lty=2)
      if(PLOT)lines(tmp$jdwarm,tmp$warmDays,lwd=2,col=colVec[4])
      cold <- c(cold,min(tmp$chillDays))

   tmp <- getHotCold(warmTemp,chillDates[y,],warmDates[y,],DFST[,'ydate'],DFST[,ambDFST])
      if(PLOT)lines(tmp$jdcold,tmp$chillDays*10,lwd=2,lty=2)
      cold <- c(cold,max(tmp$chillDays))

   tmp <- getHotCold(warmTemp,chillDates[y,],warmDates[y,],DFST[,'ydate'],DFST[,c(p3DFST,p5DFST)])
      if(PLOT)lines(tmp$jdcold,tmp$chillDays*10,lwd=2,col=colVec[2],lty=2)
      cold <- c(cold,min(tmp$chillDays))

   tmp <- getHotCold(warmTemp,chillDates[y,],warmDates[y,],HFST[,'ydate'],HFST[,ambHFST])
      if(PLOT)lines(tmp$jdcold,tmp$chillDays*10,lwd=2,col=colVec[3],lty=2)
      cold <- c(cold,min(tmp$chillDays))
      
   tmp <- getHotCold(warmTemp,chillDates[y,],warmDates[y,],HFST[,'ydate'],HFST[,c(p3HFST,p5HFST)])
      if(PLOT)lines(tmp$jdcold,tmp$chillDays*10,lwd=2,col=colVec[4],lty=2)
      cold <- c(cold,min(tmp$chillDays))

  if(PLOT){
  xlimit <- range(tmp$jdcold)
  dx <- diff(xlimit)/4/2
  xk <- xlimit[2] + 2*dx*c(1:4)

  for(k in 1:4){
    rect(xk[k] - dx,ylimit[1],xk[k] + dx,ylimit[1] - 5*cold[k],border=colVec[k],lwd=2)
    rect(xk[k] - dx,ylimit[1],xk[k] + dx,ylimit[1] - 5*cold[k+4],border=colVec[k],col=colVec[k],lwd=2)
  }
  }
   
}

  if(PLOT){
    legend('bottomleft',c('ambient Duke','elevated Duke','ambient Harvard','elevated Harvard','soil - dashed'),
            text.col=colVec,bty='n')

    text(max(jdAll)*.8,-4000,'chilling units',cex=1.2)
    text(max(jdAll)*.7,5000,'growing DD',cex=1.2)
    text(.15*max(jdAll),5000,'c) Degree days',cex=1.3)
  }
  
  list(DFAT = DF, DFSM = DFSM, HFAT = HF, HFSM = HFSM, coldHFs = coldHFs, coldDFs = coldDFs, 
       coldHFa = coldHFa, coldDFa = coldDFa, hotHF = hotHF, hotDF = hotDF)
}

#####################################
plotAdvance <- function(maxDate){

      par(mfrow=c(1,1),bty='n')
      md   <- maxDate
      xseq <- seq(-60,200,by=4)
      max2 <- hist(md[,'maxDate2'],plot=F,breaks=xseq)
      max3 <- hist(md[,'maxDate3'],plot=F,breaks=xseq)

      ddd  <- hist(md[,'bbDate'],plot=F,breaks=xseq)
      difd <- hist(md[,'maxDate3'] - md[,'bbDate'],plot=F,breaks=xseq)
      plot(max2$mids,max2$density,xlim=c(-50,200),type='s',ylim=c(0,.15),
           lwd=2,xaxt='n',xlab=' ',ylab='Density',col='grey')
      axis(1,at=(moVec[1:6]-15),labels=names(moVec)[1:6])
      axis(1,at=c(-50,-50,0))

      lines(max3$mids,max3$density,type='s',lwd=2,col='darkgreen')
      lines(ddd$mids,ddd$density,type='s',lwd=2)

      md1 <- rbind(quantile(md[md[,'site'] == 1,'maxDate2'],c(.5,.05,.95),na.rm=T),
                   quantile(md[md[,'site'] == 2,'maxDate2'],c(.5,.05,.95),na.rm=T))
      lines(md1[1,2:3],c(.06,.06),lwd=2)
      lines(md1[2,2:3],c(.063,.063),col='blue',lwd=2)
      points(md1[,1],c(.06,.063),col=c('black','blue'),pch=3)

      md2 <- rbind(quantile(md[md[,'site'] == 1,'bbDate'],c(.5,.05,.95),na.rm=T),
                   quantile(md[md[,'site'] == 2,'bbDate'],c(.5,.05,.95),na.rm=T))
      lines(md2[1,2:3],c(.07,.07),lwd=2)
      lines(md2[2,2:3],c(.073,.073),col='blue',lwd=2)
      points(md2[,1],c(.07,.073),col=c('black','blue'),pch=3)

      md3 <- rbind(quantile(md[md[,'site'] == 1,'maxDate2'] - md[md[,'site'] == 1,'bbDate'],c(.5,.05,.95),na.rm=T),
                   quantile(md[md[,'site'] == 2,'maxDate2'] - md[md[,'site'] == 2,'bbDate'],c(.5,.05,.95),na.rm=T))
      lines(md3[1,2:3],c(.06,.06),lwd=2)
      lines(md3[2,2:3],c(.063,.063),col='blue',lwd=2)
      points(md3[,1],c(.06,.063),col=c('black','blue'),pch=3)

      mg  <- which(colnames(md) %in% c('maxDate2','maxDate3'))                #include for stage 2 and 3
      md1 <- rbind(quantile(md[md[,'site'] == 1,mg],c(.5,.05,.95),na.rm=T),
                   quantile(md[md[,'site'] == 2,mg],c(.5,.05,.95),na.rm=T))

      mdd <- rbind(as.vector(t(md1)),as.vector(t(md2)),as.vector(t(md3)))
      rownames(mdd) <- c('maxdate','budbreak','difference')
      colnames(mdd) <- as.vector(outer(c('So','No'),colnames(md1),paste,sep='-'))
      mdd
}


######################################################
timeSetup <- function(){

  sampleDates <- as.numeric(colnames(pData))

  germ <- specData[,'germ']
  died <- apply(specData[,c('died','rightCensorDied')],1,max,na.rm=T) 
  names(died) <- NULL

#  lastMeasure <- matrix(as.numeric(colnames(pData)),nrow(pData),ncol(pData),byrow=T)
#  ptmp <- pData
#  ptmp[ptmp == 1] <- NA
#  lastMeasure <- apply(lastMeasure*(ptmp*0 + 1),1,max,na.rm=T)


  germYr <- findInterval(germ,c(0,yrVec))
  diedYr <- findInterval(died,c(0,yrVec))
  diedYr[is.na(specData[,'died'])] <- length(yrVec) + 1
  germYr[is.na(germYr)] <- 1
  diedYr[is.na(diedYr)] <- length(yrVec)


  dateSamp <- as.numeric(colnames(pData))
  keepi    <- rep(0,nrow(pData))
  keepy    <- matrix(NA,nrow(pData),length(yrVec))

  for(y in 2:length(yrVec)){

    wt <- which(dateSamp > yrVec[y-1] & dateSamp < yrVec[y])
    my <- apply(pData[,wt],1,max,na.rm=T)
    ww <- which(is.finite(my) & my > 2 & y > germYr)
    keepi[ww] <- keepi[ww] + 1
    keepy[ww,y] <- y
  }

  wl <- which(keepi > 0)
  kp <- keepy[wl,]

  first  <- apply(kp,1,min,na.rm=T)
  last   <- apply(kp,1,max,na.rm=T)
  diedYr <- diedYr[wl]
  germYr <- first - 1



#  last[last == length(yrVec)] <- length(yrVec) + 1
#  diedYr[last < diedYr]  <- last[last < diedYr]   #diedYr is last yr, not yr of death
#  germYr[first > (1+germYr)] <- first[first > (1+germYr)] - 1

 # wl <- which( (diedYr - germYr) > 0 & keepi > 0)

#  wl <- which(keepi > 0)

  pData    <- pData[wl,]
  specData <- specData[wl,]
#  germYr   <- germYr[wl]
#  diedYr   <- diedYr[wl]
#  germ     <- germ[wl]
#  died     <- died[wl]

  minDate <- 20   #last jd when must be in class 1
  maxDate <- 260  #first jd when must be in nclass

  endDate   <- c(0,yrVec) + maxDate
  startDate <- c(0,yrVec) + minDate

#  ws <- findInterval(min(germ,na.rm=T),startDate)
#  we <- findInterval(max(died,na.rm=T),endDate)

  sampleDates <- as.numeric(colnames(pData))
  smat <- pData
  smat[smat <= 1] <- NA
  lastDate <- max(sampleDates[apply(smat,2,sum,na.rm=T) > 0],na.rm=T)

  pData <- pData[,as.numeric(colnames(pData)) <= lastDate]
  sampleDates <- as.numeric(colnames(pData))

 # startDate <- startDate[ws:nyr]
 # endDate   <- endDate[ws:nyr]

  nyr <- length(startDate)

  jdYr <- seq(startDate[1],endDate[1],by=dt)               #days to model within a yr
  jdAll <- numeric(0)
  for(y in 1:nyr){
    jdAll <- c(jdAll,seq(startDate[y],endDate[y],by=dt))  #all days to model
  }

  jdAll <- jdAll[jdAll <= lastDate]

  sampleIndex <- match(sampleDates,jdAll)
  sint        <- findInterval(sampleDates,jdAll)
  sampleIndex[is.na(sampleIndex)] <- sint[is.na(sampleIndex)]

  wy <- rep(F,length(sampleIndex))
  for(y in 1:nyr){
    wy[sampleDates >=startDate[y] & sampleDates <= endDate[y]] <- T
  }


  pData <- pData[,wy]
  w0 <- which(apply(pData,1,sum,na.rm=T) > 0)
  pData <- pData[w0,]
  specData <- specData[w0,]
 # germYr <- germYr[w0]
 # diedYr <- diedYr[w0]
 # germ <- germ[w0]
 # died <- died[w0]

  sampleDates <- as.numeric(colnames(pData))
  sampleIndex <- findInterval(sampleDates,jdAll)

  lastDate <- max(sampleDates)
  www <- which(jdAll > lastDate)
  endDate[endDate > max(sampleDates)] <- max(sampleDates)

  pData[pData == 0] <- NA
  nclass <- 6
  pData[pData > nclass] <- NA

  yrIndex <- findInterval(sampleDates,c(0,yrVec))
  jdIndex <- findInterval(jdAll,c(0,yrVec))
  nt  <- length(jdAll)
  jcheck <- c(0,jdIndex,nyr+1)

  germTime <- findInterval(germYr-1,jcheck)
  diedTime <- findInterval(diedYr,jdIndex)

  list(pData = pData, specData = specData, jdYr = jdYr, jdAll = jdAll,germYr = germYr, diedYr = diedYr,
       sampleDates = sampleDates, sampleIndex = sampleIndex, yrIndex = yrIndex, jdIndex = jdIndex,
       germTime = germTime, diedTime = diedTime)
}


#########################################

seedOrigin <- function(){

  seed <- rep(0,nrow(specData))
  seed[specData[,'SeedOrigin'] %in% northSeed ] <- 1
  seed[specData[,'SeedOrigin'] %in% localSeed & specData[,'site'] == seedCode[1]] <- 0
  seed[specData[,'SeedOrigin'] %in% localSeed & specData[,'site'] == seedCode[2]] <- 1

  if(sum(seed,na.rm=T) < 10 | sum(1 - seed,na.rm=T) < 10){
     ws <- grep('seed',xnames)
     xnames <- xnames[-ws]
  }

  seedTable <- table(specData[,'SeedOrigin'],seed)
  if(ncol(seedTable) == 1){
    rn <- rownames(seedTable)
    wc <- as.numeric(colnames(seedTable))
    nc <- seedTable*0
    if(wc == 0)seedTable <- cbind(seedTable,nc)
    if(wc == 1)seedTable <- cbind(nc,seedTable)
  }

  w1 <- which(apply(seedTable,1,sum) > 0)
  st <- seedTable[w1,]
  if(!is.matrix(st)){
    st <- matrix(st,nrow=1)
    rownames(st) <- rownames(seedTable)[w1]
  }

  seedTable <- st
  colnames(seedTable) <- c('So','No')
  nSource <- seedTable
  nSource[nSource == 0] <- NA
  nSource <- nSource*0 + 1
  nSource <- apply(nSource,2,sum,na.rm=T)
  seedTable <- rbind(seedTable,nSource)

  list(seed = seed, xnames = xnames, sourceTable = seedTable)
}
#########################################

initialStates <- function(){

  dataMat <- h <- sstate <- matrix(NA,n,nt)
  dataMat[,sampleIndex] <- pData

  sstate[is.finite(dataMat)] <- dataMat[is.finite(dataMat)]

   yr <- c(0,yrVec)

  for(y in 1:nyr){

    wt <- which(jdIndex == y)
    wi <- which(germYr < y & diedYr > (y-1))
    wl <- which(germYr >= y | diedYr <= (y-1))
    ni <- length(wi)
    if(ni == 0)next

    dataMat[wl,wt] <- sstate[wl,wt] <- NA

    sstate[cbind(wi,wt[1])] <- h[cbind(wi,wt[1])] <- 1
    dti <- diedIndex[wi,2]

    sy  <- sstate[wi,wt]
    if(!is.matrix(sy))sy <- matrix(sy,1)
    iy  <- matrix(c(1:length(wt)),ni,length(wt),byrow=T)*(sy*0 + 1)
    wff <- which(is.finite(sy),arr.ind=T)

    tt  <- 1

    for(t in wt[-c(1,length(wt))]){


      tt <- tt + 1
      dd <- jdAll[t]
      if(tt < (length(wt)-1)){
        iyt <- iy[,(tt+1):length(wt)]
        if(!is.matrix(iyt))iyt <- matrix(iyt,1)
        wnext <- apply(iyt,1,min,na.rm=T)
      }
      wnext[wnext == Inf] <- length(wt)
      snext <- sy[cbind(c(1:ni),wnext)]


      syt <- sy[,(tt-1):tt]
      if(!is.matrix(syt))syt <- matrix(syt,1)

      snow  <- apply(syt,1,min,na.rm=T)
      snext[!is.finite(snext)] <- snow[!is.finite(snext)]
      dydx  <- (snext - snow)/(wnext-tt)
      sy[,tt] <- snow + dydx
      sstate[wi,t] <- sy[,tt]
      h[wi,t] <- 98*(sstate[wi,t] - 1)/(nclass - 1) + 1
      h[wi[dti <= t],t] <- NA
      sstate[wi[dti <= t],t] <- NA

    }
    h[wi,wt[length(wt)]] <-  h[wi,wt[length(wt)]-1]
    sstate[wi,wt[length(wt)]] <-  sstate[wi,wt[length(wt)]-1]
  }

  sstate <- round(sstate,0)

  gmat <- matrix(yrBreak[germYr],n,nt)
  dmat <- matrix(diedIndex[,2],n,nt)
  hIndex <- xIndex <- which(gmat <= tmat & dmat > tmat,arr.ind=T)
  maxht <- diedIndex[hIndex[,1],2]
  maxht[maxht > nt] <- nt
  hIndex[,2] <- hIndex[,2] + 1

  www <- which(hIndex[,2] < maxht)
  hIndex <- hIndex[www,]
  xIndex <- xIndex[www,]

  www <- which(hIndex[,2] %in% yrBreak)
  if(length(www) > 0){
    hIndex <- hIndex[-www,]
    xIndex <- xIndex[-www,]
  }

 # h[sstate == 6] <- 96
  h <- inc2H(h,hIndex,xIndex)


  fix1 <- fix6 <- numeric(0)

  for(y in 1:nyr){

    wt  <- which(jdIndex == y)
    ntt <- length(wt)
    wt  <- wt[-c(1,(ntt-1),ntt)]
    nt  <- length(wt)

    dy <- dataMat[,wt]
    dy[,nt] <- 1
    dx <- dy*0 + 1

     first2 <- dy
     first2[first2 == 1] <- 1000
     first2[is.na(first2)] <- 1000
     first2 <- apply(first2,1,which.min)

     first1   <- apply(dy,1,which.min) 
     firstobs <- apply(dx,1,which.min)
     first1[firstobs < first1] <- 1

     last1 <- dy
     last1[,nt] <- NA
     last1[last1 > 1] <- NA
     last1 <- last1*matrix(wt,n,nt,byrow=T)
     last1[is.na(last1)] <- 0
     last1 <- apply(last1,1,which.max)

     last1[last1 > first2] <- first2[last1 > first2] - 1


     dy[,nt] <- nclass
     last6 <-  which(dy == nclass,arr.ind=T)
     tab6  <- table(last6[,1],last6[,2])
     tt6 <- matrix(as.numeric(colnames(tab6)),n,ncol(tab6),byrow=T)*tab6

     d5 <- dy
     d5[d5 == nclass] <- 0
     d5[,1] <- 1
     last5 <-  which(d5 > 0,arr.ind=T)
   
     tab5  <- table(last5[,1],last5[,2])
     tab5[,ncol(tab5)] <- 0

     tt5 <- apply(tab5*matrix(as.numeric(colnames(tab5)),n,ncol(tab5),byrow=T),1,max)

     dt6 <- tt6 - matrix(tt5,n,ncol(tt6))
     colnames(dt6) <- colnames(tab6)


     dt6[dt6 <= 0] <- 999
     first6 <- apply(dt6,1,which.min)
  #   first6[first6 == 1] <- ncol(dt6)

     first6 <- tt6[cbind( c(1:n),first6)]

    fix1 <- cbind(fix1,last1)
    fix6 <- cbind(fix6,first6+1)
}

  list(dataMat = dataMat, sstate = sstate, h = h, xIndex = xIndex, hIndex = hIndex, 
       fix1 = fix1, fix6 = fix6)
}

######################################################

inc2H <- function(hh,hIndex,xIndex){

  tiny <- .0001

  inc    <- hh[hIndex] - hh[xIndex]
  www    <- which(inc <= 0)

  windex <- hIndex[www,]
  if(length(windex) == 0)return(hh)

  inc[www] <- tiny

  imat <- matrix(NA,n,nt)
  imat[hIndex] <- inc

  for(y in 1:nyr){

    wt <- which(jdIndex == y)
    if(length(wt) == 0)next

    wi <- which(germIndex[,2] < min(wt))
    if(length(wi) == 0)next

    hy <- row2Mat(hh[wi,wt])
    im <- row2Mat(imat[wi,wt[-1]])

    hy[,-1] <- hy[,1] + t(apply(im,1,cumsum))
    hh[wi,wt] <- hy
  }

  hh[hh > (100 - tiny)] <- 100 - tiny

  hh
}

###############################

date2time <- function(date){  #julian days to time index

  tt <- findInterval(date,jdAll)
  tt[tt < 1] <- min(jdAll)
  tt[tt > nt] <- max(jdAll)
  tt
}
###############################

time2date <- function(time){

  time[time < 1] <- 1
  time[time > nt] <- nt
  dd <- jdAll[time]  
  dd
}
##############################
treatments <- function(siteName,cnames){

  nn <- length(cnames)

  if(length(siteName) > 1){
     damb   <- c('G01','G04','G07','S03','S04','S08','G1','G4','G7','S3','S4','S8')
     dplus3 <- c('G03','G05','G09','S02','S05','S09','G3','G5','G9','S2','S5','S9')
     dplus5 <- c('G02','G06','G08','S01','S06','S07','G2','G6','G8','S1','S6','S7')
     dcont  <- c('G10','G11','G12','S10','S11','S12')

     hamb   <- c('G02','G04','G07','S02','S05','S07')
     hplus3 <- c('G01','G06','G08','S03','S06','S08')
     hplus5 <- c('G03','G05','G09','S01','S04','S09')
     hcont  <- c('G10','G11','G12','S10','S11','S12')

     trt <- rep(character(0),nn)
     trt[siteName == 'DF' & cnames %in% damb]   <- 'amb'
     trt[siteName == 'DF' & cnames %in% dplus3] <- '+3'
     trt[siteName == 'DF' & cnames %in% dplus5] <- '+5'
     trt[siteName == 'DF' & cnames %in% dcont]  <- 'ctl'
     trt[siteName == 'HF' & cnames %in% hamb]   <- 'amb'
     trt[siteName == 'HF' & cnames %in% hplus3] <- '+3'
     trt[siteName == 'HF' & cnames %in% hplus5] <- '+5'
     trt[siteName == 'HF' & cnames %in% hcont]  <- 'ctl'
     return(trt)
  }

        if(siteName == 'DF'){
          amb   <- c('G01','G04','G07','S03','S04','S08','G1','G4','G7','S3','S4','S8')
          plus3 <- c('G03','G05','G09','S02','S05','S09','G3','G5','G9','S2','S5','S9')
          plus5 <- c('G02','G06','G08','S01','S06','S07','G2','G6','G8','S1','S6','S7')
          cont  <- c('G10','G11','G12','S10','S11','S12')
        }
        if(siteName == 'HF'){
          amb   <- c('G02','G04','G07','S02','S05','S07')
          plus3 <- c('G01','G06','G08','S03','S06','S08')
          plus5 <- c('G03','G05','G09','S01','S04','S09')
          cont  <- c('G10','G11','G12','S10','S11','S12')
        }

        maxr    <- max(c(length(amb),length(plus3),length(plus5),length(cont)))
        codemat <- matrix(NA,maxr,4)
        codemat[,1] <- amb
        codemat[,2] <- plus3
        codemat[,3] <- plus5
        codemat[,4] <- cont

        trt <- rep(NA,length(cnames))

        for(j in 1:maxr){
                wa <- grep(codemat[j,1],cnames)  
                w3 <- grep(codemat[j,2],cnames)  
                w5 <- grep(codemat[j,3],cnames)  
                wc <- grep(codemat[j,4],cnames) 

                if(length(wa) > 0) trt[wa] <- 'amb'
                if(length(w3) > 0) trt[w3] <- '+3'
                if(length(w5) > 0) trt[w5] <- '+5'
                if(length(wc) > 0) trt[wc] <- 'ctl'
        }
        trt
}

##############################

getTempAnomaly <- function(){

   tmp   <- inProcEnvData()
   AT_DF <- tmp$DFAT
   AT_HF <- tmp$HFAT

   dfCol <- grep('AT',colnames(AT_DF))
   hfCol <- grep('AT',colnames(AT_HF))

   yrs <- sort(unique(c(AT_DF[,'yeary'],AT_HF[,'yeary'])))
   nd  <- length(yrs)

   trtDF <- treatments('DF',colnames(AT_DF)[dfCol])
   trtHF <- treatments('HF',colnames(AT_HF)[hfCol])

   anomDF <- matrix(0,365,ncol(AT_DF))
   colnames(anomDF) <- colnames(AT_DF)
   anomHF <- matrix(0,365,ncol(AT_HF))
   colnames(anomHF) <- colnames(AT_HF)

   ddata <- AT_DF[,dfCol]
   hdata <- AT_HF[,hfCol]

   yrMeanDF <- yrMeanHF <- matrix(NA,365,3)

   amb <- apply(ddata[,trtDF %in% c('ctl','amb')],1,mean,na.rm=T)
   p3  <- apply(ddata[,trtDF == '+3'],1,mean,na.rm=T)
   p5  <- apply(ddata[,trtDF == '+5'],1,mean,na.rm=T)
   meanDF <- cbind(AT_DF[,1:3],amb,p3,p5)                #across all treatments

   amb <- apply(hdata[,trtHF %in% c('ctl','amb')],1,mean,na.rm=T)
   p3  <- apply(hdata[,trtHF == '+3'],1,mean,na.rm=T)
   p5  <- apply(hdata[,trtHF == '+5'],1,mean,na.rm=T)
   meanHF <- cbind(AT_HF[,1:3],amb,p3,p5)

   sumD <- matrix(0,365,3)
   colnames(sumD) <- c('amb','p3','p5')
   sumH <- nD <- nH <- sumD

   annualAve <- matrix(NA,length(years),6)
   colnames(annualAve) <- t(outer(c('DF','HF'),colnames(sumD),paste,sep='-'))
   rownames(annualAve) <- years

   for(y in 1:length(years)){

      wd <- which(AT_DF[,'yeary'] == yrs[y])
      dfRow <- AT_DF[wd,'ydate'] - (y - 1)*365 + 1

      sumD[dfRow,1] <- sumD[dfRow,1] + apply(ddata[wd,trtDF %in% c('ctl','amb')],1,sum,na.rm=T)
      nD[dfRow,1]   <- nD[dfRow,1] + apply(ddata[wd,trtDF %in% c('ctl','amb')]*0+1,1,sum,na.rm=T)

      sumD[dfRow,2] <- sumD[dfRow,2] + apply(ddata[wd,trtDF == '+3'],1,sum,na.rm=T)
      nD[dfRow,2]   <- nD[dfRow,2] + apply(ddata[wd,trtDF == '+3']*0+1,1,sum,na.rm=T)

      sumD[dfRow,3] <- sumD[dfRow,3] + apply(ddata[wd,trtDF == '+5'],1,sum,na.rm=T)
      nD[dfRow,3]   <- nD[dfRow,3] + apply(ddata[wd,trtDF == '+5']*0+1,1,sum,na.rm=T)


      wd <- which(AT_HF[,'yeary'] == yrs[y])
      dfRow <- AT_HF[wd,'ydate'] - (y - 1)*365 + 1

      sumH[dfRow,1] <- sumH[dfRow,1] + apply(hdata[wd,trtHF %in% c('ctl','amb')],1,sum,na.rm=T)
      nH[dfRow,1]   <- nH[dfRow,1] + apply(hdata[wd,trtHF %in% c('ctl','amb')]*0+1,1,sum,na.rm=T)

      sumH[dfRow,2] <- sumH[dfRow,2] + apply(hdata[wd,trtHF == '+3'],1,sum,na.rm=T)
      nH[dfRow,2]   <- nH[dfRow,2] + apply(hdata[wd,trtHF == '+3']*0+1,1,sum,na.rm=T)

      sumH[dfRow,3] <- sumH[dfRow,3] + apply(hdata[wd,trtHF == '+5'],1,sum,na.rm=T)
      nH[dfRow,3]   <- nH[dfRow,3] + apply(hdata[wd,trtHF == '+5']*0+1,1,sum,na.rm=T)

      md <- sumD/nD
      mh <- sumH/nH

      annualAve[y,1:3] <- apply(md,2,mean,na.rm=T)
      annualAve[y,4:6] <- apply(mh,2,mean,na.rm=T)
   }

   meanAllDF <- sumD/nD
   meanAllHF <- sumH/nH

   anomDF <- ddata*0
   anomHF <- hdata*0

   for(y in 1:length(years)){

      wd <- which(AT_DF[,'yeary'] == yrs[y])
      dfRow <- AT_DF[wd,'ydate'] - (y - 1)*365 + 1

      anomDF[wd,] <- ddata[wd,] - matrix(meanAllDF[dfRow,'amb'],nrow=length(dfRow),ncol=ncol(ddata))

      wd <- which(AT_HF[,'yeary'] == yrs[y])
      dfRow <- AT_HF[wd,'ydate'] - (y - 1)*365 + 1

      anomHF[wd,] <- hdata[wd,] - matrix(meanAllHF[dfRow,'amb'],nrow=length(dfRow),ncol=ncol(hdata))
  }

  anomDF <- cbind(AT_DF[,1:3],anomDF)
  anomHF <- cbind(AT_HF[,1:3],anomHF)

  list(meanyrDF = meanAllDF, meanyrHF = meanAllHF, anomDF = anomDF, anomHF = anomHF,
       annualAve = annualAve,dfTemp = AT_DF, hfTemp = AT_HF)
}

   	 
###############################

makeXmat <- function(){
	
  xmat  <- matrix(0,nrow(xIndex),nx)
  colnames(xmat) <- xnames
  xlist <- numeric(0)

  for(k in 1:nx){
   	
    xx <- tmat*0
   	   
    if(xnames[k]== 'seed'){
   	 seed <- seedOrigin()$seed
   	 xx <- matrix(seed,n,nt)
    }
    if(xnames[k] == 'site'){
      ss <- match(specData[,'site'],siteCode) - 1
      xx <- matrix(ss,n,nt)
    }
    if(xnames[k] == 'gap')xx <- matrix(specData[,'gap'],n,nt)

    if(xnames[k] == 'h')  xx <- h
    if(xnames[k] == 'h2') xx <- h^2
    if(xnames[k] == 'int')xx <- xx + 1

    if(xnames[k] == 'chill'){

      xx <- matrix(NA,n,nt)

      tmp   <- inProcEnvData()
 #     hotD  <- tmp$hotDF
 #     hotH  <- tmp$hotHF
      coldDa <- tmp$coldDFa
      coldHa <- tmp$coldHFa
      coldDs <- tmp$coldDFs
      coldHs <- tmp$coldHFs

      colnames(coldDa) <- unlist(strsplit(colnames(coldDa),'AT'))
      colnames(coldHa) <- unlist(strsplit(colnames(coldHa),'AT'))
      colnames(coldDs) <- unlist(strsplit(colnames(coldDs),'ST'))
      colnames(coldHs) <- unlist(strsplit(colnames(coldHs),'ST'))

      yrMatch <- match(years,rownames(coldDa))
      if(is.na(yrMatch)[1]){
         coldDa <- rbind(coldDa[1,],coldDa)
         rownames(coldDa)[1] <- years[1]
      }
      yrMatch <- match(years,rownames(coldHa))
      if(is.na(yrMatch)[1]){
         coldHa <- rbind(coldHa[1,],coldHa)
         rownames(coldHa)[1] <- years[1]
      }
      yrMatch <- match(years,rownames(coldDs))
      if(is.na(yrMatch)[1]){
         coldDs <- rbind(coldDs[1,],coldDs)
         rownames(coldDs)[1] <- years[1]
      }
      yrMatch <- match(years,rownames(coldHs))
      if(is.na(yrMatch)[1]){
         coldHs <- rbind(coldHs[1,],coldHs)
         rownames(coldHs)[1] <- years[1]
      }

      for(i in 1:n){

        coldMata <- coldDa
        coldMats <- coldDs
        if(specData[i,'site'] == 'HF'){
            coldMata <- coldHa
            coldMats <- coldHs
        }
       
        icca <- match(cham[i],colnames(coldMata))
        iccs <- match(cham[i],colnames(coldMats))

        xx[i,] <- round(apply(rbind(coldMata[jdIndex,icca],coldMats[jdIndex,iccs]),2,mean),1)
     }
   }

    if(xnames[k] %in% c('AT','SM','ATXSM')){

      tmp   <- inProcEnvData()
      AT_DF <- tmp$DFAT
      AT_HF <- tmp$HFAT
      SM_DF <- tmp$DFSM
      SM_HF <- tmp$HFSM

      if(xnames[k] == 'AT'){
        DF <- AT_DF
        HF <- AT_HF
        mdf <- match(paste(cham,'AT',sep=''),colnames(DF))
        mhf <- match(paste(cham,'AT',sep=''),colnames(HF))
      }
      if(xnames[k] == 'SM'){
        DF <- SM_DF
        HF <- SM_HF
        mdf <- match(paste(cham,'SM',sep=''),colnames(DF))
        mhf <- match(paste(cham,'SM',sep=''),colnames(HF))
      }

      if(xnames[k] == 'ATXSM'){
            mt  <- match(SM_DF[,'ydate'],AT_DF[,'ydate'])
            wm     <- which(is.na(mt))
            if(length(wm) > 0){
              for(j in wm)mt[j] <- max(mt[is.finite(mt) & mt < j])
            }
            DF <- AT_DF[mt,]
            dat <- grep('AT',colnames(DF))
            DF[,dat] <- (DF[,dat] - mean(DF[,dat],na.rm=T))*
        	               (SM_DF[,dat] - mean(SM_DF[,dat],na.rm=T))
        	            
            mt  <- match(SM_HF[,'ydate'],AT_HF[,'ydate'])
            wm  <- which(is.na(mt))
            if(length(wm) > 0){
              for(j in wm)mt[j] <- max(mt[is.finite(mt) & mt < j])
            }
            HF  <- AT_HF[mt,]
            dat <- grep('AT',colnames(HF))
            HF[,dat] <- (HF[,dat] - mean(HF[,dat],na.rm=T))*
        	               (SM_HF[,dat] - mean(SM_HF[,dat],na.rm=T))   
       }
        
       mt <- findInterval(jdAll,DF[,'ydate'])
       DF <- DF[mt,]

       mt <- findInterval(jdAll,HF[,'ydate'])
       HF <- HF[mt,]

       ad <- t(DF[,mdf])
       ah <- t(HF[,mhf])        
        
       xx[site == 'DF',] <- ad[site == 'DF',]
       xx[site == 'HF',] <- ah[site == 'HF',]
     }
      
     if(xnames[k] == 'ATXseed')xx <- xx*0
     if(xnames[k] == 'ATXgap') xx <- xx*0

     xtmp <- xx[xIndex]
     if(is.na(sum(xtmp))){
        wff <- xIndex[!is.finite(xx[xIndex]),]
        if(!is.matrix(wff))wff <- matrix(wff,1,2)
        stop( paste(' ',xnames[k],' not finite, i = ',wff[,1],', t =', wff[,2],sep=' ') )
     }             
     xmat[,k] <- xx[xIndex]
     xlist <- append(xlist,list(xx))
  }


  wf  <- which(xnames %in% factors)
  wwr <- numeric(0)
  if(length(wf) > 0){
    for(kk in wf){
        stot <- table(xmat[,kk])
        if(length(stot) < 2){
          wwr <- c(wwr,kk)
          next
        }
        if(min(stot) < 10){
          wwr <- c(wwr,kk)
          next
        }
     }
   }
   if(length(wwr) > 0){
     if('seed' %in% xnames[wwr] & 'ATXseed' %in% xnames)wwr <- c(wwr,which(xnames == 'ATXseed'))
     if('gap' %in% xnames[wwr] & 'ATXgap' %in% xnames)wwr <- c(wwr,which(xnames == 'ATXgap'))
     xm   <- xmat[,-wwr]
     xnew <- numeric(0)
     for(k in c(1:ncol(xmat))[-wwr]){
       xnew <- append(xnew,list(xlist[[k]]))
     }
     xlist <- xnew
     xmat  <- xm
  }


  xnames <- colnames(xmat)
  nx     <- length(xnames)
    
  #center variables
    
  meanX  <- apply(xmat,2,mean,na.rm=T)
  rangeX <- apply(xmat,2,range,na.rm=T)
    
  for(k in 1:length(xnames)){
    	if(xnames[k] %in% c('AT','SM','chill')){
    	  xmat[,k] <- xmat[,k] - meanX[k]
    	  xx <- tmat*0
    	  xx[xIndex] <- xmat[,k]
         xlist[[k]] <- xx
      }

    	if(xnames[k] == 'ATXSM'){
    		xmat[,k] <- xmat[,'AT']*xmat[,'SM']
    		xx <- tmat*0
    	   xx[xIndex] <- xmat[,k]
         xlist[[k]] <- xx
      }
    	if(xnames[k] == 'ATXchill'){
    		xmat[,k] <- xmat[,'AT']*xmat[,'chill']
    		xx <- tmat*0
    	   xx[xIndex] <- xmat[,k]
         xlist[[k]] <- xx
      }
    	if(xnames[k] == 'ATXseed'){
    	  xmat[,k] <- xmat[,'AT']*xmat[,'seed']
    	  xx <- tmat*0
    	  xx[xIndex] <- xmat[,k]
         xlist[[k]] <- xx
       }
    	if(xnames[k] == 'ATXgap'){
    	  xmat[,k] <- xmat[,'AT']*xmat[,'gap']
    	  xx <- tmat*0
    	  xx[xIndex] <- xmat[,k]
         xlist[[k]] <- xx
       }
    }
            
    
    x3d <- array(NA,dim=c(n,nt,ncol(xmat)))
    for(j in 1:length(xlist)){
    	x3d[,,j] <- xlist[[j]]
    }
    dimnames(x3d)[[2]] <- times
    dimnames(x3d)[[3]] <- xnames
    rm(xlist)


  list(xmat = xmat, x3d = x3d, meanX = meanX, rangeX = rangeX)
}
	
##############################	 
initB <- function(){  #initialize stateSpace pars
	
   xtmp <- xmat[-missRow,]
   htmp <- log(h[hIndex[-missRow,]]) - log(h[xIndex[-missRow,]])
   V    <- solve(dt/sigma*crossprod(xtmp))
   v    <- 1/sigma*crossprod(xtmp,htmp)
   
   V%*%v
}

#################################

makeY <- function(h0,h1,hmin){

  wh <- which(is.finite(h0) & h0 > hmin & h0 < (100 - hmin) &
              is.finite(h1) & h1 > hmin & h1 < (100 - hmin))
  y  <- rep(NA,length(h0))
  y[wh] <- log( (h1[wh] - h0[wh])/dt/(1 - h0[wh]/100) )

  y
}

#################################
updateB <- function(){
	
   minGro <- log(.05)    #only fit when active growth
 #  bord <- sample(1:length(xnames))

   y <- makeY(h[xIndex],h[hIndex],.1)

   wy   <- which(is.finite(y) & y > minGro)
   xtmp <- xmat[wy,]
   y    <- y[wy]
	
   V    <- solve(1/sigma*crossprod(xtmp) + priorIVB)
   v    <- 1/sigma*crossprod(xtmp,y) + priorIVB%*%priorB
   mu   <- V%*%v

  lo <- loBA

  ws <- grep('ATXseed',xnames)
  if(length(ws) > 0){
    lo[ws] <- -bg['AT',1]
  }

  wc <- grep('ATXchill',xnames)
  if(length(wc) > 0){
    lo[wc] <- -bg['AT',1]
  }

  bg <- tnorm.mvt(bg,mu,V,lo,hiBA)

  if(length(ws) > 0)bg <- tnorm.mvt(bg,mu,V,lo,hiBA,whichSample=ws)
  if(length(wc) > 0)bg <- tnorm.mvt(bg,mu,V,lo,hiBA,whichSample=wc)

   u1 <- s1 + nrow(xtmp)/2
   u2 <- s2 + .5*crossprod(y - xtmp%*%bg)

   mm <- 1 - pgamma(1/maxS,u1,u2)
   z  <- runif(1,0,mm)
   ss <- 1/qgamma(1 - z,u1,u2)
   if(ss == 0)ss <- maxS
   
   list(bg = bg,sigma = ss)
}   



#################################
updatemissX <- function(xmat){
	
  if(length(missRow) == 0)return(xmat)
	
  hnow  <- h[xIndex[missRow,]]
  hnext <- h[hIndex[missRow,]]
 # hlast <- h[xIndex[missRow,]]

  xtmp <- xmat
  
  rangeLim <- apply(xmat,2,range,na.rm=T)
  
  for(j in 1:length(missRow)){
  	
  	  mj <- missRow[j]
  	  wm <- missX[missX[,1] == mj,2]
          lo <- rangeLim[1,wm]
          hi <- rangeLim[2,wm]
          lo[lo < (xtmp[mj-1,wm] - 5)] <- xtmp[mj-1,wm] - 5
          hi[hi < (xtmp[mj-1,wm] + 5)] <- xtmp[mj-1,wm] + 5

  	  if(length(wm) == 1){
  	  	if(xnames[wm] == 'ATXSM'){
  	  	  xtmp[mj,'ATXSM'] <- xtmp[mj,'AT']*xtmp[mj,'SM']
  	  	  next
  	  	}
  	  	if(xnames[wm] == 'seed'){
  	  	  xtmp[mj,'seed'] <- 0
  	  	  next
  	  	}
  	  }
  	  
  	  b1 <- matrix(bg[wm,],ncol=1)       #missing
  	  b2 <- matrix(bg[-wm,],ncol=1)      #not missing
  	  q  <- log(hnext[j]) - A*log(hnow[j])
  	  xj <- q - matrix(xmat[mj,-wm],nrow=1)%*%b2*dt
          if(length(wm) > 1){
            zz <- nearPD(crossprod(t(b1)))$mat  	  
  	    z  <- t(b1)%*%solve(zz)
  	    mu <- xj%*%z/dt
            xtmp[mj,wm] <- as.numeric(tnorm.mvt(mu,mu,zz*sigma*dt,rangeLim[1,wm],rangeLim[2,wm]))
          }
          if(length(wm) == 1){
            mu <- xj/b1/dt
            zz <- b1^2
            xtmp[mj,wm] <- tnorm(1,rangeLim[1,wm],rangeLim[2,wm],mu,sqrt(zz*sigma*dt))
          }  
            

          if('ATXSM' %in% xnames)  xtmp[mj,'ATXSM']   <- xtmp[mj,'AT']*xtmp[mj,'SM']
	  if('ATXseed' %in% xnames)xtmp[mj,'ATXseed'] <- xtmp[mj,'AT']*xtmp[mj,'seed']
	}
	xtmp
}	 

  


############################################

makePriorB <- function(ir,sr,breaks,h){

  n  <- length(breaks)
  mm <- matrix(cbind(seq(ir[1],ir[2],length.out=n),seq(sr[1],sr[2],length.out=n)),n,2)
  cmat <- breaks2pars(mm,breaks)
  tmp  <- pars2p(cmat,h)
  cmat
}


################################

updateObsError <- function(){
	
   omat <- matrix(0,nclass,nclass)
	
   for(j in 1:nclass){
	wj <- which(is.finite(dataMat) & sstate == j)
	oj <- table(dataMat[wj])
	omat[match(names(oj),c(1:nclass)),j] <- oj
	di <- rgamma(nclass,shape=(omat[,j] + obsPrior[,j]),scale=1)
	omat[,j] <- di/sum(di)
   }
   omat
}

######################################

procModel <- function(h0,h1,mu,dt,ss,tiny){

  H <- h0/100
 # H[H >= (1 - tiny)] <- 1 - tiny

  y <- log( (h1 - h0)/dt/(1 - H) )
  dnorm(y,mu,sqrt(ss),log=T)
}

##########################
statesAtT <- function(hhh,sss,tvec,tiny){

  ww <- which(hIndex[,2] == tvec[2])
  if(length(ww) == 0)return(list(hs = numeric(0), ss = numeric(0), wi = numeric(0)))

  hs <- hhh[,tvec]
  ss <- sss[,tvec]

  wi <- sort(unique(hIndex[ww,1]))

  hs[hs[,3] > (100 - tiny),3] <- 100 - tiny

  wl <- which(hs[wi,3] < hs[wi,2] & hs[wi,2] >= hs[wi,1])
  if(length(wl) > 0){
     hs[wi[wl],3] <- hs[wi[wl],2] + tiny
  }

  wl <- which(hs[wi,3] < hs[wi,2])
  if(length(wl) > 0)hs[wi[wl],2] <- hs[wi[wl],3] - tiny

  wl <- which(hs[wi,2] < hs[wi,1])
  if(length(wl) > 0)hs[wi[wl],1] <- hs[wi[wl],2] - tiny


  d0 <- ss[wi,2] - ss[wi,1]
  d1 <- ss[wi,3] - ss[wi,2]

  wneg0 <- which(d0 < 0)
  if(length(wneg0) > 0)ss[wi[wneg0],2] <- ss[wi[wneg0],1]

  wneg <- which(d0 > 0 & d1 < 0)
  if(length(wneg) > 0){
     tmp <- apply(ss[,-2],1,max,na.rm=T)
     ss[wi[wneg],2] <- tmp[wi[wneg]]
  }

  d1 <- ss[wi,3] - ss[wi,2]
  wneg1 <- which(d1 < 0)
  if(length(wneg1) > 0)ss[wi[wneg1],3] <- ss[wi[wneg1],2]

  list(hs = hs, ss = ss, wi = wi)
}
########################

propDis <- function(dat,sn,hh,wi){   #propose discrete states

  ni  <- length(wi)
  swi <- sn
  if(ni == 1)swi <- matrix(swi,1,3)

  swi[is.na(swi[,1]),1] <- swi[is.na(swi[,1]),2]
  swi[is.na(swi[,3]),3] <- swi[is.na(swi[,3]),2]

  wdd <- which( (swi[,2] - dat) > 1)
  swi[wdd,1] <- swi[wdd,2] - 1

  d0 <- swi[,2] - swi[,1]
  d1 <- swi[,3] - swi[,2]

  ps <- rep(0,ni)
  ps[d0 == 0 & d1 > 0] <- 1
  ps[d0 > 0 & d1 == 0] <- -1

  wneg <- which(d0 < 0 & d1 > 0)
  if(length(wneg) > 0)swi[wneg,1] <- swi[wneg,2]
  
  wneg <- which(d1 < 0)
  if(length(wneg) > 0)swi[wneg,3] <- swi[wneg,2]
 
  snow <- swi
  snew <- swi[,2] + ps
  snew[snew < snow[,1]] <- snow[snew < snow[,1],1]
  snew[snew > snow[,3]] <- snow[snew > snow[,3],3]
  swi[,2] <- snew

  pnow <- multiLogitStates(bq, snow[,2],hh,nclass)
  pnew <- multiLogitStates(bq, snew,hh,nclass)

  wdat <- which(is.finite(dat))

  psnow <- log(obsErr[cbind(snow[wdat,2],dat[wdat])])
  psnew <- log(obsErr[cbind(snew[wdat],dat[wdat])])

  inow <- which(psnow == -Inf)
  inew <- which(psnew == -Inf)

  if(length(inow) > 0){
      snow[wdat[inow],2] <- dat[wdat[inow]]
      psnow <- log(obsErr[cbind(snow[wdat,2],dat[wdat])])
  }

  if(length(inew) > 0){
    snew[wdat[inew]] <- dat[wdat[inew]]
    psnew <- log(obsErr[cbind(snew[wdat],dat[wdat])])
  }

  pnow[wdat] <- pnow[wdat] + psnow
  pnew[wdat] <- pnew[wdat] + psnew

  a <- exp(pnew - pnow)
  z <- runif(ni,0,1)
  snow[z < a,2] <- snew[z < a]

  list(ss = snow)
}

plotObsError <- function(){

    par(mfcol=c(1,1),bty='n',cex.axis=1.2,cex.lab=1.2)
 
    colE <-  colorRampPalette(c("black","white"))
    cols <- colE(5)
	
    oe <- -log10(oerror)
    
    oe[oe > 9] <- 9 
    oe <- oe + .5

    xlimit <- ylimit <- c(-.5,nclass+.5)

    plot(-10,-10,xlim=xlimit,ylim=ylimit,xlab=' ',ylab='Observed State',xaxt='n',yaxt='n')
    axis(side=3,at=(c(1:nclass)-.5),labels=c(1:nclass),tick=F,line=-2)
    axis(side=2,at=(c(1:nclass)-.5),labels=c(1:nclass),tick=F,line=-2)
    title('Inferred State')

    for(k in nclass:1){
      for(j in 1:nclass){
        polygon(c(j-1,j,j,j-1),c(k-1,k-1,k,k),border='grey',col=cols[round(oe[j,k],0)])
      }
      oe3 <- signif(oerror[k,k],4)
      text(k-.5,k-.5,oe3,col='white')
    }
}
#######################################################
propCon <- function(hs,wi,tiny){  #propose continuous

  hnew <- hs[wi,]
  if(length(wi) == 1)hnew <- matrix(hnew,1,3)

  wmm <- which(hnew[,3] > (100 - tiny))
  hnew[wmm,3]  <- 100 - tiny

  hnew[is.na(hnew[,1]),1] <- hnew[is.na(hnew[,1]),2]
  hnew[is.na(hnew[,3]),3] <- hnew[is.na(hnew[,3]),2]

  hnew[hnew[,2] > hnew[,3],2]  <- hnew[hnew[,2] > hnew[,3],3] - tiny
  hnew[hnew[,1] > hnew[,2],1]  <- hnew[hnew[,1] > hnew[,2],2] - tiny

  hnew[,2] <- tnorm(length(wi),hnew[,1],hnew[,3],hnew[,2],rexp(length(wi),2))

  hs[wi,] <- hnew
  hs
}
###########################################################
getDD <- function(threshold=0){

  xat   <- x3d[,,'AT']
  xat[xat != 0] <- xat[xat != 0] + meanX['AT']
  ddmat <- matrix(NA,n,nt)

  for(y in 1:nyr){

    wt <- which(jdIndex == y)
    if(length(wt) == 0)next

    xt <- xat[,wt]
    xt[xt < 0] <- 0

    dd <- t(apply(xt,1,cumsum))
    if(max(dd) < 50)next
    ddmat[,wt] <- dd
  }
  ddmat
}


######################

updateStates <- function(){ #for h_t+dt = h_t + exp(x%*%b)(1 - h_t/100)*dt + epsilon*sqrt(dt), sample increment
	
  predh  <- predObs <- matrix(1,n,nt)
  counth <- counts <- h*0

  dev <- 0
  tiny <- .001
  teenytiny <- .001*tiny

  hhh <- h

  for(y in 1:nyr){

    wt  <- which(jdIndex == y)
    ntt <- length(wt)
  #  wt  <- wt[-c(1,(ntt-1),ntt)]

    wt  <- wt[-c(1,ntt)]

    kk <- 0

    for(t in wt){

      kk <- kk + 1

      dd  <- jdAll[t]
      jdy <- dd
      if(y > 1)jdy <- jdy - yrVec[y-1]
  
      sstate[ddmat[,t] < 1,t] <- 1

      tvec <- c((t-1):(t+1))

      if(t %in% yrBreak)hhh[,t] <- 100

      temp <- x3d[,t,'AT']
      threshold <- rnorm(n,threshMu-meanX['AT'],threshSd)

      tmp <- statesAtT(hhh,sstate,tvec,tiny*10)
      wi  <- tmp$wi
      ni  <- length(wi)
      if(ni == 0)next
      hs  <- tmp$hs
      ss  <- tmp$ss

      if(kk == 1){
        hs[wi,1] <- 1
        hs[wi,2] <- 1
        hs[wi,3] <- 1 + tiny
        ss[wi,1:2] <- 1
      }

      ww  <- which(is.na(hs[wi,1]))
      if(length(ww) > 0)hs[wi[ww],1] <- hs[wi[ww],2] - .01

      tmp <- propDis(dataMat[wi,t],ss[wi,],hs[wi,2],wi)$ss
      ss[wi,]  <- tmp
 
      wf1 <- which(fix1[,y] > kk)
      if(length(wf1) > 0)ss[wf1,1:2] <- 1

      wf6 <- which(kk > fix6[,y])
      if(length(wf6) > 0){
         ss[wf6,2:3] <- 6
         hs[wf6,2:3] <- 100
      }

      sstate[wi,t] <- ss[wi,2]

      hnew <- propCon(hs,wi,tiny)
      hs[wi,1] <- hnew[wi,1]
      hs[wi,3] <- hnew[wi,3]

      wcold <- which(temp < threshold)      #below threshold temperature
      if(length(wcold) > 0 & kk > 1){
         hs[wcold,2] <- hnew[wcold,2] <- hs[wcold,1] + teenytiny
         ss[wcold,2] <- ss[wcold,1]
      }

      xmat0   <- x322(wi,t-1,hs[wi,1])
      xmatNow <- x322(wi,t,hs[wi,2])
      xmatNew <- x322(wi,t,hnew[wi,2])

      enow0 <- enow1 <- enew1 <- wi*0
      for(k in 1:length(xnames)){
         enow0 <- enow0 + xmat0[,k]*bg[k,1]
         enow1 <- enow1 + xmatNow[,k]*bg[k,1]
         enew1 <- enew1 + xmatNew[,k]*bg[k,1]
      }

      pnow  <- procModel(hs[wi,1],hs[wi,2],enow0,dt,sigma,tiny/100)
      pnew  <- procModel(hs[wi,1],hnew[wi,2],enow0,dt,sigma,tiny/100)

      wnext <- which(is.finite(hs[wi,3]))
      pnow[wnext] <- pnow[wnext] + procModel(hs[wi[wnext],2],hs[wi[wnext],3],enow1[wnext],dt,sigma,tiny)
      pnew[wnext] <- pnew[wnext] + procModel(hnew[wi[wnext],2],hs[wi[wnext],3],enew1[wnext],dt,sigma,tiny)

     pnow <- pnow + multiLogitStates(bq, ss[wi,2],hs[wi,2],nclass)
     pnew <- pnew + multiLogitStates(bq, ss[wi,2],hnew[wi,2],nclass)

     pnow[is.na(pnow)] <- -Inf
     pnew[is.na(pnew)] <- -Inf

     a       <- exp(pnew - pnow)
     z       <- runif(ni,0,1)
     hh      <- hs[wi,]
     hn      <- hnew[wi,]
     if(!is.matrix(hh))hh <- matrix(hh,1,3)
     if(!is.matrix(hn))hn <- matrix(hn,1,3)

     wh      <- which(z < a)
     hh[wh,2]  <- hn[wh,2]
     
     hhh[wi,tvec] <- hh

     counth[wi[wh],t] <- counth[wi[wh],t] + 1
     counts[wi[wh],t] <- counts[wi[wh],t] + 1

    #predict observations
      ddh <- exp(enow0 + rnorm(ni,0,sqrt(sigma)))*dt*(1 - predh[wi,t-1]/100)
      ddh[wi %in% wcold] <- 0

      predh[wi,t] <- predh[wi,t-1] + ddh
      tmpNow <- x3d[,t,xnames == 'AT'] + meanX['AT']
      wplus <- which(tmpNow < 0)
      if(length(wplus) > 0)predh[wplus,t] <- 1
      probs <- pars2p(bq,predh[wi,t])

      newO  <- apply(myrmultinom(1,probs),1,which.max)
      predObs[wi,t] <- newO
       
      dev <- dev - 2*sum(multiLogitStates(bq, dataMat[wi,t], predh[wi,t],nclass),na.rm=T)
   }
  }

  list(h = hhh, s = sstate,predObs = predObs, counts = counts, counth = counth, dev = dev)
}
   

x322 <- function(ii,tt,hvalues){  #individuals ii, a time tt, hvalues h[ii,tt]

  xnew <- x3d[ii,tt,]
  ni   <- length(ii)

  if(ni == 1){
    xnew <- matrix(xnew,1,nx)
    colnames(xnew) <- xnames
  }
  if('h' %in% xnames)xnew[,'h']   <- hvalues
  if('h2' %in% xnames)xnew[,'h2'] <- hvalues^2

  xnew
}
##################################################



date2Phen <- function(){   #first obs dates and DD for each stage

  ymat <- numeric(0)
  gmat <- matrix(germYr,n,nt)

  for(y in 2:nyr){

    wt <- which(jdIndex == y)
    if(length(wt) == 0)next

    for(k in 2:nclass){

     kmat  <- dataMat[,wt]
     kmat[kmat < k] <- NA
     kmat[gmat[,wt] >= y] <- NA
     kdate <- apply((kmat*0 + 1)*tmat[,wt],1,min,na.rm=T)
     wi <- which(!is.finite(kdate))
     kdate[wi] <- min(kdate,na.rm=T)

     jdate <- jdAll[kdate] - yrVec[y] + 365
     ddate <- ddmat[cbind(c(1:n),kdate)]
     jdate[wi] <- ddate[wi] <- NA
     yk <- cbind(jdate,round(ddate,0))
     colnames(yk) <- paste(c('JD','DD'),years[y],'k',k,sep='-')
     ymat <- cbind(ymat,yk)
   }
  }

  wf <- which(is.finite(ymat),arr.ind=T)[,1]
  if(length(wf) == 0)return(numeric(0))

  colF <-  colorRampPalette(c("blue","orange","red"))
  colVals <- colF(nclass)

  plot(c(0,0),c(0,0),xlim=c(20,250),ylim=c(0,1500))

  for(k in 2:nclass){

    wc <- grep(paste('k',k,sep='-'),colnames(ymat))
    wj <- grep('JD',colnames(ymat))
    wd <- grep('DD',colnames(ymat))

    jd <- intersect(wc,wj)
    dd <- intersect(wc,wd)

    points(as.vector(ymat[,jd]),as.vector(ymat[,dd]),col=colVals[k])
    
  }

  ymat
}

getTraits <- function(spec){

  traitNames <- c('xylem','succession')
  xyName <- c('ring-porous','diffuse-porous','tracheid')
  suName <- c('early','int','late')

  specs <- c("acru","acsa","acun","beal","bepo","beun","fram","list","litu","magr",
             "nysy","pipa","pire","pist","pita","piun","prse","qual","quru","quun")

  trMat <- matrix(NA,length(specs),length(traitNames))
  rownames(trMat) <- specs
  colnames(trMat) <- traitNames

  trMat[specs == 'acru',] <- c(2,2)
  trMat[specs == 'acsa',] <- c(2,3)
  trMat[specs == 'acun',] <- c(2,2)  
  trMat[specs == 'beal',] <- c(2,3)
  trMat[specs == 'bepo',] <- c(2,1)
  trMat[specs == 'beun',] <- c(2,2)
  trMat[specs == 'fram',] <- c(1,2)
  trMat[specs == 'list',] <- c(2,2)
  trMat[specs == 'litu',] <- c(2,2)
  trMat[specs == 'magr',] <- c(2,3)
  trMat[specs == 'nysy',] <- c(2,3)
  trMat[specs == 'pipa',] <- c(3,1)
  trMat[specs == 'pire',] <- c(3,1)
  trMat[specs == 'pist',] <- c(3,1)
  trMat[specs == 'pita',] <- c(3,1)
  trMat[specs == 'piun',] <- c(3,1)
  trMat[specs == 'prse',] <- c(2,2)
  trMat[specs == 'qual',] <- c(1,2)
  trMat[specs == 'quru',] <- c(1,2)
  trMat[specs == 'quun',] <- c(1,2)

  xylem      <- xyName[trMat[,'xylem']]
  succession <- suName[trMat[,'succession']]

  traitData <- as.data.frame(cbind(specs,xylem,succession))
  list(trMat = trMat, traitTable = traitData)
}

lmFit <- function(){

  wk <- grep('k-3',colnames(firstDates))
  wj <- grep('JD',colnames(firstDates))
  wd <- grep('DD',colnames(firstDates))

  wj <- intersect(wk,wj)
  wd <- intersect(wk,wd)

  aveTemp <- firstDates[,wd]/firstDates[,wj]
  www <- which(is.finite(aveTemp) & aveTemp != 0,arr.ind=T)

  y  <- firstDates[,wj][www]
  ave <- aveTemp[www]

  nyy <- length(wj)
  xxx <- numeric(0)
  xn  <- xnames

  modelForm <- '~'
  FIRST <- T

#modelForm <- '~ warm + site + seed + warm:seed'

  for(k in 1:length(xn)){
    
    kn <- xn[k]
    if(kn == 'int')xx <- rep(1,length(ave))
    if(kn == 'AT')xx <- ave
    if(kn == 'seed')   xx <- seed[www[,1]]
    if(kn == 'ATXseed')xx <- seed[www[,1]]*ave
    if(kn == 'gap')    xx <- specData[www[,1],'gap']
    xxx <- cbind(xxx,xx)
    if(kn == 'int')next
    if(FIRST) modelForm <- paste(modelForm,kn,sep=' ')
    if(!FIRST)modelForm <- paste(modelForm,' + ',kn,sep='')
    FIRST <- F
  }
  colnames(xxx) <- xn
  modelForm <- gsub('X',':',modelForm)
  modelForm <- as.formula(paste('y',modelForm,sep=' '))

  lfit <- lm(modelForm,data=as.data.frame(xxx))
  tmp  <- summary(lfit)

  list(coeffs = tmp$coefficients, n = length(y),rsq = tmp$r.squared)
}




hFix <- function(h){

    wh <- which(is.na(h[hIndex]))

    if(length(wh) > 0){
      hi <- matrix(hIndex[wh,],ncol=2)
      h1 <- hi
      h1[,2] <- hi[,2] - 1
      hh <- h[h1] + .000001
      wi <- which(is.na(hh))
      if(length(wi) > 0){
         h1[,2] <- hi[,2] + 1
         hh[wi] <- h[h1] - .000001
      }
      h[cbind(hIndex[wh,1],hIndex[wh,2])] <- hh
      message('h fix')
    }
  h
}



propSH <- function(){
	
	tmp <- initialStatus()
	ss  <- tmp$sstate
	hh  <- ss*0 + 98
	hn  <- h
	
	sb  <- ss[hIndex]
	lo  <- c(0,breakg) 
	hi  <- c(breakg,100)
	mu  <- apply(rbind(lo,hi),2,mean)
	
	hh[hIndex] <- tnorm(nrow(hIndex),lo[sb],hi[sb],mu,.4)
			
	pnow <- pnew <- matrix(NA,n,nt)	
	
	#hh[hIndex] <- exp(rnorm(nrow(hIndex),log(hh[xIndex])*A + xmat%*%bg*dt,sqrt(dt*sigma)))
	
   pnow[hIndex] <- multiLogitStates(bq,ss[hIndex],h[hIndex],nclass)
   pnew[hIndex] <- multiLogitStates(bq,ss[hIndex],hh[hIndex],nclass) 
   
   pnow[hIndex] <- pnow[hIndex] + dnorm(log(h[hIndex]),log(h[xIndex])*A + 
                                        xmat%*%bg*dt,sqrt(dt*sigma))
   pnew[hIndex] <- pnew[hIndex] + dnorm(log(hh[hIndex]),log(hh[xIndex])*A + 
                                        xmat%*%bg*dt,sqrt(dt*sigma))
	pnow <- apply(pnow,1,sum,na.rm=T)
	pnew <- apply(pnew,1,sum,na.rm=T)
	a <- exp(pnew - pnow)
	z <- runif(n,0,1)
	hn[z < a,] <- hh[z < a,]
	list(hh = hn, ss = ss)
}
			
############################################


statePars <- function(b,discStates,contStates){

  nbreak <- length(breakg)
  accept <- 0

  fixEnds <- rbinom(1,1,.5)
  
  hislope <- max(b[,2])
  if(hislope > maxSlope)b[,2] <- b[,2] + maxSlope - hislope

  if(acount > 100)pv1 <- log(.001)
  
  propvars <- exp(rnorm(nbreak,pv1,1))
  
  db    <- diff(breakg)/2
  mids  <- breakg[-nbreak] + db 
  lo    <- c(breakLims[1,1],mids)
  lo[lo < breakLims[1,]] <- breakLims[1,lo < breakLims[1,]]
  hi    <- c(mids,98)
  hi[hi > breakLims[2,]] <- breakLims[2,hi > breakLims[2,]]

  breaks <- tnorm(nbreak,lo,hi,breakg,propvars)
  
  db <- diff(b[,1])/2
  wb <- which(db < .5)
  if(length(wb) > 0){
  	for(j in wb[1]:nbreak)b[j,1] <- max(b[j-1,1] + .5,b[j,2])
  	b  <- breaks2pars(b,breakg)
  	db <- diff(b[,2])/2
  }

  db <- diff(b[,2])/2
  wb <- which(db < 0)
  if(length(wb) > 0){
  	for(j in wb[1]:nbreak)b[j,2] <- max(b[j-1,2] + .01,b[j,2])
  	b  <- breaks2pars(b,breakg)
  	db <- diff(b[,2])/2
  }

  mids <- b[-nbreak,2] + db
  lo   <- c(b[1,2]-.2,mids)
  hi   <- c(mids,0)
  
  if(acount > 100)pv2 <- log(.00001)
  propvars <- exp(rnorm(nbreak,pv2,1))
  b2 <- tnorm(nbreak,lo,hi,b[,2],propvars)
  
  bpp <- breaks2pars(cbind(b[,1],b2),breaks)
  if(fixEnds == 1)bpp[c(1,nbreak),] <- bq[c(1,nbreak),]

  nss <- (nclass-1)*2
  if(acount > 20){
     sss <- sample(1:nss,2)
     bpp[sss] <- b[sss]
  }
  if(acount > 100){
     sss <- sample(1:nss,nclass)
     bpp[sss] <- b[sss]
  }
  if(acount > 200){
     sss <- sample(1:nss,nss-1)
     bpp[sss] <- b[sss]
  }

  discStates[is.na(discStates)] <- 0

  ww <- which(discStates > 0,arr.ind=T)

  pnow <- multiLogitStates(b,discStates,contStates,nbreak+1)
  pnew <- multiLogitStates(bpp,discStates,contStates,nbreak+1)

  pnow <- sum(pnow[ww],na.rm=T) + mydmvnorm(as.vector(b),as.vector(priorBQ),priorVBQ,log=T)
  pnew <- sum(pnew[ww],na.rm=T) + mydmvnorm(as.vector(bpp),as.vector(priorBQ),priorVBQ,log=T)

  z <- runif(1,0,1)
  r <- exp(pnew - pnow)

  if(z < r){
    b <- bpp
    breakg <- breaks
    accept <- 1
  }
  list(bgg = b, breakg = breakg, accept = accept)
}

##################################################3
predLogit <- function(hh,bqgibbs,priorBQ,transpose=F,yvals=numeric(0)){          #plot logit from MCMC chains

  nh   <- length(hh)
  nd   <- ncol(bqgibbs)/2
  nsim <- nrow(bqgibbs)
  lo <- hi <- mid <- hprior <- matrix(0,nh,(nd+1))

  for(i in 1:nh){

    st <- rep(0,nsim)
    sp <- 0

    for(k in 1:nd){

      eh <- exp(bqgibbs[,k] + bqgibbs[,(k+nd)]*hh[i])
      eh[eh > 1e+5] <- 1e+5
      tk <- eh/(1 + eh) - st
      st <- st + tk
      lo[i,k]  <- quantile(tk,.025)
      hi[i,k]  <- quantile(tk,.975)
      mid[i,k] <- quantile(tk,.5)

      ep <- exp(priorBQ[k,1] + priorBQ[k,2]*hh[i])
      ep[ep > 1e+5] <- 1e+5
      tp <- ep/(1 + ep) - sp
      sp <- sp + tp
      hprior[i,k] <- tp

    }
    lo[i,nd+1]  <- 1 - quantile(st,.025)
    hi[i,nd+1]  <- 1 - quantile(st,.975)
    mid[i,nd+1] <- 1 - quantile(st,.5)

    hprior[i,nd+1] <- 1 - sp
  }

  if(!transpose){

  plot(-10,0,type='l',lwd=2,xlim=c(0,100),ylim=c(0,1.09),ylab='Probability',xlab='Latent state h')


  for(k in 1:(nd+1)){

    wk <- which(mid[,k] > .001)
    hk <- hh[wk]
    lines(hk,mid[wk,k],lwd=3,col=colStates[k])
    lines(hk,lo[wk,k],lty=2,col=colStates[k],lwd=2)
    lines(hk,hi[wk,k],lty=2,col=colStates[k],lwd=2)

    lines(hk,hprior[wk,k],lty=3,col=colStates[k],lwd=2)
    wm <- which.max(mid[,k])

    text(hh[wm],1.06,k,col=colStates[k],cex=1.1)
  }
  }

  if(transpose){

  plotSetup(c(0,1),c(0,100),xvals = c(0,1),yvals=c(' ',' '),xlabel='Probability')
 # axis(side=4,at=c(0,50,100),labels=F)
  if(length(yvals) > 0)abline(h=yvals,lwd=3,col='white')

  for(k in 1:(nd+1)){

    wk <- which(mid[,k] > .001)
    hk <- hh[wk]
    lines(mid[wk,k],hk,lwd=3,col=colStates[k])
    lines(lo[wk,k],hk,lty=2,col=colStates[k],lwd=2)
    lines(hi[wk,k],hk,lty=2,col=colStates[k],lwd=2)

    lines(hprior[wk,k],hk,lty=3,col=colStates[k],lwd=2)
    wm <- which.max(mid[,k])

 #   text(.5,hh[wm],k,col=colStates[k],cex=1.1)
    text(.5,yvals[k],k,col=colStates[k],cex=1.1)
  }
  }


  list(lo = lo, hi = hi, mid = mid)
}

#################################
predProb <- function(b,h){

    pn  <- matrix(0,length(h),3)
    pn[,1] <- invlogit(b[c(1,3)],h)
    pn[,2] <- invlogit(b[c(2,4)],h) - invlogit(b[c(1,3)],h)
    pn[,3] <- 1 - invlogit(b[c(2,4)],h)
    pn
}


#################################################3
plotHP <- function(iplot,iyr=NULL){


  xlimit <- c(0,220)
  ylimit <- c(-6,6)

  xvalues <- c(1,moVec[moVec < max(xlimit)])

  yval1 <- c('dormant','fully expanded')
  ylab1 <- 'Continuous scale h(t)'

  letters <- matrix(c('a) ','b) ','c) ','d) ','e) ','f) ','g) ','h) ','i) ','j) '),2)

  par(mfcol=c(2,(length(iplot)+1)),bty='n',mar=c(5,5,2,1),cex.axis=1.2,cex.lab=1.2)

 if(is.null(iyr)){
   iyr <- numeric(0)
   for(k in 1:length(iplot))iyr <- c(iyr,years[wy + k])
 }

 for(i in 1:length(iplot)){

  if(i > length(iplot))next

  wy <- which(years == (iyr[i] - 1))
  wg <- findInterval(germIndex[iplot[i],2],yrBreak) + 1
  wg <- wg + i - 1
  if(wg > wy){
     wy <- wg
     if(wy > length(iplot))next
  }
  wd <- findInterval(diedIndex[iplot[i],2],yrBreak) + 1

  plotyr <- years[wy + 1]

  print1 <- yrBreak[wy] 
  print2 <- yrBreak[wy + 1] 
  if(is.na(print2))print2 <- nt
  wii    <- c(print1:print2)

  wij <- wii
  jdh <- jdAll[wij] - 365*(jdIndex[wij] -1)
  hmu <- hg[iplot[i],wij]/ntot
  hse <- sqrt(hg2[iplot[i],wij]/ntot - hmu^2)

  if(!is.finite(max(hmu,na.rm=T)) | (max(hmu,na.rm=T) < 50))next

  smu <- stateSum[iplot[i],wij]/ntot

  wj  <- which(c(0,diff(jdh)) < 0 | c(0,diff(hmu)) < 0 | is.na(hmu))
  if(length(wj) > 0){
     wij <- wii[-wj]
     jdh <- jdAll[wij] - 365*(jdIndex[wij] -1)
     hmu <- hmean[iplot[i],wij]
     hse <- sqrt(hg2[iplot[i],wij]/ntot - hmu^2)
     smu <- stateSum[iplot[i],wij]/ntot
  }

  if(i > 1){
    yval1 <- rep(' ',2)
    ylab1 <- ' '
  }

  par(plt=c(.35,1,.2,.8))
  plotSetup(xvalues,c(0,100),xvals = rep(' ',length(xvalues)),
            yvals=,yval1,xlabel=' ',ylabel=ylab1)
  axis(1,xvalues[-1]-15,names(xvalues)[-1],tick=F)

  hse[1] <- 1
  hse[hse < 0] <- 0

  hse[is.na(hse)] <- mean(hse,na.rm=T)
  loh <- hmu - 1.96*hse
  loh[loh < 0] <- 0
  hhi     <- hmu + 1.96*hse
  hhi[hhi > 100] <- 100

  dat <- dataMat[iplot[i],wij]

  bqpar <- matrix(apply(bqgibbs,2,mean),nclass-1,2)

  shs <- (c(1:nclass) - 1)*100/nclass
  sds <- shs*100/max(shs)
  

  for(k in 1:nclass){
     if(max(hmu) < 90 & k == nclass)break
     pk <- exp(multiLogitStates(bqpar, rep(k,length(wij)), hmu,nclass))
     dwt <- sum(hmu*pk)/sum(pk)
     if(k == 1)dwt <- 0
     kj <- which(pk > 1e-10)
     y1 <- 10*pk[kj] + dwt
     y2 <- 0*pk[kj] + dwt
     abline(h=dwt,lwd=3,col='white')
     lines(jdh[kj],y1,col=colStates[k],lwd=2)
     lines(jdh[kj],y2,col=colStates[k],lwd=2)
     polygon(c(jdh[kj],rev(jdh[kj])),c(y1,rev(y2)),col=colStates[k],border=colStates[k])
     xpos <- max(jdh[kj])+10
     if(k > 4)xpos <- min(jdh[kj])-10
     if(i == 1)text(xpos,dwt,k,col=colStates[k],cex=1.3)
     sds[k] <- dwt
  }
  id <- paste(site[iplot[i]],treat[iplot[i]],iplot[i],sep=', ')

  text(.2*min(jdh),105,paste(letters[1,i],id,sep=''),cex=1.25,pos=4)
  lines(jdh,hmu,lwd=5,col='white')
  lines(jdh,hhi,lwd=4,col='white')
  lines(jdh,loh,lwd=4,col='white')
  lines(jdh,hmu,lwd=2)
  lines(jdh,hhi,lty=2)
  lines(jdh,loh,lty=2)

  points(jdh,sds[dat],col='white',cex=1.2,lwd=4)
  points(jdh,sds[dat],col=colStates[dat],cex=.8,lwd=3)
  points(jdh,sds[dat],col=colStates[dat],cex=.5,lwd=3)

  hm  <- diff(hmu)/dt
  sr  <- round(smu,0)
  ww  <- which(sr < nclass)
  ns  <- length(ww)

  hm <- hm[ww]
  sr <- sr[ww]

  wwg <- sample(100:nrow(bggibbs),2000,replace=T)
  gmat2 <- gmat3 <- matrix(NA,length(wwg),ns)

  for(gg in 1:length(wwg)){

     hgg <- tnorm(ns,0,100,hmu[ww],hse[ww])
     dh <- diff(hgg)/dt

     hm  <- diff(hmu)/dt

     bg <- matrix(bggibbs[wwg[gg],],ncol=1)
     bq <- matrix(bqgibbs[wwg[gg],],ncol=2)

        chk <- bq[1,1] + bq[1,2]*hgg
        dP2 <- -bq[sr,2]*exp(chk)/(1 + exp(chk))^2           #current state

        chk <- bq[2,1] + bq[2,2]*hgg
        dP3 <- -bq[sr,2]*exp(chk)/(1 + exp(chk))^2           #current state
        ba  <- matrix(bg[xnames == 'AT'],ns,1)

        if('ATXseed' %in% xnames & seed[iplot[i]] == 1){
             ba <- ba + bg[xnames == 'ATXseed']*x3d[iplot[i],wii[ww],xnames == 'AT']
        }

        if('ATXgap' %in% xnames & gap[iplot[i]] == 1){
            ba <- ba + bg[xnames == 'ATXgap']*x3d[iplot[i],wii[ww],xnames == 'AT']
        }

        gmat2[gg,] <- dP2*dh*ba
        gmat3[gg,] <- dP3*dh*ba

   }

  ci2 <- apply(log10(gmat2),2,quantile,c(.5,.025,.975),na.rm=T)
  ci3 <- apply(log10(gmat3),2,quantile,c(.5,.025,.975),na.rm=T)
  ym <- round(max(c(ci2,ci3),na.rm=T),0)
  yvalues2 <- seq((ym - 8),ym,by=2)

  par(plt=c(.35,1,.55,1))

  ylab2 <- expression(hat(gamma)[ity])
  yval2 <- 10^yvalues2

  if(i > 1){
    yval2 <- rep(' ',length(yvalues2))
    ylab2 <- ' '
  }

  plotSetup(xvalues,yvalues2,xvals = rep(' ',length(xvalues)),yvals=yval2,
            xlabel=plotyr,ylabel=ylab2)
  axis(1,xvalues[-1]-15,names(xvalues)[-1],tick=F)
 

  tmax <- 30
  tmin <- -10
  tscale <- (tmax - tmin)/(yvalues2[length(yvalues2)] - yvalues2[1])
  tvalues <- tmin + (yvalues2 - yvalues2[1])*tscale

  temp <- x3d[iplot[i],wij,'AT'] + meanX['AT']

  if(i == length(iplot)){
     axis(4,yvalues2,labels=tvalues)
     mtext(side=4,degreeLab,line=3,cex=1.1)
  }

  abline(h=yvalues2[tvalues == 0],lwd=2,col='grey')
  lines(jdh,yvalues2[1] + (temp - tmin)/tscale,lwd=2)


  text(.2*min(jdh),ym,paste(letters[2,i],'Sensitivity',sep=''),cex=1.25,pos=4)

  lines(jdh[ww],ci2[1,],lwd=3,col=colStates[2])
  lines(jdh[ww],ci2[2,],lty=2,col=colStates[2])
  lines(jdh[ww],ci2[3,],lty=2,col=colStates[2])
  lines(jdh[ww],ci3[1,],lwd=3,col=colStates[3])
  lines(jdh[ww],ci3[2,],lty=2,col=colStates[3])
  lines(jdh[ww],ci3[3,],lty=2,col=colStates[3])

}

  par(plt=c(0,.3,.2,.8))
  tmp <- predLogit(hseq,bqgibbs,priorBQ,transpose=T,yvals=sds)
  text(.3,105,'c) Species',cex=1.25)

}

predVsObsStates <- function(){

  par(mfrow=c(1,1),bty='n',cex.axis=1.2,cex.lab=1.2)

  xlimit <- ylimit <- c(.5,nclass+.5)

  wp <- which(is.finite(predvals) & is.finite(dataMat) & predvals > 0 & dataMat > 0,arr.ind=T)
  plot(jitter(dataMat[wp]),predvals[wp],col='grey',xlim=xlimit,ylim=ylimit,
         xlab='Observed',ylab='Predicted',cex=.5,,xaxt='n',yaxt='n')
  axis(side=1,at=c(1:nclass),labels=c(1:nclass))
  axis(side=2,at=c(1:nclass),labels=c(1:nclass))

  rect(xlimit[1],ylimit[1],xlimit[2],ylimit[2],col='azure2',border=NA)
  abline(v=c(1:nclass),lwd=2,col='white')
  abline(h=c(1:nclass),lwd=2,col='white')
  abline(0,1,lwd=3,col='grey',lty=2)

  points(jitter(dataMat[wp]),predvals[wp],col='grey50',cex=.5)

    bw <- .1/2

    for(j in 1:nclass){
      wj <- which(dataMat == j,arr.ind=T)
      pj <- predvals[wj]
      pj <- pj[pj != 0]
      qj <- quantile(pj,c(.5,.025,.975),na.rm=T)
      lines(c(j,j),qj[2:3],lwd=2)
      points(j,qj[1],pch=3,cex=2)
      polygon(c(j-bw,j+bw,j+bw,j-bw),c(qj[1]-bw,qj[1]-bw,qj[1]+bw,qj[1]+bw),col=1)
      qj <- quantile(pj,c(.159,.841),na.rm=T)
      lines(c(j,j),qj,lwd=6)
    }

    text(nclass-1,1,gibbsSpec)
}

plotH <- function(){
	
  hmean <- hg/ntot
  hse   <- hg2/ntot - hmean^2
  hlo   <- hmean - 1.96*hse
  hhi   <- hmean + 1.96*hse
  hlo[hlo < 0] <- 0
  hhi[hhi > 100] <- 100
  
  nh <- 6
  hs <- sample(c(1:n),nh)

  par(mfrow=c(3,2),mar=c(5,4,3,2),bty='n')
  for(j in hs){
    xj <- c(germTime[j],diedTime[j])
    xj[xj > max(jdAll)] <- max(jdAll)
    if(is.na(xj[2]))xj[2] <- max(jdAll)

    plot(jdAll,(predvals[j,] - 1)*100/(nclass-1),type='l',ylim=c(0,100),lwd=2,xaxt='n',col='brown')
    axis(1,at=yrVec)
    lines(jdAll,hmean[j,],lwd=2)

    lines(jdAll,hlo[j,],lty=2)
    lines(jdAll,hhi[j,],lty=2)
    points(jdAll,(dataMat[j,] - 1)*100/(nclass-1),col=2)
    abline(v=c(germTime[j],diedTime[j]),col=3,lwd=2)
    lines(xj,c(0,0),lwd=2,col=3)
    title(j)
  }
}

predX <- function(tplus,seedx=0,gapx=0,smoothVals=character(0),meanVals=character(0),yrPred=3){  
	
  #construct x matrix for prediction with delta tplus
  #uses year 2011

  xpred <- matrix(NA,nt,nx)
  tpred <- rep(0,nt)

  end <- yrBreak[yrPred]
  if(is.na(end)){
    end <- max(jdAll)
  }

  start <- yrBreak[yrPred - 1]

  if(end > dim(x3d)[2])end <- dim(x3d)[2]

  for(k in 1:nx){
    wk <- which(dimnames(x3d)[[3]] == xnames[k])
    xk <- x3d[,start:end,wk]        #*(tmat*0 + 1)
    xk[xk == 0] <- NA
    xxe <- apply(xk,2,quantile,.5,na.rm=T)
    xxe[is.na(xxe)] <- .0001
    xpred[start:end,k] <- xxe
    if(xnames[k] %in% smoothVals){
    	wx <- which(is.finite(xpred[,k]) & xpred[,k] != 0)
    	xx <- xpred[wx,k]
        xx <- lowess(xx,f=.01)
    	xpred[wx,k] <- xx$y
    }
  }
  colnames(xpred) <- xnames
  xpred[,xnames == 'AT'] <- xpred[,xnames == 'AT'] + tplus
  if(length(meanVals) > 0){
    xpred[,xnames %in% meanVals] <- apply(xpred,2,mean,na.rm=T)[xnames %in% meanVals]
  }
  
  xpred[,xnames == 'int']     <- 1 
  xpred[,xnames == 'ATXSM']   <- xpred[,xnames == 'AT']*xpred[,xnames == 'SM']
  if('seed' %in% xnames){
    xpred[,xnames == 'seed']    <- seedx
    xpred[,xnames == 'ATXseed'] <- xpred[,xnames == 'AT']*xpred[,xnames == 'seed']
  }
  if('gap' %in% xnames){
    xpred[,xnames == 'gap']    <- gapx
    xpred[,xnames == 'ATXgap'] <- xpred[,xnames == 'AT']*xpred[,xnames == 'gap']
  }

  xpred <- xpred[start:end,]
  tpred <- times[start:end]
  tpred <- jdAll[tpred]  - (yrPred-1)*365
  list(tpred = tpred, xpred = xpred)
}

sensk <- function(kstage,temp,seedval){

   c03 <- bqgibbs[,kstage]
   c13 <- bqgibbs[,kstage+nclass-1]

   xm <- matrix(0,nrow(bggibbs),length(xnames))
   colnames(xm) <- xnames
   xm[,1] <- 1
   xm[,'AT'] <- temp
   xm[,'h']  <- -c03/c13
   if('seed' %in% xnames)xm[,'seed'] <- seedval
   if('ATXseed' %in% xnames)xm[,'ATXseed'] <- xm[,'AT']*xm[,'seed']
   bigb <- apply(exp(xm*bggibbs),1,sum)

   bgg <- bggibbs[,'AT']
   if('ATXseed' %in% xnames)bgg <- bgg + bggibbs[,'ATXseed']*seedval  #only interaction is 'seed'
   sr <- quantile(-c13*bigb*bgg/4,c(.5,.025,.975))
   names(sr) <- paste('rate',seedval,names(sr),sep='-')

    meanTime <- -c03/bigb/c13
   st        <- c03/bigb/c13*bgg
   st <- quantile( st/meanTime ,c(.5,.025,.975))
   names(st) <- paste('time',seedval,names(st),sep='-')

   list(sRate = sr, sTime = st)
}


dpdh <- function(bqChain,kstage,hh=NULL){

   HALF <- F
   kj   <- c(1:kstage)
   cint <- bqChain[,kj]
   cslp <- bqChain[,nclass+kj-1]

   if(is.null(hh)){
     HALF <- T
     hh <- -cint[,kstage]/cslp[,kstage]
   }

   sumt <- sumk <- hh*0
   for(j in 1:kstage){

      eh <- exp(cint[,kj[j]] + cslp[,kj[j]]*hh)
      tk <- eh/(1 + eh) - sumt
      pr <- tk
      sumt <- sumt + tk

      ek <- cslp[,kj[j]]*eh/(1 + eh)^2 - sumk
      dp <- ek
      sumk <- sumk + ek 
   }

   dPdh <- -ek
   if(HALF){
     dPdh <- -cslp[,kstage]/4
   }

   list(dp = dp,pr = pr, Pr = eh/(1 + eh), dPdh = dPdh,h.5 = hh)
}



dpdt <- function(h0,dt,x,bqgibbs,bggibbs,kstage){  #x - sequence of x's

        nss <- 5000
        gs  <- sample(nrow(bqgibbs),nss,replace=T)

        kj   <- c(1:kstage)
        cint <- bqgibbs[gs,kj]
        cslp <- bqgibbs[gs,nclass+kj-1]

	ntt <- nrow(x)
	ht <- rep(h0,nss)
	pdout <- hout <- pout <- matrix(NA,ntt,3)

        ss <- sgibbs[gs]
        p50 <- 0

        DD <- rep(0,ntt)
	
	for(t in 2:ntt){

         if('h' %in% xnames) x[t,'h'] <- mean(ht)
         if('h2' %in% xnames)x[t,'h2'] <- x[t,'h']^2
         DD[t] <- DD[t-1] + max( x[t,'AT'] + meanX['AT'],0)
		
         ht   <- ht + exp(x[t,]%*%t(bggibbs[gs,]) + rnorm(length(ht),0,sqrt(ss)))*(1 - ht/100)*dt
         ht[ht > 100] <- 100
         ht[ht < h0]  <- h0

         tmp <- dpdh(bqgibbs[gs,],kstage,ht)
         pr  <- tmp$pr
         dp  <- tmp$dp

        dhdt  <- exp(x[t,]%*%t(bggibbs[gs,]) + rnorm(length(ht),0,sqrt(ss)))*(1 - ht/100)
        dpdt <- dp*dhdt 

      pdout[t,] <- quantile(pr,c(.5,.025,.975))
      pout[t,] <- quantile(dpdt,c(.5,.025,.975))
      hout[t,] <- quantile(ht,c(.5,.025,.975))
    }
    list(prob = pdout, dpdt = pout, h = hout, p50 = p50,DD = DD)
}

TeffectOnP3 <- function(bqgibbs,bggibbs){
	
	c0 <- bqgibbs[,3]
	c1 <- bqgibbs[,9]
	bT <- bggibbs[,xnames == 'AT']
	xi <- grep('ATX',xnames)
	bX <- bggibbs[,xi]
	
	hp5  <- -c0/c1
	hptq <- quantile(hp5,c(.5,.025,.975))
	dpdh <- c1
	dpdhq <- quantile(c1,c(.5,.025,.975))
	dhdt <- 1/dt * (-c0/c1)^A*(exp(meanX%*%t(bggibbs)) - 1)
	dhdtq <- quantile(dhdt,c(.5,.025,.975))
	dtdhdt <- 1/dt * (-c0/c1)^A*(bT + meanX[xi]%*%t(bggibbs[,xi]))*(exp(meanX%*%t(bggibbs)) - 1)
	dtdhdtq <- quantile(dtdhdt,c(.5,.025,.975))

	ss <- c1/4/dt * (-c0/c1)^A*(bT + meanX[xi]%*%t(bggibbs[,xi]))*(exp(meanX%*%t(bggibbs)) - 1)
	
	quantile(ss,c(.5,.025,.975))
}


dPdh <- function(c01,c11,h){  #d(1 - P)/dh
  zk <- exp(c01 + c11*h)
  -c11*zk/(1 + zk)^2
}


dhdt <- function(x,beta,hh){

  exp(x%*%beta)*(1 - hh/100)
}

prob2h <- function(c0k,c1k,p){   #h when at least stage k is reached
  (log(p/(1 - p)) - c0k)/c1k
}

b2prob <- function(c0k,c1k,c00,c10,bggibbs,xx=meanX,pk,VNAME='AT'){  
     #effect of VNAME on rate of prob P when P = pk

  if(!is.matrix(xx))xx <- matrix(xx,nrow=1)

  hh   <- prob2h(c0k,c1k,pk)
  d1   <- dPdh(c0k,c1k,hh)
  d2   <- t(dhdt(hh,xx,t(bggibbs)) )
  dpdt <- d1*d2
  wq   <- grep(VNAME,xnames)
  wi   <- intersect(wq,grep('X',xnames))
  if(length(wi) > 0)wq <- wq[!wq %in% wi]
  dd   <- bggibbs[,wq]
  if(length(wi) > 0)dd <- dd + bggibbs[,wi]*xx[wi]
  sens <- d1*dd*exp(xx%*%t(bggibbs))
  quantile(sens,c(.5,.025,.975))
}

simX <- function(xlim,ylim,seed){

  ndd   <- diff(xlim) + 1
  xpred <- matrix(1,ndd,nx)
  colnames(xpred)   <- xnames
  if('seed' %in% xnames)xpred[,'seed']    <- seed
  xpred[,'AT']      <- ylim[1] + diff(ylim)/diff(xlim)*(c(xlim[1]:xlim[2]) - xlim[1])
  if('ATXseed' %in% xnames)xpred[,'ATXseed'] <- seed*xpred[,'AT']
  list(xpred = xpred,dseq = c(xlim[1]:xlim[2]))
}

plotData <- function(){

   sites <- c('DF','HF')

   par(mfcol=c(nclass,2),bty='n',mar=c(4,4,1,1)+.1)

   for(j in 1:2){

   for(k in 1:nclass){

     ylabel <- ' '
     xlabel <- ' '
     if(k == 4 & j == 1)ylabel <- 'Density by yr'
     if(k == nclass)xlabel <- 'Year'

     mmat <- dataMat*0
     mmat[dataMat == k] <- 1
     wa <- which(specData[,'Am'] == 1)
     mmat[wa,] <- 0
     wf <- which(specData[,'site'] == sites[j])
     mmat[wf,] <- 0
     mmat <- apply(mmat,2,sum,na.rm=T)

     for(y in 1:(1+length(yrBreak))){
       sumy <- sum(mmat[jdIndex == y])
       mmat[jdIndex == y] <- mmat[jdIndex == y]/sumy/dt
     }

     plot(jdAll,mmat + .1,type='s',col='orange',ylim=c(0,.5),ylab=ylabel,xlab=xlabel,xaxt='n')
     if(k == nclass)axis(1,at=yrVec[1:length(years)]-365/2,labels=years)

     mmat <- dataMat*0
     mmat[dataMat == k] <- 1
     wa <- which(specData[,'Am'] == 0)
     mmat[wa,] <- 0
     mmat[wf,] <- 0 

     mmat <- apply(mmat,2,sum,na.rm=T)
     for(y in 1:(1 + length(yrBreak))){
       sumy <- sum(mmat[jdIndex == y])
       mmat[jdIndex == y] <- mmat[jdIndex == y]/sumy/dt
     }

     lines(jdAll,mmat,type='s',col='blue')
     if(j == 1 & k == 1)title(specList[jjj])
     abline(v=yrVec,lty=3)
     if(j == 1 & k == 1)text(30,.10,'Harvard Forest',pos=4)
     if(j == 2 & k == 1)text(30,.10,'Duke Forest',pos=4)
  }
  }
  
  dev.print(file=paste(specList[jjj],'springPhen.ps',sep='-'),
	          width=6,height=9,horizontal=F)  
}
#######################################
plotLowTemp <- function(x,threshold=0,lo=T){

  winter <- c('nov','dec','jan','feb','mar')
  moday  <- c(334,365,396,424,455)
  xlab   <- xlab2 <- ' '
  xvals  <- rep(' ',length(winter))

  ccols <- union(grep('S',colnames(x)),grep('G',colnames(x)))
  mu    <- apply(x[,ccols],1,mean,na.rm=T)
  rt    <- t(apply(x[,ccols],1,range,na.rm=T))

  temp <- x[,ccols]

  tloc <- c(.96,.80)
  tseq <- seq(0,300,length=40)
  nt   <- length(threshold)

  par(mfrow=c(4,2),mar=c(5,5,1,1),bty='n')

  for(j in 1:length(years)){

    chillyk <- numeric(0)

    par(plt=c(.2,.9,.25,.9))
    wj <- numeric(0)
    if(j > 1)wj <-which(x[,'yeary'] == years[j-1] & x[,'monthy'] %in% c('nov','dec'))
    wj <- c(wj, which(x[,'yeary'] == years[j]   & x[,'monthy'] %in% c('jan','feb','mar')) )

    if(j == length(years)){
      xlab  <- 'Year'
      xvals <- winter
      xlab2 <- 'Chilling units'
    }
    plotSetup(moday-365,c(0,1),xvals=xvals,
                 xlabel=xlab,ylabel='Fraction of treatments')
   title(years[j])

   for(k in 1:nt){
 
    x0 <- temp[wj,]

    chill <- temp[wj,] - threshold[k]
    chill[chill > 0] <- 0

    ck <- -apply(chill,2,sum,na.rm=T)
    chillyk <- rowBind(chillyk,ck,paste(j,k,sep='-'))

    x0[temp[wj,] <= threshold[k]] <- 1
    x0[temp[wj,] > threshold[k]]  <- 0
  
    xt  <- x0*0 + 1
    xlo <- apply(x0,1,sum,na.rm=T)/apply(xt,1,sum,na.rm=T)
    xd  <- as.numeric(names(xlo))

    lines(xd-365*(j-1),xlo,type='s',col=k)

    if(j == 1){
       text(moday[1]-365,tloc[k],substitute( list(tau) == list(xx),list(xx=threshold[k])),
             col=k,cex=1.5)
    }
  }

   par(plt=c(.15,.6,.25,.9))
   plotSetup(seq(0,max(tseq),by=100),c(0,.15),
                 xlabel=xlab2,ylabel='Density')
   for(k in 1:nt){
  
     tmp <- hist(chillyk[k,],breaks=tseq,plot=F)
     xx  <- tmp$mids
     yy <-  tmp$density
     lines(xx,yy,type='s',col=k)
  }

 }
}

linearModel <- function(maxSens,xn){

  cvar <- c('ddwinter','chill')    #construct ellipse for these variables

  mu <- cr <- el <- numeric(0)

  ns <- length(specList)
  for(j in 1:ns){

    wj <- intersect( grep(specList[j],rownames(maxSens)), which(is.finite(maxSens[,'ddwinter'])))
    if(length(wj) < 10)next
    y  <- maxSens[wj,'bbDate']
   # y  <- y - mean(y)
    x  <- maxSens[wj,xn]
print(cor(x)[1,2])

 #   x[,'ddwinter'] <- x[,'ddwinter'] - mean(x[,'ddwinter'],na.rm=T)
    tmp <- lm(y ~ x)
    
    int   <- y*0 + 1
    xx    <- cbind(int,x)
    covar <- solve(crossprod(xx))
    b     <- covar%*%crossprod(xx,y)
    res   <- sum((y - xx%*%b)^2)/(length(y) - ncol(xx))
    covar <- covar*res

   tmp   <- list(loc = b[cvar,], cov = covar[cvar,cvar], d2 = qchisq(.9,1) )
   tmp1  <- predict.ellipsoid(tmp,n.out=201) 
   mu    <- rowBind(mu,as.vector(b),specList[j])
   nj    <- paste(specList[j],c('x','y'),sep='-')
   colnames(tmp1) <- nj
   el    <- rbind(el,t(tmp1))
   cr    <- c(cr,covar[cvar[1],cvar[2]]/sqrt( prod(diag(covar[cvar,cvar])) ))
   names(cr)[length(cr)] <- specList[j]
  }

  colnames(mu) <- c('int',cvar)

  list(mu = mu, cor = cr, ellipse = el)
}

#####################################3
parTable <- function(){

  plotStage <- 3

  xx    <- meanX
  xx['AT'] <- 10
  if('seed' %in% xnames)   xx['seed'] <- 0
  if('ATXseed' %in% xnames)xx['ATXseed'] <- 0
  if('ATXgap' %in% xnames) xx['ATXgap'] <- 0

  hh  <- 30
  gs  <- sample(nrow(bqgibbs),2000,replace=T)
  tmp <- dpdh(bqgibbs[gs,],plotStage)
  dP  <- tmp$dPdh    #change in Pr with h at Pr = 1/2
  h.5 <- tmp$h.5     #h at Pr = 1/2
  dh  <- dhdt(matrix(xx,1),t(bggibbs[gs,]),hh=h.5)
  ba  <- bggibbs[gs,'AT']

  seedgap <- quantile(dP*dh*ba,c(.5,.025,.975))

  if('ATXseed' %in% xnames){
    ba <- bggibbs[gs,'AT'] + bggibbs[gs,'ATXseed']*xx['ATXseed']
    xs <- xx
    xs['seed'] <- 1
    xs['ATXseed'] <- xx['AT']
    dh <- dhdt(matrix(xs,1),t(bggibbs[gs,]),hh=h.5)
    seedNorthgap0 <- quantile(dP*dh*ba,c(.5,.025,.975))
    seedgap <- rbind(seedgap,seedNorthgap0)
  }
  if('ATXgap' %in% xnames){
    ba <- bggibbs[gs,'AT'] + bggibbs[gs,'ATXgap']*xx['ATXgap']
    xs <- xx
    xs['gap'] <- 1
    xs['ATXgap'] <- xx['AT']
    dh <- dhdt(matrix(xs,1),t(bggibbs[gs,]),hh=h.5)
    seedSouthgap1 <- quantile(dP*dh*ba,c(.5,.025,.975))
    seedgap <- rbind(seedgap,seedSouthgap1)
  }
  if('ATXgap' %in% xnames & 'ATXseed' %in% xnames){
    ba <- bggibbs[gs,'AT'] + bggibbs[gs,'ATXseed']*xx['ATXseed'] + bggibbs[gs,'ATXgap']*xx['ATXgap']
    xs <- xx
    xs['seed'] <- 1
    xs['ATXseed'] <- xx['AT']
    xs['gap'] <- 1
    xs['ATXgap'] <- xx['AT']
    dh <- dhdt(matrix(xs,1),t(bggibbs[gs,]),hh=h.5)
    seedNorthgap1 <- quantile(dP*dh*ba,c(.5,.025,.975))
    seedgap <- rbind(seedgap,seedNorthgap1)
  }

  bg    <- matrix(apply(bggibbs,2,median),ncol=1)
  bq    <- matrix(apply(bqgibbs,2,median),ncol=2)
  sigma <- median(sgibbs)

  meanDev <- sumdev/ntot
  devMean <- updateStates()$dev
  pd      <- meanDev - devMean
  dic     <- 2*pd + meanDev

  sig    <- sout[1,]
  parOut <- rbind(bgout[,1:3],sig,seedgap)

  n_dic      <- rep(NA,nrow(parOut))
  n_dic[1:2] <- c(n,dic)
  parOut <- cbind(parOut,n_dic)
  rownames(parOut) <- paste(specList[jspec],rownames(parOut),sep='_')
  rownames(bqout)  <- paste(specList[jspec],rownames(bqout),sep='_')

  tmp <- lmFit()
  lfit <- tmp$coeff
  rownames(lfit) <- paste(gibbsSpec,'lm',rownames(lfit),sep='_')
  ci <- cbind(lfit[,1] - 1.96*lfit[,2],lfit[,1] + 1.96*lfit[,2])
  lfit <- cbind(lfit[,1],ci,c(tmp$rsq,rep(NA,nrow(lfit)-1)))

  parOut <- rbind(parOut,lfit)

 # write.table(bqout,file=paste(filename,'pars.txt',sep='_'),quote=F,
 #           row.names=T,sep=',')
            
 # write.table(bgout,file=paste(filename,'parsG.txt',sep='_'),quote=F,
 #           row.names=T,sep=',')

 # write.table(sensOut,file=paste(filename,'pars.txt',sep='_'),quote=F,
 #           row.names=T,sep=',')

 # write.table(parOut,file=paste(filename,'allPars.txt',sep='_'),quote=F,
 #           row.names=T,sep=',')

  list(parOut = parOut, dic = dic)
}
######################################
devTime <- function(tplus,yrPred,stage,maxDay,seedx=0,gapx=0,smoothVals=character(0)){

  tmp  <- predX(0,seedx,gapx,smoothVals,yrPred)
  xamb  <- tmp$xpred
  tpred <- tmp$tpred

  tmp <- dpdt(1,dt,xamb,bqgibbs,bggibbs,stage)
  pra <- tmp$prob
  hha <- tmp$h
  ddAmb <- tmp$DD

  xele  <- xamb
  xele[,'AT'] <- xamb[,'AT'] + tplus
  if('ATXseed' %in% xnames & seedx == 1)xele[,'ATXseed'] <- xele[,'AT'] *xele[,'seed']
  if('ATXgap ' %in% xnames & gapx == 1)xele[,'ATXgap'] <- xele[,'AT'] *xele[,'gap']

  tmp <- dpdt(1,dt,xele,bqgibbs,bggibbs,stage)
  pre <- tmp$prob
  hhe <- tmp$h
  ddEl <- tmp$DD

  hout <- cbind(hha,hhe)
  pout <- cbind(pra,pre)
  colnames(hout) <- colnames(pout) <- as.vector(t(outer(c('amb','ele'),c('0.5','.025','.975'),paste,sep='-')))

  list(jd = tpred, x = cbind(xamb[,'AT'],xele[,'AT']), h = hout, pr = pout, dd = cbind(ddAmb,ddEl))
}

plotDev <- function(x,y1,y2,xl=NULL,yl=NULL,xlabel=' ',ylabel=' ',leg=NULL,ttext=NULL,colors=1,add=F,
           xlines=0,ylines=0){

  if(is.null(xl))xl <- range(x,na.rm=T)
  if(is.null(yl))yl <- range(c(y1,y2),na.rm=T)

  y11 <- y1
  y22 <- y2
  if(is.matrix(y1))y11 <- y1[,1]
  if(is.matrix(y2))y22 <- y2[,1]

  if(!add){
 
     plot(x, y11,type='l',col=colors[1],xlim=xl,lwd=3,ylim=yl,xlab=xlabel,ylab=ylabel)
     rect(xl[1],yl[1],xl[2],yl[2],col='azure2',border=NA)
     if(!is.null(xlines))abline(v=xlines,lwd=3,col='white')
     if(!is.null(ylines))abline(h=ylines,lwd=3,col='white')
     lines(x, y11,col=colors[1],lwd=3)
     lines(x,y22,col=colors[2],lwd=3)
  }
  if(add)lines(x,y11,type='l',col=colors[1],lwd=3)

  if(!is.null(ttext))title(ttext)

  if(is.matrix(y1))for(j in 2:3)lines(x,y1[,j],lty=2,col=colors[1],lwd=2)
  if(is.matrix(y2))for(j in 2:3)lines(x,y2[,j],lty=2,col=colors[2],lwd=2)

  if(!is.null(leg))legend('bottomright',leg,text.col=colors,bty='n',cex=1.4)

}

######################################################
plotDDSens <- function(bggibbs,bqgibbs,PLOT=F){

  #  colF <-  colorRampPalette(c("blue","orange","red"))

    DDdate <- bbDate <- bbDD <- maxSens2 <- maxSens3 <- maxDate2 <- maxDate3 <- gapDate <- seedDate <- 
    yearDate <- siteDate <- chillDate <- index <- trt <- ddwinter <- numeric(0)
    chamDate <- character(0)
    JDmat  <- matrix(jdAll - (jdIndex-1)*365,n,nt,byrow=T)
    firstbud <- numeric(0)

    endWinter <- 119    #jd for 31 March

    if('chill' %in% xnames){
        chillMat <- ddmat*0
        chillMat[xIndex] <- xmat[,'chill']
    }

    miny <- 1e-14

    bq <- matrix(apply(bqgibbs,2,mean),nclass-1,2)
    bg <- matrix(apply(bggibbs,2,mean),ncol=1)

    plot(c(1,1),c(1,1),xlim=c(1,1000),ylim=c(miny,1),log='xy',
          xlab='Degree days',ylab='Temperature effect on rate',bty='n')

    for(y in 2:nyr){

        wt <- which(jdIndex == y)
        ntt <- length(wt)
        if(ntt == 0)next

        wi <- which(germYr < y)
        ni <- length(wi)
        if(ni == 0)next

        hh <- row2Mat(hmean[wi,wt])
        ss <- row2Mat(smean[wi,wt])

        dd <- row2Mat(ddmat[wi,wt])
        jd <- row2Mat(JDmat[wi,wt])
        pd <- row2Mat(dataMat[wi,wt])       #observations

        di       <- findInterval(endWinter,jd[1,])  #DD for 31 March
        ddwin <- dd[,di]
        ddwin[ddwin == 0] <- NA
        ddwinter <- c(ddwinter,ddwin)
 
        chill <- chillMat[wi,wt]

        yy <- matrix(1:length(wt),ni,length(wt),byrow=T)
        yy[smean[wi,wt] < 2 | is.na(smean[wi,wt])] <- NA
        s2  <- apply(yy,1,min,na.rm=T)                         #reach state 2
        s2i <- cbind(c(1:ni),s2)
        s2i <- s2i[s2i[,2] > 1 & s2i[,2] < length(wt),]
        if(!is.matrix(s2i))s2i <- matrix(s2i,1,2)

        yy <- matrix(1:length(wt),ni,length(wt),byrow=T)
        yy[smean[wi,wt] < 3 | is.na(smean[wi,wt])] <- NA
        s3  <- apply(yy,1,min,na.rm=T)                         #reach state 3
        s3i <- cbind(c(1:ni),s3)
        s3i <- s3i[s3i[,2] > 1 & s3i[,2] < length(wt),]
        if(!is.matrix(s3i))s3i <- matrix(s3i,1,2)

        pd[pd > 3] <- 3
        pd[,1] <- 1
        pd[is.na(pd)] <- -999
        pf <- apply(pd,1,which.max)
        sf <- pd[cbind(c(1:ni),pf)]
        pf[sf != 3] <- NA      
        pf3 <- jd[1,pf]         #date first observed >= 3

        pd[pd > 2] <- 2
        pf <- apply(pd,1,which.max)
        sf <- pd[cbind(c(1:ni),pf)]
        pf[sf != 2] <- NA      
        pf2 <- jd[1,pf]         #date first observed >= 2

        dh  <- t(apply(hh,1,diff))
        hm  <- row2Mat(hh[,-1])

        sr  <- round(ss[,-1],0)
        sr[sr == nclass] <- NA
        chk <- bq[1,1] + bq[1,2]*hm
        dP2  <- -bq[sr,2]*exp(chk)/(1 + exp(chk))^2           #current state 2

        chk <- bq[2,1] + bq[2,2]*hm
        dP3  <- -bq[sr,2]*exp(chk)/(1 + exp(chk))^2           #current state 3


        ba  <- matrix(bg[xnames == 'AT'],ni,length(wt))
        if('ATXgap' %in% xnames){
           xtmp <- row2Mat(x3d[wi,wt,xnames == 'AT'])
           as1 <- ba*0
           as1[gap[wi] == 1,] <- bg[xnames == 'ATXgap']*xtmp[gap[wi] == 1,]
           ba <- ba + as1
        }
        if('ATXseed' %in% xnames){
           xtmp <- row2Mat(x3d[wi,wt,xnames == 'AT'])
           as1 <- ba*0
           as1[seed[wi] == 1,] <- bg[xnames == 'ATXseed']*xtmp[seed[wi] == 1,]
           ba <- ba + as1     
        }
        if('ATXchill' %in% xnames){
           xtmp <- row2Mat(x3d[wi,wt,xnames == 'ATXchill'])
           as1 <- ba*0
           as1 <- bg[xnames == 'ATXchill']*xtmp
           ba <- ba + as1     
        }

        sens2 <- dP2*dh*ba[,-1]
        tmp2  <- sens2
        tmp2[is.na(tmp2)] <- 0
        wy <- apply(tmp2,1,max,na.rm=T)
        wx <- apply(tmp2,1,which.max)
        wdd2 <- jd[1,wx[s2i[,1]]] - dt

        sens3 <- dP3*dh*ba[,-1]
        tmp3  <- sens3
        tmp3[is.na(tmp3)] <- 0
        wy <- apply(tmp3,1,max,na.rm=T)
        wx <- apply(tmp3,1,which.max)
        wdd3 <- jd[1,wx[s3i[,1]]] - dt

        s2i <- row2Mat(s2i[match(s3i[,1],s2i[,1]),])
        wdd2 <- wdd2[match(s3i[,1],s2i[,1])]

        bbDate  <- c(bbDate,jd[s3i])
        bbDD    <- c(bbDD,dd[s3i])
        maxSens2 <- c(maxSens2,wy[s2i[,1]])
        maxDate2 <- c(maxDate2,wdd2)
        maxSens3 <- c(maxSens3,wy[s3i[,1]])
        maxDate3 <- c(maxDate3,wdd3)
        gapDate <- c(gapDate,gap[s3i[,1]])
        seedDate <- c(seedDate,seed[s3i[,1]])
        chillDate <- c(chillDate,chill[s3i[,1]])
   #     siteDate <- c(siteDate,(match(site[s3i[,1]],siteCode) - 1))
        siteDate <- c(siteDate,site[s3i[,1]])
        chamDate <- c(chamDate,cham[s3i[,1]])
        yearDate <- c(yearDate,rep(y,length(s3i[,1])))
        firstbud <- c(firstbud,pf3[s3i[,1]])
        index    <- c(index,s3i[,1])
        trt      <- c(trt,treatments(site[s3i[,1]],cham[s3i[,1]]))

        hn <- as.vector(hh[,-1])
        xx <- as.vector(dd[,-1])
        yy <- as.vector(sens3)
        ss <- as.vector(ss[,-1])
        zz <- findInterval(ss,seq(0.01,5.99,length=100))

        wp <- which(hn < 99 & ss > 1.5 & is.finite(hn) & is.finite(ss) & yy > miny)

        gcol <- rep(1,n)
        if('gap' %in% xnames)gcol <- (gap-2)*-1

        points(xx[wp],yy[wp],col=colVals[zz[wp]],cex=.5)

        pd <- dd[s3i]
        ps <- sens3[s3i]

        wsd <- which(is.finite(pd) & is.finite(ps) & pd > 1 & ps > miny)
        if(length(wsd) > 0){
            symbols(pd[wsd],ps[wsd],circles=rep(.03,length=length(wsd)),
                    add=T,inches=F,bg=colVals[50],fg=gcol)
        }

        di <- si <- ji <- rep(NA,n)
        ji[s3i[,1]] <- jd[s3i]
        di[s3i[,1]] <- dd[s3i]
        si[s3i[,1]] <- sens3[s3i]

        dsout <- cbind(ji,di,si)
        colnames(dsout) <- paste(c('JD','DD','sens'),years[y],sep='-')
        DDdate <- cbind(DDdate,dsout)

     }
     lines(seq(1:5000),1/seq(1:5000),lty=2,lwd=2)

     stageLegend('bottomleft')

    DDvec <- as.vector(DDdate[,grep('DD',colnames(DDdate))])
    SSvec <- as.vector(DDdate[,grep('sens',colnames(DDdate))])

    tmp <- biVarMoments(log10(DDvec),log10(SSvec),1,PLOT=F)      #log10 scale
    dsMu <- tmp$mu
    dsVr <- tmp$var
    tmp  <- list(loc = dsMu, cov = dsVr, d2 = qchisq(.95,1) )
    tmp1  <- predict.ellipsoid(tmp)

    lines(10^tmp1[,1],10^tmp1[,2],type='l',lwd=2,lty=2,col=colVals[50])

    www <- which(bbDate > 25 & maxDate2 > 25 & bbDate < 200 & maxDate2 < 200)
    maxDate2 <- maxDate2[www]
    maxSens2 <- round(maxSens2[www],4)
    maxDate3 <- maxDate3[www]
    maxSens3 <- round(maxSens3[www],4)
    bbDate <- bbDate[www]
    bbDD   <- bbDD[www]
    seed   <- seedDate[www]
    gap    <- gapDate[www]
    chill  <- chillDate[www]
    site   <- siteDate[www]
    yr     <- yearDate[www]
    BB     <- firstbud[www]
    index  <- index[www]
    trt    <- match(trt[www],c('amb','+3','+5','ctl'))
    trt[trt == 4] <- 1
    ddwinter <- ddwinter[www]

    i <- index
    maxDate <- cbind(yr,i,site,trt,seed,gap,bbDD,maxDate2,maxSens2,maxDate3,maxSens3,bbDate,BB,chill,ddwinter)
    rownames(maxDate) <- paste(gibbsSpec,chamDate[www],sep='-')

    list(ddssMat = DDdate, dsMu = dsMu, dsVr = dsVr, ellipse = tmp, maxDate = maxDate,
         firstbud = firstbud)
}
###################################################################
stageLegend <- function(location='topright',back=NULL){

 # colF <-  colorRampPalette(c("blue","orange","red"))

#  colF <-  colorRampPalette(c("brown","green","darkgreen"))
  colVals <- colF(nclass)[c(1,3,6)]

  ltext <- c('dormant','budbreak','expanded')

  if(is.null(back)) legend(location,ltext,text.col=colVals,bty='n')
  if(!is.null(back))legend(location,ltext,text.col=colVals,bg=back,box.col=back)
}

########################################################3
plotstart <- function(plotfile){

  #initiates plot
   dev.off()
#  postscript(file=plotfile,width=6, height=9,horizontal=FALSE)

}

####################################################

plotend <- function(plotfile){

  #terminates plot

#  if (REMOTE)dev.off()
  dev.print(device=postscript,file=plotfile,width=7, height=10, horizontal=FALSE)

}


getDelta <- function(xy,treat,siten,year){

  tcode <- c('+3','+5','amb')

  sc  <- match(siten,c('DF','HF'))
  tc  <- match(treat,tcode)

  muDF <- apply(meanDF[1:130,],2,mean)   #spring mean values
  muHF <- apply(meanHF[1:130,],2,mean)

  names(muDF)[names(muDF) == 'p3'] <- '+3'
  names(muDF)[names(muDF) == 'p5'] <- '+5'
  names(muHF)[names(muHF) == 'p3'] <- '+3'
  names(muHF)[names(muHF) == 'p5'] <- '+5'

  mmj <- muDF

  temp <- rep(0,length(treat))

  for(y in 1:nyr){
    ytemp <- annualAve[y,]
    temp[siten == 'DF' & treat == 'amb' & year == y] <- ytemp['DF-amb']
    temp[siten == 'DF' & treat == '+3' & year == y]  <- ytemp['DF-amb'] + 3
    temp[siten == 'DF' & treat == '+5' & year == y]  <- ytemp['DF-amb'] + 5
    temp[siten == 'HF' & treat == 'amb' & year == y] <- ytemp['HF-amb']
    temp[siten == 'HF' & treat == '+3' & year == y]  <- ytemp['HF-amb'] + 3
    temp[siten == 'HF' & treat == '+5' & year == y]  <- ytemp['HF-amb'] + 5
  }
  muT <- tapply(temp,list(siten,treat,year),mean)
  muD <- tapply(xy[,2],list(siten,treat,year),mean)

  wwt <- which(is.na(muT),arr.ind=T)
  if(length(wwt) > 0){
    for(kw in 1:nrow(wwt)){
       w <- wwt[kw,]
       if(w[1] == 1 & w[2] == 1)muT[w[1],w[2],w[3]] <- annualAve[w[3]+1,'DF-p3']
       if(w[1] == 2 & w[2] == 1)muT[w[1],w[2],w[3]] <- annualAve[w[3]+1,'HF-p3']
       if(w[1] == 1 & w[2] == 2)muT[w[1],w[2],w[3]] <- annualAve[w[3]+1,'DF-p5']
       if(w[1] == 2 & w[2] == 2)muT[w[1],w[2],w[3]] <- annualAve[w[3]+1,'HF-p5']
       if(w[1] == 1 & w[2] == 3)muT[w[1],w[2],w[3]] <- annualAve[w[3]+1,'DF-amb']
       if(w[1] == 2 & w[2] == 3)muT[w[1],w[2],w[3]] <- annualAve[w[3]+1,'HF-amb']
    }
  }

  ndim <- dim(muT)

  dD <- dT <- muT*0

  for(j in 1:ndim[1]){         #departures
    for(y in 1:ndim[3]){
      for(k in 1:ndim[2]){
        dD[j,k,y] <- muD[j,k,y] - muD[j,'amb',y]
        dT[j,k,y] <- muT[j,k,y] - muT[j,'amb',y]
      }
   }
  }

  muR <- vrR <- matrix(NA,ndim[1],ndim[2])
  rownames(muR) <- rownames(vrR) <- sort(unique(siten))
  colnames(muR) <- colnames(vrR) <- c('+3','+5','amb')
  

  for(j in 1:ndim[1]){ 
     for(k in 1:3){
         r        <- dD[j,k,]/dT[j,k,]
         muR[j,k] <- mean(r,na.rm=T)
         njk <- sum(r*0 + 1,na.rm=T)
         v   <- var(dD[j,k,],na.rm=T)
         vrR[j,k] <- max(1/dT[j,k,]^2*v/njk,na.rm=T)
       }
  }

  srR <- sqrt(vrR)

  if(nrow(muR) == 1){
       newRow <- rep(NA,3)
       if(rownames(muR)[1] == 'DF'){
            HF   <- newRow
            muR <- rbind(muR,HF)
            srR <- rbind(srR,HF)
       }
       if(rownames(muR)[1] == 'HF'){
            DF   <- newRow
            muR <- rbind(DF,muR)
            srR <- rbind(DF,srR)
       }
    }

    en <- outer(rownames(muR),colnames(muR),paste,sep='')
    ev <- as.vector(muR)
    sv <- as.vector(srR)
    names(ev) <- names(sv) <- en

  list(erat = ev, srat = sv)
}

get95 <- function(mvec,svec){   #vector of means and sds

  lo <- mvec - 1.96*svec
  hi <- mvec + 1.96*svec

  rbind(mvec,lo,hi)
}


getOutFiles <- function(path,NEWFILE){

  if(NEWFILE){
    filekkk <- list.files(path=pathkk,pattern='_gibbs_')
    if(length(filekkk) == 0)return(character(0))
    fk <- filekkk
    fk <- sub('.RData','',fk)
    tmp <- matrix(unlist(strsplit(fk,'_')),ncol=3,byrow=T)
    fspec <- tmp[,1]
    fng   <- as.numeric(tmp[,3])
    fss   <- unique(fspec)
    fkeep    <- character(0)
    for(j in 1:length(fss)){
      ws <- which(fspec == fss[j])
      fkeep <- c(fkeep,filekkk[ws[which.max(fng[ws])]])
    }
    fdel <- filekkk[!filekkk %in% fkeep]
    if(length(fdel) > 0)file.remove(paste(pathkk,fdel,sep='/'))
    filekkk <- fkeep
  }

  if(!NEWFILE){
   filekkk <- list.files(path=pathkk,pattern='RData')
   if(length(filekkk) == 0)return(character(0))
  }

  filekkk
}
############################################################
plotDDbyTrt <- function(){

  ns <- length(specList)

  dd1 <- dd2 <- doy1 <- doy2 <- deldd <- deldoy <- perDegree <- numeric(0)
  ee1x <- ee1y <- ee2x <- ee2y <- numeric(0)
  cr1  <- cr2 <- numeric(0)

  for(j in 1:ns){

     wj   <- grep(specList[j],rownames(maxSens))
     if(length(wj) < 120)next
     dd   <- maxSens[wj,'bbDD']
     doy  <- maxSens[wj,'BB']
     site <- maxSens[wj,'site']
     
     trt  <- maxSens[wj,'trt']
     deg  <- (trt - 1)*2 + 1
     deg[deg == 1] <- 0
     trt[trt == 3] <- 2
     seed <- maxSens[wj,'seed']

     tmp  <- siteByTrt(dd,site,trt,specList[j],maxx = 700)
     dd1  <- rbind(dd1,tmp$mu[[1]])
     dd2  <- rbind(dd2,tmp$mu[[2]])
     deldd <- rbind(deldd,tmp$delta[[2]])

     tmp  <- siteByTrt(doy,site,trt,specList[j],maxx = 150)
     doy1 <- rbind(doy1,tmp$mu[[1]])
     doy2 <- rbind(doy2,tmp$mu[[2]])
     deldoy <- rbind(deldoy,tmp$delta[[2]])

     tmp  <- siteByTrt(doy,site,deg,specList[j],maxx = 150)
     perDegree <- rbind(perDegree,tmp$perDegree)

     doyy <- doy
     doyy[doyy > 150] <- NA

     w1 <- which(site == 1)
     if(length(w1) > 20){
       tmp <- biVarMoments(doyy[w1],deg[w1],wt = 1)
       cr  <- tmp$var[2]/sqrt(tmp$var[1]*tmp$var[4])
       cr1 <- rowBind(cr1,c(tmp$mu,cr),specList[j])
       
       tmp  <- list(loc = tmp$mu, cov = tmp$var, d2 = qchisq(.50,1) )
       tmp1 <- predict.ellipsoid(tmp)
       ee1x <- rowBind(ee1x,tmp1[,1],specList[j])
       ee1y <- rowBind(ee1y,tmp1[,2],specList[j])
     }
     w2 <- which(site == 2)
     if(length(w2) > 20){
       tmp <- biVarMoments(doyy[w2],deg[w2],wt = 1)
       cr  <- tmp$var[2]/sqrt(tmp$var[1]*tmp$var[4])
       cr2 <- rowBind(cr2,c(tmp$mu,cr),specList[j])
       tmp  <- list(loc = tmp$mu, cov = tmp$var, d2 = qchisq(.50,1) )
       tmp1 <- predict.ellipsoid(tmp)
       ee2x <- rowBind(ee2x,tmp1[,1],specList[j])
       ee2y <- rowBind(ee2y,tmp1[,2],specList[j])
     }

  }

  par(mfcol=c(2,2),bty='n',plt=c(.2,1,.3,.95))
  tmp <- matrix(unlist(strsplit(rownames(dd1),'-')),ncol=2,byrow=T)[,1]
  scol <- c('black','blue')
  
  ylab <- expression(paste('Advance per',~degree~C,sep=''))

  plotSetup(c(0,200,400,600),c(-4,-2,0,2),xlabel=' ',ylabel= ylab,
            endFactor=c(.05,.15),fc='wheat')
  for(j in 1:nrow(dd1)){
     if(length(grep('1',rownames(dd1)[j])) > 0){  #site 1
        dj <- dd1[j,]
        sc <- scol[1]
        sp <- 4
        xt <- dj[1] + max( c(dj[2],0) ,na.rm=T)
     }
     if(length(grep('2',rownames(dd2)[j])) > 0){
        dj <- dd1[j,]
        sc <- scol[2]
        sp <- 2
        xt <- dj[1] - max( c(dj[2],0) ,na.rm=T)
     }
     aj <- perDegree[j,]
     points(dj[1],aj[1],pch=3,col=sc)
     lines(c(dj[1],dj[1]),c(aj[1] - aj['se'],aj[1] + aj['se']),col=sc,lwd=2 )
     lines(c(dj[1] - dj[2],dj[1] + dj[2]),c(aj[1],aj[1]), col=sc,lwd=2)
     text(xt,aj[1],tmp[j],col=sc,pos=sp)
  }

  par(bty='n',plt=c(.2,.6,.3,.95))
  #ellipse response
  ylab <- expression(paste('Proportionate advance per',~degree~C,sep=''))
  plotSetup(c(0,1,2,3,4,5),c(80,100,120,140),ylabel='DOY',
             xlabel=expression(~degree~C),endFactor=c(.05,.15),fc='wheat')
  for(j in 1:nrow(ee1x))lines(ee1y[j,],ee1x[j,],lwd=2)
  for(j in 1:nrow(ee2x))lines(ee2y[j,],ee2x[j,],lwd=2,col='blue')

  plotSetup(c(0,1),c(80,100,120,140),ylabel='DOY',
             xlabel='Correlation',endFactor=c(.05,.15),fc='wheat')


  sscol <- rep(1,nrow(dd1))
  sscol[grep('2',rownames(dd1))] <- 2
  ps <- sscol*0
  ps[sscol == 1] <- 4
  ps[sscol == 2] <- 2
  points(dd1[,1],perDegree[,1]/doy1[,1],col=scol[sscol])
  text(dd1[,1],perDegree[,1]/doy1[,1],tmp,col=scol[sscol],pos=ps)

}



siteByTrt <- function(xx,site,trt,spec,maxx = Inf){

     xx[xx > maxx] <- NA

     mu  <- by(xx,list(site = site, trt = trt),mean,na.rm=T)
     rn  <- paste('site',rownames(mu),sep='_')
     cn  <- paste('trt',colnames(mu),sep='_')
     mu  <- matrix(mu,nrow=length(rownames(mu)))
     rownames(mu) <- rn
     colnames(mu) <- cn

     sdv <- by(xx,list(site = site, trt = trt),sd,na.rm=T)
     sdv <- matrix(sdv,nrow=length(rownames(sdv)))
     rownames(sdv) <- rn
     colnames(sdv) <- cn

     n  <- table(site,trt)
     se <- sdv/sqrt(n)
     n1 <- n2 <- 0
     if('1' %in% rownames(n))n1 <- sum(n[rownames(n) == '1',])
     if('2' %in% rownames(n))n2 <- sum(n[rownames(n) == '2',])

     muD <- sdD <- seD <- matrix(0,nrow(mu),1)
     rownames(muD) <- rownames(sdD) <- rownames(seD) <- rownames(mu)

     w0  <- which(trt > 0)
     w1  <- which(trt == 0)

     if(length(w0) > 0 & length(w1) > 0){
       m0   <- by(xx[-w0],list(site = site[-w0]),mean,na.rm=T)
       mm   <- rep(0,length(trt))
       if('1' %in% rownames(m0))mm[site == 1] <- m0[rownames(m0) == '1']
       if('2' %in% rownames(m0))mm[site == 2] <- m0[rownames(m0) == '2']
      
       muD <- by( (xx - mm)[w0]/trt[w0],list(site = site[w0]),mean,na.rm=T)
       rn   <- paste('site',rownames(muD),sep='_')
       muD  <- matrix(muD,nrow=length(rownames(muD)))
       rownames(muD) <- rn

       sdD  <- by( (xx - mm)[w0]/trt[w0],list(site = site[w0]),sd,na.rm=T)
       seD  <- sdD/table(site[w0])
       sdD  <- matrix(sdD,nrow=length(rownames(sdD)))
       rownames(sdD) <- rn
       seD  <- matrix(seD,nrow=length(rownames(seD)))
       rownames(seD) <- rn
     }

     if(length(grep('1',rownames(mu))) == 0){
       v <- matrix(NA,1,ncol(mu),dimnames=list('site_1',colnames(mu)))
       mu  <- rbind(v,mu)
       sdv <- rbind(v,sdv)
       se  <- rbind(v,se)
     }
    if(length(grep('2',rownames(mu))) == 0){
       v <- matrix(NA,1,ncol(mu),dimnames=list('site_1',colnames(mu)))
       mu  <- rbind(mu,v)
       sdv <- rbind(sdv,v)
       se  <- rbind(se,v)
     }
     if(length(grep('1',rownames(muD))) == 0){
       muD <- rbind(NA,muD)
       sdD <- rbind(NA,sdD)
       seD <- rbind(NA,seD)
       rownames(muD)[1] <- rownames(sdD)[1] <- rownames(seD)[1] <- 'site_1'
    }
    if(length(grep('2',rownames(mu))) == 0){
       muD <- rbind(muD,NA)
       sdD <- rbind(sdD,NA)
       seD <- rbind(seD,NA)
       rownames(muD)[2] <- rownames(sdD)[2] <- rownames(seD)[2] <- 'site_2'
    }


     if(n1 < 5)muD[1,] <- sdD[1,] <- seD[1,] <- NA
     if(n2 < 5)muD[2,] <- sdD[2,] <- seD[2,] <- NA

    perDegree <- cbind(muD,sdD,seD)
    colnames(perDegree) <-  c('mu','sd','se')
    rownames(perDegree) <- paste(spec,rownames(perDegree),sep='-')

     kt <- c(0,3,5)
     slist <- dlist <- numeric(0)
     for(k in 1:ncol(mu)){
       s1 <- cbind(mu[,k],sdv[,k],se[,k])
       rownames(s1) <- paste(spec,rownames(s1),sep='-')
       colnames(s1) <- c('mu','sd','se')
       kn <- colnames(mu)[k]
       slist <- append(slist,list(s1))
   #    names(slist)[k] <- kn

       dmu <- (s1[,1] - slist[[1]][,1])/kt[k]
       dsd <- sqrt(s1[,2]^2 + slist[[1]][,2]^2)/kt[k]
       dse <- sqrt(s1[,3]^2 + slist[[1]][,3]^2)/kt[k]
       delta <- cbind(dmu,dsd,dse)
       colnames(delta) <- colnames(s1)
       rownames(delta) <- rownames(s1)

       dlist <- append(dlist,list(delta))
  #     names(dlist)[k] <- kn
     }

     list(mu = slist, dlist = dlist, perDegree = perDegree, n = n)
}








