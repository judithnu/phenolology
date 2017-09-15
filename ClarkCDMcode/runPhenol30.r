
#sftp jimclark@dscr-login-06.oit.duke.edu
#ssh jimclark@dscr-login-06.oit.duke.edu

#qfile:
#qsub runBB highprio


rm(list=ls())


INITALL <- F     #if F, read in previous output

mergeACSA <- T       #merge sugar maples?
lastYr    <- 2012

models  <- c(2:3)


ngg    <- 15000
burn0  <- 8000


load('bbModel.RData')
source('bbFunctions.r')
source('functionLibrary.r')

vars      <- c('int','AT','chill','ATXchill','seed','ATXseed','gap','ATXgap')

modelList <- vrlist <- character(0)
for(k in 2:length(vars)){
  modelList <- append(modelList,list(vars[1:k]))
}
var2 <- c('int','AT','gap','ATXgap','chill')
for(k in 3:length(var2)){
   modelList <- append(modelList,list(var2[1:k]))
}



nmodel <- length(modelList)


spec2mod <- c("acru","acsa","beal","bepa","fram","list","litu","magr","nysy",
              "pipa","pire","pist","pita","qual","quru")

#spec2mod <- 'nysy'

nclass <- 6

kkk <- 1
jjj <- 2

dicModel <- matrix(NA,length(specList),(length(vars)-1))
rownames(dicModel) <- specList
colnames(dicModel) <- vars[-1]


for(kkk in models){

  for(jjj in 1:length(specList)){

   INIT <- INITALL
   ng   <- ngg

   gibbsSpec <- specList[jjj]

   if(!gibbsSpec %in% spec2mod)next

   xnames    <- modelList[[kkk]]
   nx        <- length(xnames)

   totalg <- 0

   pathName <- 'model'
   for(k in 2:nx)pathName <- paste(pathName,xnames[k],sep='-')

   if(!pathName %in% list.files())dir.create(pathName)
   fname    <- paste(gibbsSpec,'gibbs',sep='_')
   fullname <- paste(pathName,'/',fname,sep='')

   if(!INIT){
      files <- list.files(path=pathName,pattern=fname)
      kfile <- grep(gibbsSpec,files)
      if(length(kfile) == 0)INIT <- T

      if(length(kfile) > 1){
        fk <- files[kfile]
        fs <- unlist(strsplit(fk,'.RData'))
        fs <- as.numeric(matrix(unlist(strsplit(fs,'_')),ncol=3,byrow=T)[,3])
        wm <- which.max(fs)
        wl <- kfile[-wm]

        fkeep <- fk[wm]
        infile <- paste(pathName,fkeep,sep='/')
        file.remove(paste(pathName,fk[wl],sep='/'))
        kfile <- kfile[wm]
     }

      if(length(kfile) == 1){
        fk <- files[kfile]
        infile <- paste(pathName,fk,sep='/')
        load(infile)
        INIT <- F
      }
   }

   ff <- paste(gibbsSpec,'.Rdata',sep='')
   if(!INIT)save(sstate, x3d, multiLogitStates, years, yrBreak,yrIndex,yrVec,sprCols,springDates,sampleDates,sampleIndex,hIndex,xIndex,jdIndex,jdAll,germTime,fallDates,dtvec,diedTime,chambers,ddmat,phenData,meanX, file=ff)

  }

}


   modelName <- modelList[[kkk]]

   burnin <- burn0

   if(INIT)source('phenInit.r')


   source('phenolGibbs.r')

   dicModel[jjj,kkk] <- round(dic,0)

   newname <- paste(fullname,'_',totalg,'.RData',sep='')

   stuff2save <- ls()[(!ls() %in% c('jjj','kkk','ngg','ng','INIT'))]


   save(list = stuff2save, file = newname)

   print(newname)
   print(bg)





