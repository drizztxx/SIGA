nullsimul.siga <- function(simul.times,
                           exp,
                           CPUs,
                           clnum=5,
                           objF="objF.logrank",
                           popSize=50,
                           n.folds=3,
                           MAXGEN=50,
                           varScreen=c(2,5)
                           ){
  library(parallel)
  if (clnum*1*CPUs > detectCores()) 
    stop(sprintf("Core numbers %s >= detect max numbers %s",clnum*1*CPUs,detectCores()))
  cl3 <- makeCluster(getOption("cl.cores", CPUs))
  clusterExport(cl3, c("exp","clnum","objF","popSize","MAXGEN","n.folds","varScreen"),envir = environment())
  aa <- clusterEvalQ(cl3,source("~/init.R"))
  aa <- clusterEvalQ(cl3,source("~/simul.Cox.R"))
  aa <- clusterEvalQ(cl3,library(survival))
  res <- parLapply(cl3,1:simul.times,function(s){
    testdata <- simul.Cox(n=400,exp=exp, 
              weib.shape = 4.34,weib.scale =1.2,censor.rate=0.1,seed=s)
    cox1 <- coxph(Surv(Ts,event)~z,data=testdata)
    sum<-summary(cox1)
    overall.hr <- sum$coefficients[2]
    cv<-CrossValidation(data.cv=testdata,method="parSIGA",n.folds = n.folds,CPUs = 1,
      methodargslist = list(trt = "z=1",cutoff = 10,varScreen = varScreen,
      resptyp="tte",clnum=clnum,
      ctrl = "z=0",time="Ts",status="event",
     objFunction = objF,popSize = popSize,MAXGEN = MAXGEN,verbose = F))
    return(list(cv=cv,hr=overall.hr))
  })
  stopCluster(cl3)
  return(res)
}