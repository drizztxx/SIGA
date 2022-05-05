SIGA <- function(ds,
                 trt,
                 ctrl,
                 cutoff = 10,
                 varScreen = c(2,10),
                 resptyp="tte",
                 clnum=5,
                 time,
                 status,
                 objFunction,
                 popSize = 50,
                 MAXGEN = 50,
                 n.boot=10,
                 bootres=NULL,
                 freq=0.5,
                 multiplicity=0.25,
                 CPUs=1){
  siga.argslist <- list(trt = trt,cutoff = cutoff,varScreen = varScreen,resptyp=resptyp,clnum=clnum,
                        ctrl = ctrl,time=time,status=status,objFunction = objFunction,popSize = popSize,MAXGEN = MAXGEN)
  ## Get original dataset
  newdata <- ds
  trtvar <- strsplit(trt,'=')[[1]][1]
  trtval <- strsplit(trt,"=")[[1]][2]
  ctrlval <- strsplit(ctrl,"=")[[1]][2]
  newdata[,trtvar] <- factor(newdata[,trtvar],levels=c(ctrlval,trtval))
  ## Calculate overall hr
  overall.fit <- coxph(as.formula(sprintf("Surv(%s, %s) ~ %s",
                                          time,status,trtvar)),data=newdata)
  overall.hr <- summary(overall.fit)$coefficients[2]
  
  ## Apply resampling
  if (is.null(bootres)){
    bootres <- boot.SIGA(data.boot = ds,method="parSIGA",n.boot = n.boot,CPUs = CPUs,
                         methodargslist = siga.argslist)
  }
  vars <- freqvars(bootres,freq=freq)
  meanhr <- mean(bootres$testd,na.rm = T)
  if (meanhr >= 1 || meanhr > overall.hr*(1-multiplicity) || length(vars) == 0)
    best.sub="no subgroup"
  else {
    siga.argslist$data <- str2lang("ds")
    siga.argslist$varused <- vars
    sigares <- do.call("parSIGA",siga.argslist)
    Best <- Best.SIGA(sigares)
    bestsub <- Best$best.subgroup
    best.sub=bestsub
  }
  dsn <- deparse(substitute(ds))
  res <- list(bestsub=best.sub,bootres=bootres,overall.hr=overall.hr,meanhr=meanhr,
              varused=vars,dsn=dsn,call=match.call())
  class(res) <- "SIGA"
  return(res)
}