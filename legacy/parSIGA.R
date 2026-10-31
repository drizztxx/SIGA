parSIGA <- function(data,
                    trt="treat=='E'",
                    ctrl="treat=='C'",
                    resp="res==1",
                    status="status",
                    time="time",
                    varScreen=c(1,5),
                    varused=NULL,
                    popSize=100,
                    resptyp=c("binary","tte"),
                    clnum=1,
                    cutoff=10,
                    objFunction=c("objF.glm2","objF.glm","objF.difference","objF.difference2","objF.cox","objF.cox2","objF.logrank","objF.logrank2"),
                    MAXGEN = 400,        ## max Number of generations
                    GGAP = 0.9,          ## Generation gap, how many new individuals are created
                    SEL_F = 'sus',       ## Name of selection function
                    XOV_F = 'reclin',     ## Name of recombination function for individuals
                    MUT_F = 'mutbga',       ## Name of mutation function for individuals
                    ...
){
  # data=er.data;trt = "treat=1";cutoff = 10;clnum=1;resptyp ="tte";varScreen = c(1,5)
  #  ctrl = "treat=2";time="t.rfs";status="e.rfs";objFunction = "objF.cox2";MAXGEN = 1
  #  GGAP = 0.9       ## Generation gap, how many new individuals are created
  #  SEL_F = 'sus'    ## Name of selection function
  #  XOV_F = 'reclin'    ## Name of recombination function for individuals
  #  MUT_F = 'mutbga'     ## Name of mutation function for individuals
  #  varused=NULL
  #  popSize=100
  library(parallel)
  library(randomForest)
  library(randomForestSRC)
  library(gatbxr)
  library(survival)
  resptyp <- match.arg(resptyp)
  if (clnum > detectCores()) 
    stop(sprintf("Core numbers %s >= detect max numbers %s",clnum,detectCores()))
  cl <- makeCluster(getOption("cl.cores", clnum))
  if (length(deparse(substitute(data))) > 1)
    stop("Data cannot be deparsed")
  if (!exists(deparse(substitute(data))))
    assign(deparse(substitute(data)),dynGet(deparse(substitute(data))))
  clusterExport(cl,c("SIGA.binary","SIGA.tte",objFunction,"decodsubset"))
  clusterExport(cl, deparse(substitute(data)),envir = environment())
  aa <- clusterEvalQ(cl,library(randomForest))
  aa <- clusterEvalQ(cl,library(randomForestSRC))
  aa <- clusterEvalQ(cl,library(gatbxr))
  aa <- clusterEvalQ(cl,library(survival))
  if (resptyp == "binary"){
    trtvar <- strsplit(trt,'=')[[1]][1]
    resvar <- strsplit(resp,'=')[[1]][1]
    if (!is.null(varused)){
      col.vim <- varused
    } else if (varScreen[1]) { ## Screen the top important variables
      top <- varScreen[2]
      screen.res <- VariableScreen(data=data,trt=trt,ctrl=ctrl,resp=resp,varScreen = varScreen,resptyp = "binary")
      col.vim <- screen.res$col.vim[1:top]
    } else {
      top <- NCOL(data) - 2
      col.vim <- colnames(data)[!colnames(data) %in% c(trtvar,resvar)]
    }
    res <- system.time(pool <- parLapplyLB(cl,1:clnum, SIGA.binary,data=data,trt = trt,ctrl = ctrl, resp = resp, varused = col.vim, popSize=popSize, GGAP=GGAP,
                                         SEL_F=SEL_F,XOV_F=XOV_F,MUT_F=MUT_F,
                              objFunction = objFunction,MAXGEN = MAXGEN,verbose = FALSE,cutoff=cutoff,family="binomial"))
  } else if (resptyp == "tte"){
    trtvar <- strsplit(trt,'=')[[1]][1]
    if (!is.null(varused)){
      col.vim <- varused
    } else if (varScreen[1]) { ## Screen the top important variables
      top <- varScreen[2]
      screen.res <- VariableScreen(data=data,trt=trt,ctrl=ctrl,status = status,time=time,varScreen = varScreen,resptyp = "tte")
      col.vim <- screen.res$col.vim[1:top]
    } else {
      top <- NCOL(data) - 3
      col.vim <- colnames(data)[!colnames(data) %in% c(trtvar,time,status)]
    }
    res <- system.time(pool <- parLapplyLB(cl,1:clnum, SIGA.tte,data=data,trt = trt,ctrl = ctrl, status = status, time = time, varused=col.vim, popSize=popSize, GGAP=GGAP,
                              SEL_F=SEL_F,XOV_F=XOV_F,MUT_F=MUT_F,cutoff=cutoff,
                              objFunction = objFunction,MAXGEN = MAXGEN,verbose = FALSE))
  }
  stopCluster(cl)
  Chrom <- lapply(1:MAXGEN, function(x) {do.call(rbind, lapply(pool, function(y) y$Chrom[[x]]))})
  Objv <- lapply(1:MAXGEN, function(x) {unlist(lapply(pool, function(y) y$Objv[[x]]))})
  Subset <- lapply(1:MAXGEN, function(x) {unlist(lapply(pool, function(y) y$Subset[[x]]))})
  Best <- apply(do.call(rbind,lapply(pool, function(x) x$Best)),2,min)
  elapsed <- res[3]
  dsn <- deparse(substitute(data))
  if (!exists("screen.res")) screen.res <- NULL
  final <- list(Chrom=Chrom,Objv=Objv,Subset=Subset,Best=Best,varused=col.vim,dsn=dsn,elapsed=elapsed,
                screen.res=screen.res,call=match.call())
  class(final) <- "SIGA"
  return(final)
}
