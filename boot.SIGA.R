boot.SIGA <- function(data.boot,
                            method = c("parSIGA","SIGA.tte","SIGA.binary"),
                            methodargslist,
                            n.boot = 10,
                            seed = 123,
                            CPUs = 1){
  method = match.arg(method)
  methodargslist$data <- str2lang("data.train")
  if (method == "parSIGA")
    outcome.type <- methodargslist$resptyp
  else if (method == "SIGA.tte")
    outcome.type <- "tte"
  if (method %in% c("parSIGA","SIGA.tte","SIGA.binary") & CPUs > 1) {
    cl2 <- makeCluster(getOption("cl.cores", min(n.boot,CPUs)))
    clusterExport(cl2,c("SIGA.binary","SIGA.tte","parSIGA", "objF.logrank","objF.logrank2",
                        methodargslist$objFunction,"objF.cox","decodsubset","predict.SIGA","Best.SIGA"))
    clusterExport(cl2,c("seed", "data.boot","method","methodargslist"),envir=environment())
    aa <- clusterEvalQ(cl2,library(randomForest))
    aa <- clusterEvalQ(cl2,library(randomForestSRC))
    aa <- clusterEvalQ(cl2,library(gatbxr))
    aa <- clusterEvalQ(cl2,library(survival))
    system.time(cv.list <- parLapply(cl2,1:n.boot,function(i){
      set.seed(seed+i)
      index <- as.logical(rbinom(nrow(data.boot),1,0.632))
      data.train <- data.boot[index,]
      data.test <- data.boot[!index,]
      try1 <- try(trainres <- do.call(method,args=methodargslist))
      if (!inherits(try1,"try-error")){
        testres <- predict(trainres,data.test)
        return(list(trainres=trainres,testres=testres,vars=trainres$varused,objv=min(trainres$Best)))
      } else {
        return(NULL)
      }
    }))
    stopCluster(cl2)
  } else if (method %in% c("parSIGA","SIGA.tte","SIGA.binary") & CPUs == 1) {
    system.time(cv.list <- lapply(1:n.boot,function(i){
      set.seed(seed+i)
      index <- as.logical(rbinom(nrow(data.boot),1,0.632))
      data.train <- data.boot[index,]
      data.test <- data.boot[!index,]
      try1 <- try(trainres <- do.call(method,args=methodargslist))
      if (!inherits(try1,"try-error")){
        testres <- predict(trainres,data.test)
        return(list(trainres=trainres,testres=testres,vars=trainres$varused,objv=min(trainres$Best)))
      } else {
        return(NULL)
      }
    }))
  }
  cv.list <- cv.list[!unlist(lapply(cv.list,is.null))]
  cv.subs <- unlist(lapply(cv.list,function(x) x$testres$bestsub))
  cv.trainp <- unlist(lapply(cv.list,function(x) x$testres$trainobjv$objv))
  cv.traind <- unlist(lapply(cv.list,function(x) x$testres$trainobjv$diff))
  cv.trainprop <- unlist(lapply(cv.list,function(x) x$testres$trainobjv$sub.proption))
  cv.testp <- unlist(lapply(cv.list,function(x) x$testres$testobjv$objv))
  cv.testd <- unlist(lapply(cv.list,function(x) x$testres$testobjv$diff))
  if (outcome.type %in% c("tte","survival","s") ) cv.testd[cv.testd > 100 | cv.testd < 0.001] <- NA ## reset too large HR value to NA
  cv.testprop <- unlist(lapply(cv.list,function(x) x$testres$testobjv$sub.proption))
  if (method %in% c("parSIGA","SIGA.tte","SIGA.binary")) {
    vars <- apply(do.call(rbind,lapply(cv.list,function(x) x$vars)),1,paste,collapse=" ")
    objv <- unlist(lapply(cv.list,function(x) x$objv))
    res <- data.frame(subs=cv.subs,trainp=cv.trainp,testp=cv.testp,traind=cv.traind,testd=cv.testd,trainprop=cv.trainprop,testprop=cv.testprop,vars=vars,objv=objv)
  } else res <- data.frame(subs=cv.subs,trainp=cv.trainp,testp=cv.testp,traind=cv.traind,testd=cv.testd,trainprop=cv.trainprop,testprop=cv.testprop)
  return(res)
}
