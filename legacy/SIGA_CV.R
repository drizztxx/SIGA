CrossValidation <- function(data.cv = null,
                            method = c("parSIGA","SIGA.tte","SIGA.binary","SIGA","SIDES","PRIM"),
                            methodargslist,
                            n.folds = 10,
                            seed = 123,
                            CPUs = 1,
                            ...){
  method = match.arg(method)
  if (method == "SIGA")
    methodargslist$ds <- str2lang("data.train")
  else methodargslist$data <- str2lang("data.train")
  if(is.null(data.cv)) stop("data.cv cannot be null")
  n <- NROW(data.cv)
  idx_fix <- rep_len(1:n.folds,n)
  set.seed(seed)
  idx <- sample(idx_fix)
  if (method %in% c("parSIGA","SIGA"))
    outcome.type <- methodargslist$resptyp
  else if (method == "SIGA.tte")
    outcome.type <- "tte"
  else if (method == "SIDES")
    outcome.type <- methodargslist$data_set_parameters$outcome_variable_type
  else if (method == "PRIM")
    outcome.type <- methodargslist$type
  
  if (method %in% c("parSIGA","SIGA.tte","SIGA.binary","SIGA") & CPUs > 1) {
    cl2 <- makeCluster(getOption("cl.cores", min(n.folds,CPUs)))
    clusterExport(cl2,c("SIGA.binary","SIGA.tte","parSIGA", "objF.logrank","objF.logrank2","SIGA","boot.SIGA","freqvars",
                        methodargslist$objFunction,"objF.cox","decodsubset","predict.SIGA","Best.SIGA"))
    clusterExport(cl2,c("idx", "data.cv","method","methodargslist"),envir=environment())
    if (...length() > 0){
      for (ssss in 1:...length()){
        assign(...names()[ssss],...elt(ssss))
      }
      clusterExport(cl2,...names(),envir=environment())
    }
    aa <- clusterEvalQ(cl2,library(randomForest))
    aa <- clusterEvalQ(cl2,library(randomForestSRC))
    aa <- clusterEvalQ(cl2,library(gatbxr))
    aa <- clusterEvalQ(cl2,library(survival))
    system.time(cv.list <- parLapply(cl2,1:n.folds,function(i){
      data.train <- data.cv[!idx==i,]
      data.test <- data.cv[idx==i,]
      trainres <- do.call(method,args=methodargslist)
      testres <- predict(trainres,data.test)
      if (method == "SIGA")
        return(list(trainres=trainres,testres=testres,vars=trainres$varused,objv=NA,bootres=trainres$bootres))
      else return(list(trainres=trainres,testres=testres,vars=trainres$varused,objv=min(trainres$Best)))
    }))
    stopCluster(cl2)
  } else if (method %in% c("parSIGA","SIGA.tte","SIGA.binary","SIGA") & CPUs == 1) {
    if (...length() > 0){
      for (ssss in 1:...length()){
        assign(...names()[ssss],...elt(ssss))
      }
    }
    system.time(cv.list <- lapply(1:n.folds,function(i){
      data.train <- data.cv[!idx==i,]
      data.test <- data.cv[idx==i,]
      trainres <- do.call(method,args=methodargslist)
      testres <- predict(trainres,data.test)
      if (method == "SIGA")
        return(list(trainres=trainres,testres=testres,vars=trainres$varused,objv=NA,bootres=trainres$bootres))
      else return(list(trainres=trainres,testres=testres,vars=trainres$varused,objv=min(trainres$Best)))
    }))
  } else if (method == "SIDES") {
    SIDES.trt <- sprintf("%s=%s",methodargslist$data_set_parameters$treatment_variable_name,
                         unique(data.cv[,methodargslist$data_set_parameters$treatment_variable_name])[unique(data.cv[,methodargslist$data_set_parameters$treatment_variable_name]) != methodargslist$data_set_parameters$treatment_variable_control_value])
    SIDES.ctrl <- sprintf("%s=%s",methodargslist$data_set_parameters$treatment_variable_name,
                          methodargslist$data_set_parameters$treatment_variable_control_value)
    cv.list <- lapply(1:n.folds,function(i){
      data.train <- data.cv[!idx==i,]
      data.test <- data.cv[idx==i,]
      methodargslist$data_set_parameters$data_set <- data.train
      methodargslist$algorithm_parameters$min_subgroup_size <- floor(nrow(data.train)*0.1)
      methodargslist$algorithm_parameters$random_seed <- seed+i
      res <- SIDES(methodargslist)
      sub <- get_top_subgroup(res,1)
      if (is.null(sub))
        SIDES.sub <- "no subgroup"
      else {
        subs <- lapply(unlist(strsplit(sub,";")),function(x){
          sign <- paste(unlist(strsplit(x,"[[:alnum:],]",perl = T)),collapse = "")
          if (sign != "=")
            return(x)
          else {
            xx <- unlist(strsplit(x,"="))
            return(sprintf("%s %%in%% c(%s)",xx[1],xx[2]))
          }
        })
        SIDES.sub <- paste(sprintf("(%s)",unlist(subs)),collapse = " & ")
      }
      trainobjv <- objF.logrank(subgroups=SIDES.sub, 
                           data=data.train,
                           trt = SIDES.trt,
                           ctrl = SIDES.ctrl,status=methodargslist$data_set_parameters$outcome_censor_name,
                           time=methodargslist$data_set_parameters$outcome_variable_name,
                           cutoff=0,
                           int=FALSE)
      testobjv <- objF.logrank(subgroups=SIDES.sub, 
                           data=data.test,
                           trt = SIDES.trt,
                           ctrl = SIDES.ctrl,status=methodargslist$data_set_parameters$outcome_censor_name,
                           time=methodargslist$data_set_parameters$outcome_variable_name,
                           cutoff=0,
                           int=FALSE)
      testres <- list(bestsub=SIDES.sub,trainobjv=trainobjv,testobjv=testobjv)
      return(list(testres=testres))
    })
  } else if (method == "PRIM") {
    PRIM.trt <- sprintf("%s=%s",methodargslist$trtvar,
                         unique(data.cv[,methodargslist$trtvar])[unique(data.cv[,methodargslist$trtvar]) != methodargslist$trtref])
    PRIM.ctrl <- sprintf("%s=%s",methodargslist$trtvar,methodargslist$trtref)
    cv.list <- lapply(1:n.folds,function(i){
      data.train <- data.cv[!idx==i,]
      data.test <- data.cv[idx==i,]
      methodargslist$data <- data.train
      res <- do.call("prim.train",methodargslist)
      if (is.null(res))
        PRIM.sub <- "no subgroup"
      else PRIM.sub <- paste(sprintf("(%s)",apply(res,1,paste,collapse=" ")),collapse = " & ")
      trainobjv <- objF.logrank(subgroups=PRIM.sub, 
                                data=data.train,
                                trt = PRIM.trt,
                                ctrl = PRIM.ctrl,status=methodargslist$censorvar,
                                time=methodargslist$yvar,
                                cutoff=0,
                                int=FALSE)
      testobjv <- objF.logrank(subgroups=PRIM.sub, 
                               data=data.test,
                               trt = PRIM.trt,
                               ctrl = PRIM.ctrl,status=methodargslist$censorvar,
                               time=methodargslist$yvar,
                               cutoff=0,
                               int=FALSE)
      testres <- list(bestsub=PRIM.sub,trainobjv=trainobjv,testobjv=testobjv)
      return(list(testres=testres))
    })
  }
  
  cv.subs <- unlist(lapply(cv.list,function(x) x$testres$bestsub))
  cv.trainp <- unlist(lapply(cv.list,function(x) x$testres$trainobjv$objv))
  cv.traind <- unlist(lapply(cv.list,function(x) x$testres$trainobjv$diff))
  cv.trainprop <- unlist(lapply(cv.list,function(x) x$testres$trainobjv$sub.proption))
  cv.testp <- unlist(lapply(cv.list,function(x) x$testres$testobjv$objv))
  cv.testd <- unlist(lapply(cv.list,function(x) x$testres$testobjv$diff))
  if (outcome.type %in% c("tte","survival","s") ) cv.testd[cv.testd > 100 | cv.testd < 0.001] <- NA ## reset too large HR value to NA
  cv.testprop <- unlist(lapply(cv.list,function(x) x$testres$testobjv$sub.proption))
  if (method == "SIGA")
    cv.bootres <- lapply(cv.list,function(x) x$bootres)
  if (method %in% c("parSIGA","SIGA.tte","SIGA.binary")) {
    vars <- apply(do.call(rbind,lapply(cv.list,function(x) x$vars)),1,paste,collapse=" ")
    objv <- unlist(lapply(cv.list,function(x) x$objv))
    res <- data.frame(subs=cv.subs,trainp=cv.trainp,testp=cv.testp,traind=cv.traind,testd=cv.testd,trainprop=cv.trainprop,testprop=cv.testprop,vars=vars,objv=objv)
  } else res <- data.frame(subs=cv.subs,trainp=cv.trainp,testp=cv.testp,traind=cv.traind,testd=cv.testd,trainprop=cv.trainprop,testprop=cv.testprop)
  if (method == "SIGA")
    finalres <- list(res=res,cv.bootres=cv.bootres)
  else finalres <- list(res=res)
  return(finalres)
}
