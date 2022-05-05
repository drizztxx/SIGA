predict.SIGA <- function(sigaobject,newdata=NULL){
  if (!inherits(sigaobject,"SIGA")) stop ("Un-legal SIGA object")
  if(is.null(sigaobject$bestsub)){
    Best <- Best.SIGA(sigaobject)
    bestsub <- Best$best.subgroup
  } else bestsub <- sigaobject$bestsub
  ori.args <- as.list(sigaobject$call)
  if (!inherits(try(get(sigaobject$dsn),T),"try-error"))
    data <- get(sigaobject$dsn)
  else if (!inherits(try(dynGet(sigaobject$dsn),T),"try-error"))
    data <- dynGet(sigaobject$dsn)
  if (ori.args$resptyp == "tte" || sigaobject$resptyp == "tte"){
    trtvar <- strsplit(ori.args$trt,'=')[[1]][1]
    trainobjv <- do.call("objF.logrank",list(subgroups=bestsub, 
                                    data=data[,c(sigaobject$varused,ori.args$time,ori.args$status,trtvar)],
                                    trt = ori.args$trt,
                                    ctrl = ori.args$ctrl,status=ori.args$status,time=ori.args$time,cutoff=0,
                                    int=FALSE))
    ObjV <- do.call("objF.logrank",list(subgroups=bestsub, 
                                              data=newdata[,c(sigaobject$varused,ori.args$time,ori.args$status,trtvar)],
                                              trt = ori.args$trt,
                                              ctrl = ori.args$ctrl,status=ori.args$status,time=ori.args$time,cutoff=0,
                                              int=FALSE))
    ObjV2 <- do.call("objF.logrank",list(subgroups=sprintf("!(%s)",bestsub), 
                                    data=newdata[,c(sigaobject$varused,ori.args$time,ori.args$status,trtvar)],
                                    trt = ori.args$trt,
                                    ctrl = ori.args$ctrl,status=ori.args$status,time=ori.args$time,cutoff=0,
                                    int=FALSE))
  } else if (ori.args$resptyp == "binary" || sigaobject$resptyp == "binary") {
    trtvar <- strsplit(ori.args$trt,'=')[[1]][1]
    resvar <- strsplit(ori.args$resp,'=')[[1]][1]
    trainobjv <- do.call("objF.glm2",list(subgroups=bestsub, 
                                     data=data[,c(sigaobject$varused,resvar,trtvar)],
                                     ifelse(is.null(ori.args$cutoff),0,ori.args$cutoff),ori.args$trt,ori.args$ctrl,ori.args$resp,resvar,trtvar,FALSE,
                                     family=ifelse(is.null(ori.args$family),"binomial",ori.args$family) ))
    ObjV <- do.call("objF.glm2",list(subgroups=bestsub, 
                                              data=newdata[,c(sigaobject$varused,resvar,trtvar)],
                                              ifelse(is.null(ori.args$cutoff),0,ori.args$cutoff),ori.args$trt,ori.args$ctrl,ori.args$resp,resvar,trtvar,FALSE,
                                              family=ifelse(is.null(ori.args$family),"binomial",ori.args$family) ))
    ObjV2 <- do.call("objF.glm2",list(subgroups=sprintf("!(%s)",bestsub), 
                                     data=newdata[,c(sigaobject$varused,resvar,trtvar)],
                                     ifelse(is.null(ori.args$cutoff),0,ori.args$cutoff),ori.args$trt,ori.args$ctrl,ori.args$resp,resvar,trtvar,FALSE,
                                     family=ifelse(is.null(ori.args$family),"binomial",ori.args$family) ))
  }
  return(list(bestsub=bestsub,trainobjv=trainobjv,testobjv=ObjV,testobjv2=ObjV2))
}