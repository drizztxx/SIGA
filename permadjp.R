permadjp <- function(sigaobject,nperm=1000,subgroup=NULL,seed=123){
  if (!inherits(sigaobject,"SIGA")) stop ("Un-legal SIGA object")
  if (is.null(subgroup)){
    Best <- Best.SIGA(sigaobject)
    bestsub <- Best$best.subgroup
    trainobjv <- Best$ObjV
  } else bestsub <- subgroup
  ori.args <- as.list(sigaobject$call)
  if (!inherits(try(get(sigaobject$dsn),T),"try-error"))
    data <- get(sigaobject$dsn)
  else if (!inherits(try(dynGet(sigaobject$dsn),T),"try-error"))
    data <- dynGet(sigaobject$dsn)
  else stop(sprintf("%s cannot be get",sigaobject$dsn))
  if (!is.null(sigaobject$resptyp))
    resptyp <- sigaobject$resptyp
  else if (!is.null(ori.args$resptyp))
    resptyp <- ori.args$resptyp
  if (resptyp == "tte") {
    trtvar <- strsplit(ori.args$trt,'=')[[1]][1]
    var.fix <- c(trtvar,ori.args$time,ori.args$status) 
    var.perm <- sigaobject$varused
  } else if (resptyp == "binary") {
    trtvar <- strsplit(ori.args$trt,'=')[[1]][1]
    resvar <- strsplit(ori.args$resp,'=')[[1]][1]
    var.fix <- c(trtvar,resvar) 
    var.perm <- sigaobject$varused
  }
  permdata.list <- lapply(1:nperm,function(x){
    set.seed(seed+x)
    perm.row <- sample(1:nrow(data))
    perm <- data[perm.row,var.perm]
    df <- data.frame(data[,var.fix],perm)
    return(df)
  })
  permdata.list[[nperm+1]] <- data
  if (resptyp == "tte"){
    permobjv.list <- lapply(permdata.list,function(df){
      ObjV <- do.call("objF.cox",list(subgroups=bestsub, 
                                      data=df,
                                      trt = ori.args$trt,
                                      ctrl = ori.args$ctrl,status=ori.args$status,time=ori.args$time,cutoff=0,
                                      int=FALSE))
      return(ObjV)
      })
  } else if (resptyp == "binary") {
    permobjv.list <- lapply(permdata.list,function(df){
      ObjV <- do.call("objF.glm",list(subgroups=bestsub, 
                                       data=df,
                                       0,ori.args$trt,ori.args$ctrl,ori.args$resp,resvar,trtvar,
                                       family=ifelse(is.null(ori.args$family),"binomial",ori.args$family)))
      return(ObjV)
      })
  }
  perm.pvalue <- unlist(lapply(permobjv.list,function(x) x$objv))
  train.rank <- rank(perm.pvalue)[nperm+1]
  adjpvalue <- train.rank/(nperm+1)
  return(list(adjp=adjpvalue,all.pvalues=perm.pvalue))
}