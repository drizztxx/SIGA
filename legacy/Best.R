Best.SIGA <- function(sigaobject){
  if (!inherits(sigaobject,"SIGA"))
    stop("Un-legal SIGA object")
  minobjv <- min(sigaobject$Best)
  if (minobjv == 1)
    stop("minimum objective value equals to 1")
  minind <- unlist(sigaobject$Objv) == minobjv
  all.subs <- unique(unlist(sigaobject$Subset)[minind])
  ori.args <- as.list(sigaobject$call)
  if (!inherits(try(get(sigaobject$dsn),T),"try-error"))
    data <- get(sigaobject$dsn)
  else if (!inherits(try(dynGet(sigaobject$dsn),T),"try-error"))
    data <- dynGet(sigaobject$dsn)
  else stop(sprintf("%s cannot be get",sigaobject$dsn))
  if (ori.args$resptyp == "tte" || sigaobject$resptyp == "tte"){
    trtvar <- strsplit(ori.args$trt,'=')[[1]][1]
    ObjV <- do.call(ori.args$objFunction,list(subgroups=all.subs, 
                                     data=data[,c(sigaobject$varused,ori.args$time,ori.args$status,trtvar)],
                                     trt = ori.args$trt,
                                     ctrl = ori.args$ctrl,status=ori.args$status,time=ori.args$time,cutoff=ifelse(is.null(ori.args$cutoff),10,ori.args$cutoff) ))
  } else {
    trtvar <- strsplit(ori.args$trt,'=')[[1]][1]
    resvar <- strsplit(ori.args$resp,'=')[[1]][1]
    ObjV <- do.call(ori.args$objFunction,list(subgroups=all.subs, 
                                     data=data[,c(sigaobject$varused,resvar,trtvar)],
                                     ifelse(is.null(ori.args$cutoff),10,ori.args$cutoff),ori.args$trt,ori.args$ctrl,ori.args$resp,resvar,trtvar,
                                     family=ifelse(is.null(ori.args$family),"binomial",ori.args$family) ))
  }
  all.matrix <- cbind(all.subs,do.call(cbind,ObjV))
  all.matrix <- all.matrix[all.matrix[,4] == max(all.matrix[,4]),]
  all.subs.split <- do.call(rbind,lapply(strsplit(all.subs,"&"),function(x){
    split <- do.call(rbind,strsplit(trimws(x)," "))
    colnames(split) <- c("var","sign","value")
    nn <- nrow(split)
    return(data.frame(split,nn))
  }))
  all.subs.split <- all.subs.split[all.subs.split[,4] == min(all.subs.split[,4]),]
  all.subs.split <- all.subs.split[order(all.subs.split[,1]),]
  all.subs.split <- unique(all.subs.split)
  best.subdata <- subset(data, eval(parse(text = all.subs[1])))
  tocomb <- unlist(lapply(unique(as.character(all.subs.split$var)),function(x){
    ##x<-unique(all.subs.split$var)[5]
    subs <- all.subs.split[all.subs.split$var == x,]
    if (NROW(unique(subs$sign)) > 1){
      subs.min <- subs[subs$sign == ">=",]
      subs.max <- subs[subs$sign == "<=",]
      min.value <- min(as.numeric(as.character(best.subdata[,x])),na.rm = T)
      max.value <- max(as.numeric(as.character(best.subdata[,x])),na.rm = T)
      sub.desc <- sprintf("%s >= %s | %s <= %s",x,min.value,x,max.value)
      return(sub.desc)
    } else if (NROW(unique(subs$sign)) == 1) {
      if (unique(subs$sign) == ">=") {
        subs.min <- subs[subs$sign == ">=",]
        min.value <- min(as.numeric(as.character(best.subdata[,x])),na.rm = T)
        sub.desc <- sprintf("%s >= %s",x,min.value)
      } else if (unique(subs$sign) == "<=") {
        subs.max <- subs[subs$sign == "<=",]
        max.value <- max(as.numeric(as.character(best.subdata[,x])),na.rm = T)
        sub.desc <- sprintf("%s <= %s",x,max.value)
      } else if (unique(subs$sign) == "=="){
        if (NROW(unique(best.subdata[,x])) == 1) {
          sub.desc <- sprintf("%s == '%s'",x,as.character(unique(best.subdata[,x])))
        } else sub.desc <- NULL
      } else if (unique(subs$sign) == "%in%") {
        combb <- paste("\'",as.character(unique(best.subdata[,x])), "\'",sep="",collapse = ",")
        cutoff <- paste0("c(",combb,")")
        sub.desc <- paste(x,cutoff,sep=" %in% ",collapse = " & ")
      }
      return(sub.desc)
    }
  }))
  best.subgroup <- paste(sprintf("(%s)",tocomb),collapse = " & ")
  best.subdata2 <- subset(data, eval(parse(text = best.subgroup)))
  if (!identical(best.subdata,best.subdata2)) {
    warning("Datasets are not identical, see all.matrix for more details")
    return(list(best.subgroup=best.subgroup,ObjV=list(objv=ObjV$objv[1],diff=ObjV$diff[1],sub.size=ObjV$sub.size[1],sub.proption=ObjV$sub.proption[1]),all.matrix=all.matrix))
  } else { ## Try to shrinkage the nubmer of variables
    shrinkage <- do.call(rbind,lapply(1:NROW(tocomb),function(i){
      combned <- combn(tocomb,i)
      do.call(rbind,lapply(1:NCOL(combned),function(xx){
        x <- combned[,xx]
        best.cand <- paste(sprintf("(%s)",x),collapse = " & ")
        best.subdata3 <- subset(data, eval(parse(text = best.cand)))
        if (identical(best.subdata,best.subdata3)) {
          selected <- 1
        } else selected <- 0
        return(cbind(best.cands=best.cand,selected=selected))
      }))
    }))
    orivarnum <- as.numeric(as.character(all.subs.split[1,"nn"]))
    newvarnum <- NROW(unlist(strsplit(trimws(shrinkage[which.max(shrinkage[,2]),1]),"&")))
    if (newvarnum < orivarnum ) {
      cat(sprintf("The number of variables has been shrinked from %s to %s",orivarnum,newvarnum))
      return(list(best.subgroup=shrinkage[which.max(shrinkage[,2]),1],
                  ObjV=list(objv=ObjV$objv[1],diff=ObjV$diff[1],sub.size=ObjV$sub.size[1],sub.proption=ObjV$sub.proption[1]),
                  all.matrix=all.matrix))
    } else {
      return(list(best.subgroup=best.subgroup,
                  ObjV=list(objv=ObjV$objv[1],diff=ObjV$diff[1],sub.size=ObjV$sub.size[1],sub.proption=ObjV$sub.proption[1]),
                  all.matrix=all.matrix))
    }
  }
}
