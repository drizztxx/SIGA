nullres.calc <- function(cvres,cutoff=0.1){
  cv.matrix <- do.call(rbind,lapply(cvres,function(x){
    cv <- x$cv
    hr <- x$hr
    meanhr <- mean(cv$testd,na.rm = T)
    if (meanhr >= 1 || meanhr > hr*(1-cutoff))
      res <- 0
    else res <- 1
    trainprop <- mean(as.numeric(gsub("[%]","",cv$trainprop)),na.rm = T)
    testprop <- mean(as.numeric(gsub("[%]","",cv$testprop)),na.rm = T)
    return(data.frame(res=res,hr=hr,meanhr=meanhr,trainprop=trainprop,testprop=testprop))
  }))
  rate <- sum(cv.matrix$res)/length(cvres)
  return(list(success.rate=rate,cv.matrix=cv.matrix))
}
