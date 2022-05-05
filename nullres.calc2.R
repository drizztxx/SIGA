

nullres.calc2 <- function(cvres){
  vars <- do.call(rbind,lapply(cvres,function(x){
    aa <- x$cv$subs
    aa.split <- do.call(rbind,lapply(strsplit(aa,"&"),function(x){
      split <- as.data.frame(do.call(rbind,strsplit(trimws(gsub("[()]","",x))," ")))
      colnames(split) <- c("var","sign","value")
      res <- split["var"]
      return(res)
    }))
    return(aa.split)
  }))
  return(table(vars)/NROW(vars))
}

test <- nullres.calc2(nullsimul1)
test2 <- nullres.calc2(nullsimul2)
