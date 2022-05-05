freqvars <- function(bootres,freq=0.5){
  subs <- bootres$subs
  vars <- do.call(rbind,lapply(subs,function(x){
    aa <- x
    if (!is.character(aa))
      aa <- as.character(aa)
    aa.split <- do.call(rbind,lapply(strsplit(aa,"[&|]"),function(x){
      split <- as.data.frame(do.call(rbind,strsplit(trimws(gsub("[()]","",x)),"[ ]+")))
      colnames(split) <- c("var","sign","value")
      res <- split["var"]
      return(unique(res))
    }))
    return(aa.split)
  }))
  varsfreq<-table(vars)[order(table(vars),decreasing = T)]/NROW(subs)
  return(names(varsfreq)[varsfreq>=freq])
}
