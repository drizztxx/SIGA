## Decoding the generated Chrom according to the genetic representation rule.
## The genetic representation rule is as below:
## 1. Value of gene(B) is ranged from 0-4 (uniform distribution).
## 2. For continuous predictor, if B <= 2, it will not be selected (50% possibility);
##    if 2 < B <=3, then the subgroup will be defined as <= (B - 2) quantile;
##    if 3 < B <=4, then the subgroup will be defined as >= (B - 3) quantile.
## 3. For binary predictor, B <= 2 indicates the variable will not be selected;
##    and 2 < B <= 3 indicates the first category(level of factor) will be used to define the subgroup 
##    whilst 3 < B <= 4 indicates the second one will be used.
## 4. For multi-categorical predictor (i=1 to C), it will not be selected if B <= 2;
##    else there will be (2^C - 2) subgroups to be choosed and each subgroup is one possible
##    combination of the category for the predictor. The value of B (2-4) will be equally divided to
##    (2^C - 2) blocks and each block stands for one category combination.
decodsubset <- function(data,Chrom){
  if (!is.data.frame(data))
    data <- as.data.frame(data)
  if (ncol(data) != ncol(Chrom)) 
    stop("Error: number of columns not identical")
  columns <- colnames(data)
  var.class <- sapply(data,class)
  var.numeric <- names(var.class)[var.class %in% c("numeric","integer")]
  if(!anyNA(var.numeric)) data.numeric <- data[,var.numeric]
  var.numeric.idx <- columns %in% var.numeric
  var.factor <- names(var.class)[var.class == "factor"]
  if (NROW(var.factor) == 0) {
    levels <- nlevels(data[,var.factor])
  } else if (NROW(var.factor) == 1) {
    levels <- nlevels(data[,var.factor])
  } else if (NROW(var.factor) > 1) {
    levels <- sapply(data[,var.factor],nlevels)}
  var.binary <- var.factor[levels == 2]
  if(!anyNA(var.binary)) data.binary <- data[,var.binary]
  var.binary.idx <- columns %in% var.binary
  var.mc <- var.factor[levels > 2]
  if(!anyNA(var.mc)) data.mc <- data[,var.mc]
  var.mc.idx <- columns %in% var.mc
  
  subset.desc <- vector()
  decod.list <- lapply(1:NROW(Chrom),function(x){
    if (!anyNA(var.numeric) & length(var.numeric) > 0) { ## Decoding numeric variables
      Chrom.sub <- Chrom[x,var.numeric.idx]
      vars.plus <- var.numeric[ifelse(Chrom.sub > 3, T, F)]
      vars.minus <- var.numeric[ifelse(Chrom.sub > 2 & Chrom.sub <= 3, T, F)]
      if (NROW(vars.plus) >= 2) {
        quantiles <- Chrom.sub[ifelse(Chrom.sub > 3, T, F)] - 3
        cutoff <- diag(apply(data.numeric[,vars.plus],2,quantile,prob=quantiles,na.rm=T))
        condition1 <- paste(vars.plus,cutoff,sep=" >= ",collapse = " & ")
      } else if(NROW(vars.plus) == 1) {
        quantiles <- Chrom.sub[ifelse(Chrom.sub > 3, T, F)] - 3
        cutoff <- quantile(data.numeric,prob=quantiles,na.rm=T)
        condition1 <- paste(vars.plus,cutoff,sep=" >= ",collapse = " & ")
      } else {
        condition1 <- NULL
      }
      if (NROW(vars.minus) >= 2) {
        quantiles <- Chrom.sub[ifelse(Chrom.sub > 2 & Chrom.sub <= 3, T, F)] - 2
        cutoff <- diag(apply(data.numeric[,vars.minus],2,quantile,prob=quantiles,na.rm=T))
        condition2 <- paste(vars.minus,cutoff,sep=" <= ",collapse = " & ")
      } else if(NROW(vars.minus) == 1) {
        quantiles <- Chrom.sub[ifelse(Chrom.sub > 2 & Chrom.sub <= 3, T, F)] - 2
        cutoff <- quantile(data.numeric,prob=quantiles,na.rm=T)
        condition2 <- paste(vars.minus,cutoff,sep=" <= ",collapse = " & ")
      } else {
        condition2 <- NULL
      }
    } else {condition1 <- NULL; condition2 <- NULL}
    if (!anyNA(var.binary) & length(var.binary) > 0) { ## Decoding binary variables
      Chrom.sub <- Chrom[x,var.binary.idx]
      vars.plus <- var.binary[ifelse(Chrom.sub > 3, T, F)]
      vars.minus <- var.binary[ifelse(Chrom.sub > 2 & Chrom.sub <= 3, T, F)]
      if (NROW(vars.plus) >= 2) {
        level <- sapply(data.binary[,vars.plus],levels)
        cutoff <- paste("\'",level[2,], "\'",sep="")
        condition3 <- paste(vars.plus,cutoff,sep=" == ",collapse = " & ")
      } else if(NROW(vars.plus) == 1) {
        level <- levels(data.binary[vars.plus])
        cutoff <- paste("\'",level[2], "\'",sep="")
        condition3 <- paste(vars.plus,cutoff,sep=" == ",collapse = " & ")
      } else {
        condition3 <- NULL
      }
      if (NROW(vars.minus) >= 2) {
        level <- sapply(data.binary[,vars.minus],levels)
        cutoff <- paste("\'",level[1,], "\'",sep="")
        condition4 <- paste(vars.minus,cutoff,sep=" == ",collapse = " & ")
      } else if(NROW(vars.minus) == 1) {
        level <- levels(data.binary[vars.minus])
        cutoff <- paste("\'",level[1], "\'",sep="")
        condition4 <- paste(vars.minus,cutoff,sep=" == ",collapse = " & ")
      } else {
        condition4 <- NULL
      }
    } else {condition3 <- NULL; condition4 <- NULL}
    if (!anyNA(var.mc) & length(var.mc) > 0) { ## Decoding multi-category variables
      Chrom.sub <- Chrom[x,var.mc.idx]
      vars <- var.mc[Chrom.sub > 2]
      chrom <- Chrom.sub[Chrom.sub > 2]
      if (NROW(vars) >= 2) {
        cutoff <- sapply(1:NROW(vars),function(xx){
          nlevel <- nlevels(data.mc[,vars[xx]])
          level <- levels(data.mc[,vars[xx]])
          chrom.sub <- chrom[xx]
          ncomb <- 2^nlevel - 2
          comb.list <- unlist(sapply(1:(nlevel - 1), function(y) combn(level,y,simplify = F)),recursive = F)
          idx <- ceiling((chrom.sub - 2)*(ncomb/2)) 
          comb <- comb.list[[idx]]
          combb <- paste("\'",comb, "\'",sep="",collapse = ",")
          res <- paste0("c(",combb,")")
          return(res)
        })
        condition5 <- paste(vars,cutoff,sep=" %in% ",collapse = " & ")
      } else if(NROW(vars) == 1) {
        nlevel <- nlevels(data.mc[vars])
        level <- levels(data.mc[vars])
        chrom.sub <- chrom
        ncomb <- 2^nlevel - 2
        comb.list <- unlist(sapply(1:(nlevel - 1), function(y) combn(level,y,simplify = F)),recursive = F)
        idx <- ceiling((chrom.sub - 2)*(ncomb/2)) 
        comb <- comb.list[[idx]]
        combb <- paste("\'",comb, "\'",sep="",collapse = ",")
        cutoff <- paste0("c(",combb,")")
        condition5 <- paste(vars,cutoff,sep=" %in% ",collapse = " & ")
      } else {
        condition5 <- NULL
      }
    } else {condition5 <- NULL}
    subset.condition <- paste(c(condition1,condition2,condition3,condition4,condition5),collapse = " & ")
    subset.desc[x] <<- subset.condition
    if (subset.condition != "") {
      subgroup <- subset(data, eval(parse(text = get("subset.condition"))))
    } else {
      subset.desc[x] <<- "no subgroup"
      subgroup <- data
    }
    return(subgroup)
  })
  return(list(subgroup.list=decod.list,subgroup.desc=subset.desc))
}
