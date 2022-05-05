
VariableScreen <- function(data,
                           trt="treat=='E'",
                           ctrl="treat=='C'",
                           resp="res==1",
                           status="status",
                           time="time",
                           varScreen=c(1,5),
                           resptyp=c("binary","tte"),
                           ...){
  resptyp = match.arg(resptyp)
  if (resptyp == "binary"){
    trtvar <- strsplit(trt,'=')[[1]][1]
    resvar <- strsplit(resp,'=')[[1]][1]
    newdata <- data
    newdata[,trtvar] <- factor(newdata[,trtvar],labels = c(1,-1))
    newdata[,resvar] <- as.factor(newdata[,resvar])
    select1 <- colnames(newdata)[!colnames(newdata) %in% trtvar]
    select2 <- colnames(newdata)[!colnames(newdata) %in% c(trtvar,resvar)]
    if (varScreen[1] == 1){ #trigger RF method
      library(randomForest)
      rf1 <- randomForest(as.formula(paste(resvar,'.',sep='~')),data=subset(newdata[newdata[,trtvar]==1,],select = select1),
                          ntree = 2000,importance=T,na.action = na.omit)
      vim1 <- importance(rf1,type=1,scale=F)
      rf2 <- randomForest(as.formula(paste(resvar,'.',sep='~')),data=subset(newdata[newdata[,trtvar]==-1,],select = select1),
                          ntree = 2000,importance=T,na.action = na.omit)
      vim2 <- importance(rf2,type=1,scale=F)
      vim <- vim1 - vim2
      col.vim <- rownames(vim)
      col.vim <- col.vim[order(abs(vim),decreasing = T)]
    } else if (varScreen[1] == 2){#trigger P valule method, for binary endpoint, logistic regression will be fitted.
      p.list <- lapply(select2,function(x){
        frmla <- as.formula(sprintf("%s ~ %s + %s", resvar, trtvar,x))
        model1 <- glm(frmla,family = binomial("logit"), data = newdata)
        frmla2 <- as.formula(sprintf("%s ~ %s * %s", resvar, trtvar,x))
        model2 <- glm(frmla2,family = binomial("logit"), data = newdata)
        test <- anova(model1,model2)
        pval <- test[2,4]
        if (is.na(pval)) pval <- 1
        return(pval)
      })
      vim <- unlist(p.list)
      names(vim) <- select2
      col.vim <- select2[order(vim)]
    }
  } else if  (resptyp == "tte"){
    trtvar <- strsplit(trt,'=')[[1]][1]
    newdata <- data
    newdata[,trtvar] <- factor(newdata[,trtvar],labels = c(1,-1))
    select1 <- colnames(newdata)[!colnames(newdata) %in% trtvar]
    select2 <- colnames(newdata)[!colnames(newdata) %in% c(trtvar,time,status)]
    if (varScreen[1] == 1){ #trigger RFsrc method
      library(randomForestSRC)
      options(rf.cores=1)
      rf1 <- rfsrc(as.formula(sprintf("Surv(%s, %s) ~ .",time,status)),data=subset(newdata[newdata[,trtvar]==1,],select = select1)
                   ,ntree = 2000,importance = T,na.action = "na.omit",block.size = 1)
      vim1 <- rf1$importance
      rf2 <- rfsrc(as.formula(sprintf("Surv(%s, %s) ~ .",time,status)),data=subset(newdata[newdata[,trtvar]==-1,],select = select1)
                   ,ntree = 2000,importance = T,na.action = "na.omit",block.size = 1)
      vim2 <- rf2$importance
      vim <- vim1 - vim2
      col.vim <- names(vim)
      col.vim <- col.vim[order(abs(vim),decreasing = T)]
    } else if (varScreen[1] == 2){#trigger P valule method, for survival endpoint, Cox regression will be fitted.
      p.list <- lapply(select2,function(x){
        frmla <- as.formula(sprintf("Surv(%s, %s) ~ %s + %s ", time,status,trtvar, x))
        try1 <- try(model1 <- coxph(frmla,data = newdata),silent = T)
        frmla2 <- as.formula(sprintf("Surv(%s, %s) ~ %s * %s ", time,status,trtvar, x))
        try2 <- try(model2 <- coxph(frmla2,data = newdata),silent = T)
        if (any(inherits(try1,"try-error"),inherits(try2,"try-error")))
          pval <- 1
        else {
          test <- anova(model1,model2)
          pval <- test[2,4]
          if (is.na(pval)) pval <- 1
        }
        return(pval)
      })
      vim <- unlist(p.list)
      names(vim) <- select2
      col.vim <- select2[order(vim)]
    }
  }
  return(list(col.vim=col.vim,vim=vim))
}
