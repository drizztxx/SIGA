## Caculate the treatment effect (that is the odds ratio between the two groups) for every subgroups
objF.OR <- function(subgroups,data,cutoff=10,trt,ctrl,resp){
  npop <- NROW(data)
  objv <- lapply(subgroups,function(xx){
    if (!grepl("no subgroup",xx)) {
      sub <- subset(data, eval(parse(text = get("xx"))))
    } else sub <- data
    sub.treat <- subset(sub,eval(parse(text = get("trt"))))
    resp.treat <- NROW(subset(sub.treat,eval(parse(text = get("resp")))))
    if (resp.treat == NROW(sub.treat)) {
      odds.treat <- resp.treat
    } else odds.treat <- resp.treat/(NROW(sub.treat) - resp.treat)
    sub.control <- subset(sub,eval(parse(text = get("ctrl"))))
    resp.control <- NROW(subset(sub.control,eval(parse(text = get("resp")))))
    if (resp.control == NROW(sub.control)) {
      odds.control <- resp.control
    } else odds.control <- resp.control/(NROW(sub.control) - resp.control)
    or <- odds.treat/odds.control
    if (is.infinite(or) || is.na(or)) or <- 0
    nsub <- NROW(sub)
    prop <- nsub*100/npop
    return(c(or=-or,sub.size=nsub,sub.proption=prop))
  })
  objv2 <- do.call(rbind,objv)
  objv3 <- apply(objv2,1,function(vc){
    discount <- vc[3]/cutoff
    if (discount == 0) {
      objv <- vc[1]
    } else if (vc[3] >= cutoff) {
      objv <- vc[1]
    } else if (vc[1] >= 0) {
      objv <- vc[1]/discount
    } else if (vc[1] < 0) {
      objv <- vc[1]*discount
    }
    return(objv)
  })
  return(list(objv=objv3,or=objv2[,1],sub.size=objv2[,2],sub.proption=sprintf("%s%%",objv2[,3])))
}

## Caculate the treatment effect (that is the difference of response rate between the two groups) for every subgroups
objF.difference2 <- function(subgroups,data,cutoff=10,trt,ctrl,resp){
  npop <- NROW(data)
  sub.ctrl <- subset(data, eval(parse(text = get("ctrl"))))
  sub.ctrl.res <- subset(sub.ctrl, eval(parse(text = get("resp"))))
  rr.target <- NROW(sub.ctrl.res)/NROW(sub.ctrl)
  objv <- lapply(subgroups,function(xx){
    if (!grepl("no subgroup",xx)) {
      sub <- subset(data, eval(parse(text = get("xx"))))
    } else sub <- data
    sub.treat <- subset(sub,eval(parse(text = get("trt"))))
    rr.treat <- NROW(subset(sub.treat,eval(parse(text = get("resp")))))/NROW(sub.treat)
    sub.control <- subset(sub,eval(parse(text = get("ctrl"))))
    rr.control <- NROW(subset(sub.control,eval(parse(text = get("resp")))))/NROW(sub.control)
    rr.difference <- rr.treat - rr.control
    if (is.nan(rr.difference)) rr.difference <- -1
    nsub <- 2*min(NROW(sub.treat),NROW(sub.control))
    prop <- nsub*100/npop
    return(c(diff=-rr.difference,sub.size=nsub,sub.proption=prop,rr.control=rr.control))
  })
  objv2 <- do.call(rbind,objv)
  objv3 <- apply(objv2,1,function(vc){
    newobjv <- ifelse(vc[4] >= rr.target,vc[1],(rr.target + vc[1] - vc[4]))
    newobjv <- ifelse(is.na(newobjv),1,newobjv)
    discount <- vc[3]/cutoff
    if (discount == 0 | vc[3] >= cutoff) {
      objv <- newobjv
    } else if (newobjv >= 0) {
      objv <- newobjv/discount
    } else if (newobjv < 0) {
      objv <- newobjv*discount
    }
    return(objv)
  })
  return(list(objv=objv3,diff=objv2[,1],sub.size=objv2[,2],sub.proption=sprintf("%s%%",objv2[,3])))
}

## Caculate the treatment effect (that is the difference of response rate between the two groups) for every subgroups
objF.difference <- function(subgroups,data,cutoff=10,trt,ctrl,resp){
  npop <- NROW(data)
  sub.ctrl <- subset(data, eval(parse(text = get("ctrl"))))
  sub.ctrl.res <- subset(sub.ctrl, eval(parse(text = get("resp"))))
  rr.target <- NROW(sub.ctrl.res)/NROW(sub.ctrl)
  objv <- lapply(subgroups,function(xx){
    if (!grepl("no subgroup",xx)) {
      sub <- subset(data, eval(parse(text = get("xx"))))
    } else sub <- data
    sub.treat <- subset(sub,eval(parse(text = get("trt"))))
    rr.treat <- NROW(subset(sub.treat,eval(parse(text = get("resp")))))/NROW(sub.treat)
    sub.control <- subset(sub,eval(parse(text = get("ctrl"))))
    rr.control <- NROW(subset(sub.control,eval(parse(text = get("resp")))))/NROW(sub.control)
    rr.difference <- rr.treat - rr.control
    if (is.nan(rr.difference)) rr.difference <- -1
    nsub <- 2*min(NROW(sub.treat),NROW(sub.control))
    prop <- nsub*100/npop
    return(c(diff=-rr.difference,sub.size=nsub,sub.proption=prop,rr.control=rr.control))
  })
  objv2 <- do.call(rbind,objv)
  objv3 <- apply(objv2,1,function(vc){
    newobjv <- vc[1]
    discount <- vc[3]/cutoff
    if (discount == 0 | vc[3] >= cutoff) {
      objv <- newobjv
    } else if (newobjv >= 0) {
      objv <- newobjv/discount
    } else if (newobjv < 0) {
      objv <- newobjv*discount
    }
    return(objv)
  })
  return(list(objv=objv3,diff=objv2[,1],sub.size=objv2[,2],sub.proption=sprintf("%s%%",objv2[,3])))
}

## Caculate the treatment effect (that is the pvalue caculated through the glm model) for every subgroups
objF.glm <- function(subgroups,data,cutoff=10,trt,ctrl,resp,resvar,trtvar,...){
  npop <- NROW(data)
  objv <- lapply(subgroups,function(xx){
    if (!grepl("no subgroup",xx)) {
      sub <- subset(data, eval(parse(text = get("xx"))))
    } else sub <- data
    sub.treat <- subset(sub,eval(parse(text = get("trt"))))
    rr.treat <- NROW(subset(sub.treat,eval(parse(text = get("resp")))))/NROW(sub.treat)
    sub.control <- subset(sub,eval(parse(text = get("ctrl"))))
    rr.control <- NROW(subset(sub.control,eval(parse(text = get("resp")))))/NROW(sub.control)
    rr.difference <- rr.treat - rr.control
    if (is.nan(rr.difference)) rr.difference <- -1
    nsub <- 2*min(NROW(sub.treat),NROW(sub.control))
    prop <- nsub*100/npop
    if (nsub == 0 | prop < cutoff){
      pval <- 1
    } else {
      sum <- summary(glm(formula = as.formula(paste(resvar,trtvar,sep='~')), data = sub, ...))
      zval <- sum$coefficients[2,3]
      pval <- pnorm(zval, lower.tail = F)
    }
    return(c(pval=pval,diff=rr.difference,sub.size=nsub,sub.proption=prop,rr.control=rr.control))
  })
  objv2 <- do.call(rbind,objv)
  return(list(objv=objv2[,1],diff=objv2[,2],sub.size=objv2[,3],sub.proption=sprintf("%s%%",objv2[,4])))
}

## Caculate the treatment effect (that is the pvalue caculated through the glm model) for every subgroups
## Enhanced to make sure there will be an interactive effect: if the treatment effect in the complementory
## set of choosed subgroup less than that in the choosed subgroup, then the p-value will set to 1
objF.glm2 <- function(subgroups,data,cutoff=10,trt,ctrl,resp,resvar,trtvar,int=T,...){
  npop <- NROW(data)
  objv <- lapply(subgroups,function(xx){
    if (!grepl("no subgroup",xx)) {
      sub <- subset(data, eval(parse(text = get("xx"))))
      sub2 <- subset(data, !eval(parse(text = get("xx")))) #the complementory set
    } else {sub <- data; sub2 <- data}
    sub.treat <- subset(sub,eval(parse(text = get("trt"))))
    sub.treat2 <- subset(sub2,eval(parse(text = get("trt"))))
    rr.treat <- NROW(subset(sub.treat,eval(parse(text = get("resp")))))/NROW(sub.treat)
    rr.treat2 <- NROW(subset(sub.treat2,eval(parse(text = get("resp")))))/NROW(sub.treat2)
    sub.control <- subset(sub,eval(parse(text = get("ctrl"))))
    sub.control2 <- subset(sub2,eval(parse(text = get("ctrl"))))
    rr.control <- NROW(subset(sub.control,eval(parse(text = get("resp")))))/NROW(sub.control)
    rr.control2 <- NROW(subset(sub.control2,eval(parse(text = get("resp")))))/NROW(sub.control2)
    rr.difference <- rr.treat - rr.control
    rr.difference2 <- rr.treat2 - rr.control2
    if (is.nan(rr.difference)) rr.difference <- -1
    nsub <- 2*min(NROW(sub.treat),NROW(sub.control))
    nsub2 <- 2*min(NROW(sub.treat2),NROW(sub.control2))
    prop <- nsub*100/npop
    if (int & (nsub == 0 | prop < cutoff | nsub2 == 0)){
      pval <- 1
    } else {
      sum <- summary(glm(formula = as.formula(paste(resvar,trtvar,sep='~')), data = sub, ...))
      sum2 <- summary(glm(formula = as.formula(paste(resvar,trtvar,sep='~')), data = sub2, ...))
      zval <- sum$coefficients[2,3]
      zval2 <- sum2$coefficients[2,3]
      pval <- pnorm(zval, lower.tail = F)
      pval2 <- pnorm(zval2, lower.tail = F)
      if (int & pval2 < pval) pval <- 1
    }
    return(c(pval=pval,diff=rr.difference,sub.size=nsub,sub.proption=prop,sub.treat=NROW(sub.treat),
             sub.control=NROW(sub.control)))
  })
  objv2 <- do.call(rbind,objv)
  return(list(objv=objv2[,1],diff=objv2[,2],sub.size=objv2[,3],sub.proption=sprintf("%s%%",objv2[,4]),sub.treat=objv2[,5],
              sub.control=objv2[,6]))
}

## Caculate the treatment effect (that is the HR and pvalue caculated through the cox model) for time to event endpoints
## Enhanced to make sure there will be an interactive effect: if the treatment effect in the complementory
## set of choosed subgroup less than that in the choosed subgroup, then the p-value will set to 1
objF.cox <- function(subgroups,data,cutoff=10,trt="treat==1",ctrl="treat==2",status="e.rfs",time="t.rfs",int=T){
  require(survival)
  trtvar <- strsplit(trt,"=")[[1]][1]
  trtval <- strsplit(trt,"=")[[1]][2]
  ctrlval <- strsplit(ctrl,"=")[[1]][2]
  data$newtreat <- factor(data[,trtvar],levels = c(trtval, ctrlval))
  npop <- NROW(data)
  objv <- lapply(subgroups,function(xx){
    if (!grepl("no subgroup",xx)) {
      sub <- subset(data, eval(parse(text = get("xx"))))
      sub2 <- subset(data, !eval(parse(text = get("xx")))) #the complementory set
    } else {sub <- data; sub2 <- data}
    sub.treat <- subset(sub,eval(parse(text = sprintf("%s == %s",trtvar,trtval))))
    sub.control <- subset(sub,eval(parse(text = sprintf("%s == %s",trtvar,ctrlval))))
    nsub <- 2*min(NROW(sub.treat),NROW(sub.control))
    prop <- nsub*100/npop
    if (int & (nsub == 0 | prop < cutoff)){
      pval <- 1
      hr <- 1
    } else if (NROW(sub) == 0 | NROW(sub2) == 0) {
      pval <- 1
      hr <- 1
    } else {
      try1 <- try(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub))
      try2 <- try(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub2))
      try3 <- try(summary(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub)),T)
      try4 <- try(summary(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub2)),T)
      if (!int & (inherits(try1,"try-error") | inherits(try3,"try-error") )){
        pval <- 1;hr <- 1
      } else if (int & any(inherits(try1,"try-error"),inherits(try2,"try-error"),inherits(try3,"try-error"),inherits(try4,"try-error"))){
        pval <- 1;hr <- 1
      } else {
        coxfit1 <- coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub)
        coxfit2 <- coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub2)
        sum <- summary(coxfit1)
        sum2 <- summary(coxfit2)
        zval <- sum$coefficients[4]
        zval2 <- sum2$coefficients[4]
        if (!is.na(zval)) pval <- pnorm(zval, lower.tail = F) else pval <- 1
        if (!is.na(zval2)) pval2 <- pnorm(zval2, lower.tail = F) else pval2 <- 1
        hr <- sum$conf.int[2]
        if (int & pval2 < pval) {pval <- 1;hr<-1} 
      }
    }
    return(c(pval=pval,diff=hr,sub.size=nsub,sub.proption=prop,sub.treat=NROW(sub.treat),
             sub.control=NROW(sub.control)))
  })
  objv2 <- do.call(rbind,objv)
  return(list(objv=objv2[,1],diff=objv2[,2],sub.size=objv2[,3],sub.proption=sprintf("%s%%",objv2[,4]),sub.treat=objv2[,5],
              sub.control=objv2[,6]))
}

## Caculate the objective value and estimated treatment effect (that is the HR caculated through the cox model)
## of subgroup in training data
## Objective value is calculated via: 2(1-Phi(|Z+ - Z-|/2^0.5)), where Phi is the CDF of standard normal distribution and
## Z+, Z- are the test statistics for testing the one-sided hypothesis of treatment efficacy in subgroup positive and
## negative population respectively.
objF.cox2 <- function(subgroups,data,cutoff=10,trt="treat==1",ctrl="treat==2",status="e.rfs",time="t.rfs"){
  require(survival)
  trtvar <- strsplit(trt,"=")[[1]][1]
  trtval <- strsplit(trt,"=")[[1]][2]
  ctrlval <- strsplit(ctrl,"=")[[1]][2]
  data$newtreat <- factor(data[,trtvar],levels = c(trtval, ctrlval))
  npop <- NROW(data)
  objv <- lapply(subgroups,function(xx){
    if (!grepl("no subgroup",xx)) {
      sub <- subset(data, eval(parse(text = get("xx"))))
      sub2 <- subset(data, !eval(parse(text = get("xx")))) #the complementory set
    } else {sub <- data; sub2 <- data}
    sub.treat <- subset(sub,eval(parse(text = sprintf("%s == %s",trtvar,trtval))))
    sub.control <- subset(sub,eval(parse(text = sprintf("%s == %s",trtvar,ctrlval))))
    nsub <- 2*min(NROW(sub.treat),NROW(sub.control))
    prop <- nsub*100/npop
    if (nsub == 0 | prop < cutoff){
      pval <- 1
      hr <- 1
    } else if (NROW(sub) == 0 | NROW(sub2) == 0) {
      pval <- 1
      hr <- 1
    } else {
      try1 <- try(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub))
      try2 <- try(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub2))
      try3 <- try(summary(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub)),T)
      try4 <- try(summary(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub2)),T)
      if (any(inherits(try1,"try-error"),inherits(try2,"try-error"),inherits(try3,"try-error"),inherits(try4,"try-error"))){
        pval <- 1;hr <- 1
      } else {
        coxfit1 <- coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub)
        coxfit2 <- coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub2)
        sum <- summary(coxfit1)
        sum2 <- summary(coxfit2)
        zval <- sum$coefficients[4]
        zval2 <- sum2$coefficients[4]
        if(is.na(zval) | is.na(zval2))
          pval <- 1
        else pval <- 2*(1 - pnorm((zval - zval2)/sqrt(2)))
        hr <- sum$conf.int[2]
      }
    }
    return(c(pval=pval,diff=hr,sub.size=nsub,sub.proption=prop,sub.treat=NROW(sub.treat),
             sub.control=NROW(sub.control)))
  })
  objv2 <- do.call(rbind,objv)
  return(list(objv=objv2[,1],diff=objv2[,2],sub.size=objv2[,3],sub.proption=sprintf("%s%%",objv2[,4]),sub.treat=objv2[,5],
              sub.control=objv2[,6]))
}

## Caculate the treatment effect (that is the HR and pvalue caculated through the one-sided logrank test) for time to event endpoints
## Enhanced to make sure there will be an interactive effect: if the treatment effect in the complementory
## set of choosed subgroup less than that in the choosed subgroup, then the p-value will set to 1
objF.logrank <- function(subgroups,data,cutoff=10,trt="treat==1",ctrl="treat==2",status="e.rfs",time="t.rfs",int=T){
  require(survival)
  trtvar <- strsplit(trt,"=")[[1]][1]
  trtval <- strsplit(trt,"=")[[1]][2]
  ctrlval <- strsplit(ctrl,"=")[[1]][2]
  data$newtreat <- factor(data[,trtvar],levels = c(trtval, ctrlval))
  npop <- NROW(data)
  objv <- lapply(subgroups,function(xx){
    if (!grepl("no subgroup",xx)) {
      sub <- subset(data, eval(parse(text = get("xx"))))
      sub2 <- subset(data, !eval(parse(text = get("xx")))) #the complementory set
    } else {sub <- data; sub2 <- data}
    sub.treat <- subset(sub,eval(parse(text = sprintf("%s == %s",trtvar,trtval))))
    sub.control <- subset(sub,eval(parse(text = sprintf("%s == %s",trtvar,ctrlval))))
    nsub <- 2*min(NROW(sub.treat),NROW(sub.control))
    prop <- nsub*100/npop
    if (int & (nsub == 0 | prop < cutoff)){
      pval <- 1
      hr <- 1
    } else if (NROW(sub) == 0 | NROW(sub2) == 0) {
      pval <- 1
      hr <- 1
    } else {
      try1 <- try(coxfit1 <- coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub),T)
      try2 <- try(coxfit2 <- coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub2),T)
      try3 <- try(sum <- summary(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub)),T)
      try4 <- try(sum2 <- summary(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub2)),T)
      if (!int & (inherits(try1,"try-error") | inherits(try3,"try-error") )){
        pval <- 1;hr <- 1
      } else if (int & any(inherits(try1,"try-error"),inherits(try2,"try-error"),inherits(try3,"try-error"),inherits(try4,"try-error"))){
        pval <- 1;hr <- 1
      } else {
        zval <- sum$coefficients[4]
        if (inherits(try2,"try-error") | inherits(try4,"try-error"))
          zval2 <- NA
        else {
          zval2 <- sum2$coefficients[4]  
        }
        if (!is.na(zval)) {
          zval.lr <- sqrt(sum$sctest[1])*(-1)^(!zval>0)
          pval <- pnorm(zval.lr, lower.tail = F)
        } else pval <- 1
        if (!is.na(zval2)) {
          zval.lr2 <- sqrt(sum2$sctest[1])*(-1)^(!zval2>0)
          pval2 <- pnorm(zval.lr2, lower.tail = F)
        } else pval2 <- 1
        hr <- sum$conf.int[2]
        if (int & pval2 < pval) {pval <- 1;hr<-1} 
      }
    }
    return(c(pval=pval,diff=hr,sub.size=nsub,sub.proption=prop,sub.treat=NROW(sub.treat),
             sub.control=NROW(sub.control)))
  })
  objv2 <- do.call(rbind,objv)
  return(list(objv=objv2[,1],diff=objv2[,2],sub.size=objv2[,3],sub.proption=sprintf("%s%%",objv2[,4]),sub.treat=objv2[,5],
              sub.control=objv2[,6]))
}

## Caculate the objective value and estimated treatment effect (that is the HR caculated through the cox model)
## of subgroup in training data
## Objective value is calculated via: 2(1-Phi(|Z+ - Z-|/2^0.5)), where Phi is the CDF of standard normal distribution and
## Z+, Z- are the test statistics for testing the one-sided hypothesis of treatment efficacy in subgroup positive and
## negative population respectively.
objF.logrank2 <- function(subgroups,data,cutoff=10,trt="treat==1",ctrl="treat==2",status="e.rfs",time="t.rfs"){
  require(survival)
  trtvar <- strsplit(trt,"=")[[1]][1]
  trtval <- strsplit(trt,"=")[[1]][2]
  ctrlval <- strsplit(ctrl,"=")[[1]][2]
  data$newtreat <- factor(data[,trtvar],levels = c(trtval, ctrlval))
  npop <- NROW(data)
  objv <- lapply(subgroups,function(xx){
    if (!grepl("no subgroup",xx)) {
      sub <- subset(data, eval(parse(text = get("xx"))))
      sub2 <- subset(data, !eval(parse(text = get("xx")))) #the complementory set
    } else {sub <- data; sub2 <- data}
    sub.treat <- subset(sub,eval(parse(text = sprintf("%s == %s",trtvar,trtval))))
    sub.control <- subset(sub,eval(parse(text = sprintf("%s == %s",trtvar,ctrlval))))
    nsub <- 2*min(NROW(sub.treat),NROW(sub.control))
    prop <- nsub*100/npop
    if (nsub == 0 | prop < cutoff){
      pval <- 1
      hr <- 1
    } else if (NROW(sub) == 0 | NROW(sub2) == 0) {
      pval <- 1
      hr <- 1
    } else {
      try1 <- try(coxfit1 <- coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub),T)
      try2 <- try(coxfit2 <- coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub2),T)
      try3 <- try(sum <- summary(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub)),T)
      try4 <- try(sum2 <- summary(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat",time,status)), data = sub2)),T)
      if (any(inherits(try1,"try-error"),inherits(try2,"try-error"),inherits(try3,"try-error"),inherits(try4,"try-error"))){
        pval <- 1;hr <- 1
      } else {
        zval <- sum$coefficients[4]
        zval2 <- sum2$coefficients[4]
        if(is.na(zval) | is.na(zval2))
          pval <- 1
        else {
          zval.lr <- sqrt(sum$sctest[1])*(-1)^(!zval>0)
          zval.lr2 <- sqrt(sum2$sctest[1])*(-1)^(!zval2>0)
          pval <- 1 - pnorm((zval.lr - zval.lr2)/sqrt(2))
        }
        hr <- sum$conf.int[2]
      }
    }
    return(c(pval=pval,diff=hr,sub.size=nsub,sub.proption=prop,sub.treat=NROW(sub.treat),
             sub.control=NROW(sub.control)))
  })
  objv2 <- do.call(rbind,objv)
  return(list(objv=objv2[,1],diff=objv2[,2],sub.size=objv2[,3],sub.proption=sprintf("%s%%",objv2[,4]),sub.treat=objv2[,5],
              sub.control=objv2[,6]))
}

## Caculate the objective value of subgroup in training data
## Objective value is calculated through cox model: H(t) = H0(t)*exp(B*z), where B is the binary predictor indicating if
## signature is positive and z is the treatment indicator. Object value is the p value of B*z interaction term.
objF.cox3 <- function(subgroups,data,cutoff=10,trt="treat==1",ctrl="treat==2",status="e.rfs",time="t.rfs"){
  require(survival)
  trtvar <- strsplit(trt,"=")[[1]][1]
  trtval <- strsplit(trt,"=")[[1]][2]
  ctrlval <- strsplit(ctrl,"=")[[1]][2]
  data$newtreat <- factor(data[,trtvar],levels = c(ctrlval, trtval))
  npop <- NROW(data)
  objv <- lapply(subgroups,function(xx){
    if (!grepl("no subgroup",xx)) {
      sub <- subset(data, eval(parse(text = get("xx"))))
      sub2 <- subset(data, !eval(parse(text = get("xx")))) #the complementory set
      I <- with(data, eval(parse(text = get("xx"))))
      I[is.na(I)] <- FALSE
      predictor.true <- factor(I,levels=c(F,T))
    } else {sub <- data; sub2 <- data;predictor.true <- factor(rep(F,nrow(data)),levels=c(F,T))}
    sub.treat <- subset(sub,eval(parse(text = sprintf("%s == %s",trtvar,trtval))))
    sub.control <- subset(sub,eval(parse(text = sprintf("%s == %s",trtvar,ctrlval))))
    nsub <- 2*min(NROW(sub.treat),NROW(sub.control))
    prop <- nsub*100/npop
    if (nsub == 0 | prop < cutoff){
      pval <- 1
      hr <- 1
    } else if (NROW(sub) == 0 | NROW(sub2) == 0) {
      pval <- 1
      hr <- 1
    } else {
      newdata <- data.frame(data,predictor.true=predictor.true)
      coxfit <- try(coxph(as.formula(sprintf("Surv(%s, %s) ~ newtreat*predictor.true",time,status)), data = newdata),T)
      if (inherits(coxfit,"try-error")){
        pval <- 1;hr <- 1
      } else {
        sum <- try(summary(coxfit),T)
        if (inherits(sum,"try-error")) {
          pval <- 1;hr <- 1
        }else {
          zval <- sum$coefficients[3,4]
          pval <- sum$coefficients[3,5]
          if(is.na(pval) | is.na(zval))
            pval <- 1
          else if (zval > 0 && pval < 0.5)
            pval <- 1 - pval
          hr <- sum$conf.int[3,1]
        }
      } 
    }
    return(c(pval=pval,diff=hr,sub.size=nsub,sub.proption=prop,sub.treat=NROW(sub.treat),
             sub.control=NROW(sub.control)))
  })
  objv2 <- do.call(rbind,objv)
  return(list(objv=objv2[,1],diff=objv2[,2],sub.size=objv2[,3],sub.proption=sprintf("%s%%",objv2[,4]),sub.treat=objv2[,5],
              sub.control=objv2[,6]))
}