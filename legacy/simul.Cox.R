simul.Cox <- function(n,exp,weib.shape,weib.scale,censor.par=NULL,extra.noise=0,seed=123){
  library(mvtnorm)
  exp<-str2expression(exp)
  set.seed(seed)
  x1 <- rnorm(n);
  x<-rmvnorm(n,sigma = matrix(c(1,0.5,0.5,1), ncol=2))
  x2<-x[,1];x3<-x[,2]
  x4<-rexp(n)
  x5<-rbinom(n,1,0.5)
  x6<-rbinom(n,1,0.5)
  # x6<-apply(rmultinom(n,1,rep(0.1,10)),2,which.max)
  sigma <- matrix(rep(0.5,16), ncol=4)
  diag(sigma) <- 1
  x<-rmvnorm(n,sigma = sigma)
  x7<-x[,1];x8<-x[,2];x9 <- x[,3];x10<-x[,4]
  z<-rbinom(n,1,0.5)
  u<-runif(n,0,1)
  if (extra.noise){
    if (extra.noise %% 10) stop(sprintf("extra noise %s should be divisible by 10",extra.noise))
    stks <- extra.noise/10
    extra.x <- do.call(cbind,lapply(1:stks,function(stk){
      set.seed(seed+stk)
      sigma2 <- matrix(rep(0.5,100), ncol=10)
      diag(sigma2) <- 1
      xx <- rmvnorm(n,sigma = sigma2)
      colnames(xx) <- paste("x",(10*stk+1):(10*stk+10),sep="")
      return(xx)
    }))
  }
  RR <- exp(eval(exp))
  Ts <- weib.scale*(-1*log(u)*RR)^(1/weib.shape)
  if (!is.null(censor.par)) {
    Tc <- runif(n,0,censor.par)
    event <- ifelse(Tc > Ts, 1, 0)
  } else event <- rep(1,n)
  if (!extra.noise)
    data <- data.frame(Ts=Ts,event=event,z=z,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,
                       x6=x6,x7=x7,x8=x8,x9=x9,x10=x10)
  else data <- data.frame(Ts=Ts,event=event,z=z,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,
                            x6=x6,x7=x7,x8=x8,x9=x9,x10=x10,extra.x)
  data$x5 <- as.factor(data$x5)
  data$x6 <- as.factor(data$x6)
  data$z <- factor(data$z,levels = c(1,0))
  return(data)
}


