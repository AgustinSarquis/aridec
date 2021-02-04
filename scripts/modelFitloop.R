# load libraries
library(SoilR)
library(FME)
library(aridec)

# load single entry
db=loadEntries()
deco=db[["Arriaga2007"]]

# load necessary functions
# One pool function
onepFit=function(timeSeries, initialCarbon){
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  if(n < 3) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::OnepModel(t=tt,k=pars[1], C0=initialCarbon, In=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=Ct))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  inipars=c(-1*initialCarbon/complete[1,2])
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower= -Inf, upper=0)
  bestMod=Func(pars=Fit$par)
  print(paste("Best fit parameter: ",Fit$par))
  AIC=2-2*log(Fit$ms)
  print(paste("AIC = ",AIC))
  SoilRmodel=SoilR::OnepModel(t=tt,k=Fit$par[1], C0=initialCarbon, In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}

# Two pool parallel function
twoppFit=function(timeSeries, initialCarbon, inipars=c(1, 0.5, 0.5)){
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  if(n < 5) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::TwopParallelModel(t=tt,ks=pars[1:2], C0=initialCarbon*c(pars[3], 1-pars[3]), In=0, gam=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,3), upper=c(Inf, Inf, 1))
  bestMod=Func(pars=Fit$par)
  print(paste(c("k1=", "k2=", "proportion of C0 in pool 1="),Fit$par))
  AIC=(2*length(Fit$par))-2*log(Fit$ms)
  print(paste("AIC = ",AIC))
  SoilRmodel=SoilR::TwopParallelModel(t=tt,ks=Fit$par[1:2], C0=initialCarbon*c(Fit$par[3], 1-Fit$par[3]), In=0, gam=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}

# Two Pool Series Function
twopsFit=function(timeSeries, initialCarbon, inipars=c(1, 0.5, 0.5, 0.3)){
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  if(n < 5) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::TwopSeriesModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], C0=initialCarbon*c(pars[4], 1-pars[4]), In=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,4), upper=c(Inf, Inf, 1,1))
  bestMod=Func(pars=Fit$par)
  print(paste(c("k1=", "k2=", "a21=", "Proportion of C0 in pool 1="),Fit$par))
  AIC=(2*length(Fit$par))-2*log(Fit$ms)
  print(paste("AIC = ",AIC))
  SoilRmodel=SoilR::TwopSeriesModel(t=tt,ks=Fit$par[1:2], a21=Fit$par[1]*Fit$par[3], C0=initialCarbon*c(Fit$par[4], 1-Fit$par[4]), In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}

# Two Pool Model with Feedback Function
twopfFit=function(timeSeries, initialCarbon, inipars=c(1, 0.5, 0.5, 0.5, 0.3)){
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  if(n < 6) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::TwopFeedbackModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], a12=pars[2]*pars[4],C0=initialCarbon*c(pars[5], 1-pars[5]), In=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0, 5), upper=c(Inf, Inf, 1, 1, 1))
  bestMod=Func(pars=Fit$par)
  print(paste(c("k1=", "k2=", "a21=", "a12=", "Proportion of C0 in pool 1="),Fit$par))
  AIC=(2*length(Fit$par))-2*log(Fit$ms)
  print(paste("AIC = ",AIC))
  SoilRmodel=SoilR::TwopFeedbackModel(t=tt,ks=Fit$par[1:2], a21=Fit$par[1]*Fit$par[3], a12=Fit$par[2]*Fit$par[4],
                                      C0=initialCarbon*c(Fit$par[5], 1-Fit$par[5]), In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}

# Three Pool Parallel Function
threeppFit=function(timeSeries, initialCarbon, inipars=c(1, 0.5, 0.5, 0.5, 0.5)){
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  if(n < 6) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::ThreepParallelModel(t=tt,ks=pars[1:3], C0=initialCarbon*c(pars[4], pars[5], 1-sum(pars[4:5])), In=0, gam1=0, gam2=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,5), upper=c(3, 3, 3, 1, 1))
  bestMod=Func(pars=Fit$par)
  print(paste(c("k1=", "k2=", "k3=","proportion of C0 in pool 1=", "Proportion of C0 in pool 2="),Fit$par[1:5]))
  AIC=(2*length(Fit$par))-2*log(Fit$ms)
  print(paste("AIC = ",AIC))
  SoilRmodel=SoilR::ThreepParallelModel(t=tt,ks=Fit$par[1:3], C0=initialCarbon*c(Fit$par[4], Fit$par[5], 1-sum(Fit$par[4:5])), In=0, gam1=0, gam2=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}

# Three Pool series Function
threepsFit=function(timeSeries, initialCarbon, inipars=c(1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)){
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  if(n < 8) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::ThreepSeriesModel(t=tt,ks=pars[1:3], a21=pars[1]*pars[4], a32=pars[2]*pars[5], C0=initialCarbon*c(pars[6], pars[7], 1-sum(pars[6:7])), In=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,7), upper=c(Inf, Inf, Inf, 1, 1, 1,1))
  bestMod=Func(pars=Fit$par)
  print(paste(c("k1=", "k2=", "k3=","alpha21", "alpha32","proportion of C0 in pool 1=", "Proportion of C0 in pool 2="),Fit$par[1:7]))
  AIC=(2*length(Fit$par))-2*log(Fit$ms)
  print(paste("AIC = ",AIC))
  SoilRmodel=SoilR::ThreepSeriesModel(t=tt,ks=Fit$par[1:3], a21=Fit$par[1]*Fit$par[4], a32=Fit$par[2]*Fit$par[5], C0=initialCarbon*c(Fit$par[6], Fit$par[7], 1-sum(Fit$par[6:7])), In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}

##########################################################################################################################################################################################
M1output=list()
M1list=list()
n=ncol(deco$timeSeries)
for (i in 2:n){
  names=paste(colnames(deco$timeSeries[i]))
  M1=onepFit(deco$timeSeries[,c(1,i)], initialCarbon = 100)
  M1list[[names]]=M1output[[i]]=list("k"=M1$FMEmodel$par, "AIC"=2-2*log(M1$FMEmodel$ms), "SA"=-1/M1$FMEmodel$par)
   }

# Run a Two Pool Parallel Model
M2=twoppFit(deco$timeSeries[,c(1,2)], initialCarbon = 100) 
k=M2$FMEmodel$par[c(1,2)]
A2=diag(-k)
u2=matrix(c(M2$FMEmodel$par[3], 1-M2$FMEmodel$par[3]), ncol=1)
SA2=systemAge(A=A2, u=u2)
TT2=transitTime(A=A2, u=u2)

# Run a two Pool Series Model
M3=twopsFit(deco$timeSeries[,c(1,2)], initialCarbon = 100)
k=M3$FMEmodel$par[c(1,2)]
A3=diag(-k)
A3[2,1]=M3$FMEmodel$par[3]
u3=matrix(c(M3$FMEmodel$par[4], 1-M3$FMEmodel$par[4]), ncol=1)
SA3=systemAge(A=A3, u=u3)
TT3=transitTime(A=A3, u=u3)

# Run a Two Pool Model with feedback
M4=twopfFit(deco$timeSeries[,c(1,2)], initialCarbon = 100)
k=M4$FMEmodel$par[c(1,2)]
A4=diag(-k)
A4[2,1]=M4$FMEmodel$par[3]
A4[1,2]=M4$FMEmodel$par[4]
u4=matrix(c(M4$FMEmodel$par[5], 1-M4$FMEmodel$par[5]), ncol=1)
SA4=systemAge(A=A4, u=u4)
TT4=transitTime(A=A4, u=u4)

# Run a Three Pool Parallel Model
M5=threeppFit(deco$timeSeries[,c(1,2)], initialCarbon = 100)
k=M5$FMEmodel$par[c(1,2,3)]
A5=diag(-k)
u5=matrix(c(M5$FMEmodel$par[4], M5$FMEmodel$par[5], 1-M5$FMEmodel$par[4]-M5$FMEmodel$par[5]), ncol=1)
SA5=systemAge(A=A5, u=u5)
TT5=transitTime(A=A5, u=u5)

# Run a Three Pool Series Model
M6=threepsFit(deco$timeSeries[,c(1,2)], initialCarbon = 100)
k=M6$FMEmodel$par[c(1,2,3)]
A6=diag(-k)
A6[2,1]=M6$FMEmodel$par[4]
A6[3,2]=M6$FMEmodel$par[5]
u6=matrix(c(M6$FMEmodel$par[6], M6$FMEmodel$par[7], 1-M6$FMEmodel$par[6]-M6$FMEmodel$par[7]), ncol=1)
SA6=systemAge(A=A6, u=u6)
TT6=transitTime(A=A6, u=u6)




