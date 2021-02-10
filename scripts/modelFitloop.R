# load libraries
library(SoilR)
library(FME)
library(aridec)

# load single entry
db=loadEntries()
deco=db[["Arriaga2007"]]

# load necessary functions
# One pool function
onepFit=function(timeSeries, initialCarbon=100){
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
  SoilRmodel=SoilR::OnepModel(t=tt,k=Fit$par[1], C0=initialCarbon, In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}

# Two pool parallel function
twoppFit=function(timeSeries, initialCarbon=100, inipars=c(1, 0.5, 0.5)){
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
  SoilRmodel=SoilR::TwopParallelModel(t=tt,ks=Fit$par[1:2], C0=initialCarbon*c(Fit$par[3], 1-Fit$par[3]), In=0, gam=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}

# Two Pool Series Function
twopsFit=function(timeSeries, initialCarbon=100, inipars=c(0.5, 0.25, 0.25, 0.15)){
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
  SoilRmodel=SoilR::TwopSeriesModel(t=tt,ks=Fit$par[1:2], a21=Fit$par[1]*Fit$par[3], C0=initialCarbon*c(Fit$par[4], 1-Fit$par[4]), In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}

# Two Pool Model with Feedback Function
twopfFit=function(timeSeries, initialCarbon=100, inipars=c(1, 0.5, 0.5, 0.5, 0.3)){
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
  SoilRmodel=SoilR::TwopFeedbackModel(t=tt,ks=Fit$par[1:2], a21=Fit$par[1]*Fit$par[3], a12=Fit$par[2]*Fit$par[4],
                                      C0=initialCarbon*c(Fit$par[5], 1-Fit$par[5]), In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}

# Three Pool Parallel Function
threeppFit=function(timeSeries, initialCarbon=100, inipars=c(0.25, 0.125, 0.125, 0.125, 0.125)){
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
  SoilRmodel=SoilR::ThreepParallelModel(t=tt,ks=Fit$par[1:3], C0=initialCarbon*c(Fit$par[4], Fit$par[5], 1-sum(Fit$par[4:5])), In=0, gam1=0, gam2=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}

# Three Pool series Function
threepsFit=function(timeSeries, initialCarbon=100, inipars=c(0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25)){
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
  SoilRmodel=SoilR::ThreepSeriesModel(t=tt,ks=Fit$par[1:3], a21=Fit$par[1]*Fit$par[4], a32=Fit$par[2]*Fit$par[5], C0=initialCarbon*c(Fit$par[6], Fit$par[7], 1-sum(Fit$par[6:7])), In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}

##########################################################################################################################################################################################

Arriaga2007=list()
M1list=list()
M1output=list()
M2list=list()
M2output=list()
M3list=list()
M3output=list()
n=ncol(deco$timeSeries)
for (i in 2:n){
  names=paste(colnames(deco$timeSeries[i]))
  M1=onepFit(deco$timeSeries[,c(1,i)])
  M1list[[names]]=M1output[[i]]=list("k"=M1$FMEmodel$par, "AIC"=2-2*log(M1$FMEmodel$ms), "System age"=-1/M1$FMEmodel$par)
  M2=twoppFit(deco$timeSeries[,c(1,i)])
  M2list[[names]]=M2output[[i]]=list("k1"=-(M2$FMEmodel$par[1]),
                                     "k2"=-(M2$FMEmodel$par[2]),
                                     "Proportion of C in pool 1"=M2$FMEmodel$par[3],
                                     "AIC"=(2*length(M2$FMEmodel$par))-2*log(M2$FMEmodel$ms),
                                     "System age"=systemAge(A=diag(-(M2$FMEmodel$par[c(1,2)])),u=matrix(c(M2$FMEmodel$par[3], 1-M2$FMEmodel$par[3]), ncol=1)),
                                     "Transit Time"=transitTime(A=diag(-(M2$FMEmodel$par[c(1,2)])),u=matrix(c(M2$FMEmodel$par[3], 1-M2$FMEmodel$par[3]), ncol=1))
                                     )
  M3=twopsFit(deco$timeSeries[,c(1,i)]) # le cambie los inipars
  M3list[[names]]=M3output[[i]]=list("k1"=-(M3$FMEmodel$par[1]),
                                     "k2"=-(M3$FMEmodel$par[2]),
                                     "a21"=(M3$FMEmodel$par[3]),
                                     "Proportion of C in pool 1"=M3$FMEmodel$par[4],
                                     "AIC"=(2*length(M3$FMEmodel$par))-2*log(M3$FMEmodel$ms),
                                     "System age"=systemAge(A=matrix(c(-M3$FMEmodel$par[1],M3$FMEmodel$par[3],0,-M3$FMEmodel$par[2]),ncol=2),u=matrix(c(M3$FMEmodel$par[4], 1-M3$FMEmodel$par[4]), ncol=1)),
                                     "Transit time"=transitTime(A=matrix(c(-M3$FMEmodel$par[1],M3$FMEmodel$par[3],0,-M3$FMEmodel$par[2]),ncol=2),u=matrix(c(M3$FMEmodel$par[4], 1-M3$FMEmodel$par[4]), ncol=1))
                                     )
  Arriaga2007=list("One pool"=M1list,"Two pools parallel"=M2list,"Two pools series"=M3list)
}

M4list=list()
M4output=list()
for (i in 2:n){ # no funciona: "Error in 1:s : el resultado seria un vector muy largo"
  names=paste(colnames(deco$timeSeries[i]))
  M4=twopfFit(deco$timeSeries[,c(1,i)])
  M4list[[names]]=M4output[[i]]=list("k1"=-(M4$FMEmodel$par[1]),
                                   "k2"=-(M4$FMEmodel$par[2]),
                                   "a21"=(M4$FMEmodel$par[3]),
                                   "a12"=(M4$FMEmodel$par[4]),
                                   "Proportion of C in pool 1"=M4$FMEmodel$par[5],
                                   "AIC"=(2*length(M4$FMEmodel$par))-2*log(M4$FMEmodel$ms),
                                   "System age"= systemAge(A=matrix(c(-M4$FMEmodel$par[1],M4$FMEmodel$par[3],M4$FMEmodel$par[4],-M4$FMEmodel$par[2]),ncol=2),u=matrix(c(M4$FMEmodel$par[5], 1-M4$FMEmodel$par[5]), ncol=1)),
                                   "Transit time"=transitTime(A=matrix(c(-M4$FMEmodel$par[1],M4$FMEmodel$par[3],M4$FMEmodel$par[4],-M4$FMEmodel$par[2]),ncol=2),u=matrix(c(M4$FMEmodel$par[5], 1-M4$FMEmodel$par[5]), ncol=1))
)
}

# M5 me daba valores de C0 para los pools 1 y 2 cercanos a 1. Al cambiar inipars se soluciono en casi todos los casos
M5list=list()
M5output=list()
for (i in 2:n){ 
 names=paste(colnames(deco$timeSeries[i]))
 M5=threeppFit(deco$timeSeries[,c(1,i)])
 M5list[[names]]=M5output[[i]]=list("k1"=-(M5$FMEmodel$par[1]),
                                   "k2"=-(M5$FMEmodel$par[2]),
                                   "k3"=-(M5$FMEmodel$par[3]),
                                   "Proportion of C in pool 1"=M5$FMEmodel$par[4],
                                   "Proportion of C in pool 2"=M5$FMEmodel$par[5],
                                   "AIC"=(2*length(M5$FMEmodel$par))-2*log(M5$FMEmodel$ms),
                                   "System age"=systemAge(A=diag(-(M5$FMEmodel$par[c(1:3)])),u=matrix(c(M5$FMEmodel$par[4],M5$FMEmodel$par[5], 1-M5$FMEmodel$par[4]-M5$FMEmodel$par[5]), ncol=1)),
                                   "Transit Time"=transitTime(A=diag(-(M5$FMEmodel$par[c(1:3)])),u=matrix(c(M5$FMEmodel$par[4],M5$FMEmodel$par[5], 1-M5$FMEmodel$par[4]-M5$FMEmodel$par[5]), ncol=1))
)
}

# M6 me genera 6 modelos donde la proporcion total de C es mayor que 1. La baje los inipars porque no corria, pero si se los bajo mas tampoco corre. 
M6output=list()
M6list=list()
for (i in 2:n){ 
  names=paste(colnames(deco$timeSeries[i]))
  M6=threepsFit(deco$timeSeries[,c(1,i)])
  M6list[[names]]=M6output[[i]]=list("k1"=-(M6$FMEmodel$par[1]),
                                     "k2"=-(M6$FMEmodel$par[2]),
                                     "k3"=-(M6$FMEmodel$par[3]),
                                     "a21"=M6$FMEmodel$par[4],
                                     "a32"=M6$FMEmodel$par[5],
                                     "Proportion of C in pool 1"=M6$FMEmodel$par[6],
                                     "Proportion of C in pool 2"=M6$FMEmodel$par[7],
                                     "AIC"=(2*length(M6$FMEmodel$par))-2*log(M6$FMEmodel$ms),
                                     "System age"=systemAge(A=matrix(c(-M6$FMEmodel$par[1], M6$FMEmodel$par[4], 0, 0, -M6$FMEmodel$par[2], M6$FMEmodel$par[5], 0, 0, -M6$FMEmodel$par[3]), ncol=3) ,u=matrix(c(M6$FMEmodel$par[6],M6$FMEmodel$par[7], 1-M6$FMEmodel$par[6]-M6$FMEmodel$par[7]), ncol=1)),
                                     "Transit Time"=transitTime(A=matrix(c(-M6$FMEmodel$par[1], M6$FMEmodel$par[4], 0, 0, -M6$FMEmodel$par[2], M6$FMEmodel$par[5], 0, 0, -M6$FMEmodel$par[3]), ncol=3),u=matrix(c(M6$FMEmodel$par[6],M6$FMEmodel$par[7], 1-M6$FMEmodel$par[6]-M6$FMEmodel$par[7]), ncol=1))
  )
  }





