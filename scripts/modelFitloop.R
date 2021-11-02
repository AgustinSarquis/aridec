# load libraries
library(SoilR)
library(FME)
library(aridec)
library(dplyr)

# load single entry
db=loadEntries()
entry=db[["Arriaga2007"]]
Ct=entry$timeSeries[-1]

# load necessary functions
# One pool function
onepFit=function(Ct, time, C0=100){
  complete=data.frame(time,Ct)
  n=nrow(complete)
  if(n < 3) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::OnepModel(t=tt,k=pars[1], C0=C0, In=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=Ct))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  inipars=c(-1*C0/complete[1,2])
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower= -Inf, upper=0)
  bestMod=Func(pars=Fit$par)
  SoilRmodel=SoilR::OnepModel(t=tt,k=Fit$par[1], C0=C0, In=0)
  #AIC=2-2*log(Fit$ms)
  return(list(k=Fit$par))#, TransitTime=1/(Fit$par), AIC=AIC))
}

# Two pool parallel function
twoppFit=function(Ct, time, C0=100, inipars=c(1, 0.5, 0.5)){
  complete=data.frame(time,Ct)
  n=nrow(complete)
  if(n < 5) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::TwopParallelModel(t=tt,ks=pars[1:2], C0=C0*c(pars[3], 1-pars[3]), In=0, gam=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,3), upper=c(Inf, Inf, 1))
  bestMod=Func(pars=Fit$par)
  SoilRmodel=SoilR::TwopParallelModel(t=tt,ks=Fit$par[1:2], C0=C0*c(Fit$par[3], 1-Fit$par[3]), In=0, gam=0)
  SystemAge=systemAge(A=diag(-(Fit$par[c(1,2)])),u=matrix(c(Fit$par[3], 1-Fit$par[3]), ncol=1), a=0, q=0.5) # como q=0.5, mean es la mediana?
  TransitTime=transitTime(A=diag(-(Fit$par[c(1,2)])),u=matrix(c(Fit$par[3], 1-Fit$par[3]), ncol=1), a=0, q=0.5)
  AIC=(2*length(Fit$par))-2*log(Fit$ms)
  return(list(k1=Fit$par[1], k2=Fit$par[2], Pool1C=Fit$par[3], SystemAge=SystemAge, TransitTime=TransitTime, AIC=AIC))
}

# Two Pool Series Function
twopsFit=function(Ct, time, C0=100, inipars=c(0.5, 0.25, 0.25, 0.15)){
  complete=data.frame(time,Ct)
  n=nrow(complete)
  if(n < 5) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::TwopSeriesModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], C0=C0*c(pars[4], 1-pars[4]), In=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,4), upper=c(Inf, Inf, 1,1))
  bestMod=Func(pars=Fit$par)
  SoilRmodel=SoilR::TwopSeriesModel(t=tt,ks=Fit$par[1:2], a21=Fit$par[1]*Fit$par[3], C0=C0*c(Fit$par[4], 1-Fit$par[4]), In=0)
  SystemAge=systemAge(A=matrix(c(-Fit$par[1],Fit$par[3],0,-Fit$par[2]),ncol=2),u=matrix(c(Fit$par[4], 1-Fit$par[4]), ncol=1), a=0, q=0.5)
  TransitTime=transitTime(A=matrix(c(-Fit$par[1],Fit$par[3],0,-Fit$par[2]),ncol=2),u=matrix(c(Fit$par[4], 1-Fit$par[4]), ncol=1), a=0, q=0.5)
  AIC=(2*length(Fit$par))-2*log(Fit$ms)
  return(list(k1=Fit$par[1], k2=Fit$par[2], a21=Fit$par[3], Pool1C=Fit$par[4], SystemAge=SystemAge, TransitTime=TransitTime, AIC=AIC))
}

# Two Pool Model with Feedback Function
twopfFit=function(Ct, time, C0=100, inipars=c(1, 0.5, 0.5, 0.5, 0.3)){
  complete=data.frame(time,Ct)
  n=nrow(complete)
  if(n < 6) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::TwopFeedbackModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], a12=pars[2]*pars[4],C0=C0*c(pars[5], 1-pars[5]), In=0)
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
                                      C0=C0*c(Fit$par[5], 1-Fit$par[5]), In=0)
  AIC=(2*length(Fit$par))-2*log(Fit$ms)
  return(list(k1=-(Fit$par[1]), k2=-(Fit$par[2]), a21=(Fit$par[3]), a12=(Fit$par[4]), Pool1C=Fit$par[5], AIC=AIC))
}

# Three Pool Parallel Function
 threeppFit=function(Ct, time, C0=100, inipars=c(0.25, 0.125, 0.125, 0.125, 0.125)){
  complete=data.frame(time,Ct)
  n=nrow(complete)
  if(n < 6) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::ThreepParallelModel(t=tt,ks=pars[1:3], C0=C0*c(pars[4], pars[5], 1-sum(pars[4:5])), In=0, gam1=0, gam2=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,5), upper=c(3, 3, 3, 1, 1))
  bestMod=Func(pars=Fit$par)
  SoilRmodel=SoilR::ThreepParallelModel(t=tt,ks=Fit$par[1:3], C0=C0*c(Fit$par[4], Fit$par[5], 1-sum(Fit$par[4:5])), In=0, gam1=0, gam2=0)
  SystemAge=systemAge(A=diag(-(Fit$par[c(1:3)])),u=matrix(c(Fit$par[4],Fit$par[5], 1-Fit$par[4]-Fit$par[5]), ncol=1), a=0, q=0.5)
  TransitTime=transitTime(A=diag(-(Fit$par[c(1:3)])),u=matrix(c(Fit$par[4],Fit$par[5], 1-Fit$par[4]-Fit$par[5]), ncol=1), a=0, q=0.5)
  AIC=(2*length(Fit$par))-2*log(Fit$ms)
  return(list(k1=(Fit$par[1]), k2=(Fit$par[2]), k3=(Fit$par[3]), Pool1C=Fit$par[4], Pool2C=Fit$par[5],
    SystemAge=SystemAge, TransitTime=TransitTime, AIC=AIC))
}

# Three Pool series Function
threepsFit=function(Ct, time, C0=100, inipars=c(0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25)){
  complete=data.frame(time,Ct)
  n=nrow(complete)
  if(n < 8) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::ThreepSeriesModel(t=tt,ks=pars[1:3], a21=pars[1]*pars[4], a32=pars[2]*pars[5], C0=C0*c(pars[6], pars[7], 1-sum(pars[6:7])), In=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,7), upper=c(Inf, Inf, Inf, 1, 1, 1,1))
  bestMod=Func(pars=Fit$par)
  SoilRmodel=SoilR::ThreepSeriesModel(t=tt,ks=Fit$par[1:3], a21=Fit$par[1]*Fit$par[4], a32=Fit$par[2]*Fit$par[5], C0=C0*c(Fit$par[6], Fit$par[7], 1-sum(Fit$par[6:7])), In=0)
  SystemAge=systemAge(A=matrix(c(-Fit$par[1], Fit$par[4], 0, 0, -Fit$par[2], Fit$par[5], 0, 0, -Fit$par[3]), ncol=3) ,u=matrix(c(Fit$par[6],Fit$par[7], 1-Fit$par[6]-Fit$par[7]), ncol=1), a=0, q=0.5)
  TransitTime=transitTime(A=matrix(c(-Fit$par[1], Fit$par[4], 0, 0, -Fit$par[2], Fit$par[5], 0, 0, -Fit$par[3]), ncol=3),u=matrix(c(Fit$par[6],Fit$par[7], 1-Fit$par[6]-Fit$par[7]), ncol=1), a=0, q=0.5)
  AIC=(2*length(Fit$par))-2*log(Fit$ms)
  return(list(k1=-(Fit$par[1]), k2=-(Fit$par[2]), k3=-(Fit$par[3]), a21=Fit$par[4], a32=Fit$par[5], Pool1C=Fit$par[6], Pool2C=Fit$par[7],
     SystemAge=SystemAge, TransitTime=TransitTime, AIC=AIC))
}

##########################################################################################################################################################################################

onePloop= function(Ct, time) {lapply(Ct, onepFit, time)}

entry=db[[107]]
df=entry$timeSeries

df=mutate(df, Time=Time*7) # transform weeks to days
df=mutate(df, Time=Time*30) # transform months to days
df=mutate(df, Time=Time*365) # transform years to days

Ct=df[-1]
time=df[,1]
write.csv(onePloop(Ct, time),"~/k x days/k107.csv")

# for entries with NAs
entrynoNA=entry$timeSeries[,c(1,5)]
tSnoNA=entrynoNA[complete.cases(entrynoNA),]
write.csv(onepFit(time=tSnoNA$Time, Ct=tSnoNA[,2]),"~/k171DR.csv")

onePoutput=lapply(Ct, onepFit, time=entry$timeSeries[,1])
twoPPoutput=lapply(Ct, twoppFit, time=entry$timeSeries[,1])
twoPSoutput=lapply(Ct, twopsFit, time=entry$timeSeries[,1])
twoPFoutput=lapply(Ct, twopfFit, time=entry$timeSeries[,1]) # agregando transit time y system age el vector se hace muy grande e imposible de procesar
#SystemAge=systemAge(A=matrix(c(-Fit$par[1],Fit$par[3],Fit$par[4],-Fit$par[2]),ncol=2),u=matrix(c(Fit$par[5], 1-Fit$par[5]), ncol=1), a=0, q=0.5)
#TransitTime=transitTime(A=matrix(c(-Fit$par[1],Fit$par[3],Fit$par[4],-Fit$par[2]),ncol=2),u=matrix(c(Fit$par[5], 1-Fit$par[5]), ncol=1), a=0, q=0.5)
threePPoutput=lapply(Ct, threeppFit, time=entry$timeSeries[,1])
threePSoutput=lapply(Ct, threepsFit, time=entry$timeSeries[,1]) # modifique los inipars
