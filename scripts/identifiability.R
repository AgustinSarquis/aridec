# load libraries
library(SoilR)
library(FME)
library(aridec)

# load single entry
db=loadEntries()
entry=db[["Arriaga2007"]]
Ct=entry$timeSeries[-1]

# 2 pool parallel model 
inipars=c(1, 0.5, 0.5)
Ct=entry$timeSeries[,1:2]
colnames(Ct)=c("time", "Ct")
tt=seq(from=0, to=tail(entry$timeSeries[,1],1), length.out = 500)
  
Func=function(pars){
    mod=SoilR::TwopParallelModel(t=tt,ks=pars[1:2], C0=100*c(pars[3], 1-pars[3]), In=0, gam=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }

costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=Ct, x="time")) 
  }

costFunc(inipars)$model # suma de cuadrados residuales

Sfun <- sensFun(costFunc, inipars)  # funciones de sensibilidad
summary(Sfun)
plot(Sfun)
pairs(Sfun)

# cuales parametros pueden ser estimados simultaneamente?
collin(Sfun) # colinearidad

# cual es el set de parametros optimos con el mejor ajuste?
Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,3), upper=c(Inf, Inf, 1))
bestMod=Func(pars=Fit$par)

# cual es la distribucion de probabilidad de los parametros? 
MCMC <- modMCMC(f = costFunc, p = Fit$par, niter = 10000, var0 = Fit$var_ms_unweighted, wvar0 = 0.1, updatecov = 50, lower=rep(0,3), upper=c(Inf, Inf, 1))
summary(MCMC)
plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000)

###########################################################################################
# 2 pool series model
inipars=c(0.5, 0.25, 0.25, 0.15)
Ct=entry$timeSeries[,1:2]
colnames(Ct)=c("time", "Ct")
tt=seq(from=0, to=tail(entry$timeSeries[,1],1), length.out = 500)

Func=function(pars){
  mod=SoilR::TwopSeriesModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], C0=100*c(pars[4], 1-pars[4]), In=0)
  Ct=SoilR::getC(mod)
  return(data.frame(time=tt, Ct=rowSums(Ct)))
}

costFunc=function(pars){
  output=Func(pars)
  return(modCost(model=output, obs=Ct, x="time")) 
}

costFunc(inipars)$model # suma de cuadrados residuales

Sfun <- sensFun(costFunc, inipars)  # funciones de sensibilidad
summary(Sfun)
plot(Sfun)
pairs(Sfun)

# cuales parametros pueden ser estimados simultaneamente?
collin(Sfun) # colinearidad

# cual es el set de parametros optimos con el mejor ajuste?
Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,4), upper=c(Inf, Inf, 1,1))
bestMod=Func(pars=Fit$par)

# cual es la distribucion de probabilidad de los parametros? 
MCMC <- modMCMC(f = costFunc, p = Fit$par, niter = 10000, var0 = Fit$var_ms_unweighted, wvar0 = 0.1, updatecov = 50, lower=rep(0,4), upper=c(Inf, Inf, 1,1))
summary(MCMC)
plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000)

###########################################################################################
# 2 pool model with feedback
inipars=c(1, 0.5, 0.5, 0.5, 0.3)
Ct=entry$timeSeries[,1:2]
colnames(Ct)=c("time", "Ct")
tt=seq(from=0, to=tail(entry$timeSeries[,1],1), length.out = 500)

Func=function(pars){
  mod=SoilR::TwopFeedbackModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], a12=pars[2]*pars[4],C0=100*c(pars[5], 1-pars[5]), In=0)
  Ct=SoilR::getC(mod)
  return(data.frame(time=tt, Ct=rowSums(Ct)))
}

costFunc=function(pars){
  output=Func(pars)
  return(modCost(model=output, obs=Ct, x="time")) 
}

costFunc(inipars)$model # suma de cuadrados residuales

Sfun <- sensFun(costFunc, inipars)  # funciones de sensibilidad
summary(Sfun)
plot(Sfun)
pairs(Sfun)

# cuales parametros pueden ser estimados simultaneamente?
collin(Sfun) # colinearidad

# cual es el set de parametros optimos con el mejor ajuste?
Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,5), upper=c(Inf, Inf, 1,1,1))
bestMod=Func(pars=Fit$par)

# cual es la distribucion de probabilidad de los parametros? 
MCMC <- modMCMC(f = costFunc, p = Fit$par, niter = 10000, var0 = Fit$var_ms_unweighted, wvar0 = 0.1, updatecov = 50, lower=rep(0,5), upper=c(Inf, Inf, 1,1,1))
summary(MCMC)
plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000)

#############################################################################################
# 3 pool parallel model 
inipars=c(0.25, 0.125, 0.125, 0.125, 0.125)
Ct=entry$timeSeries[,1:2]
colnames(Ct)=c("time", "Ct")
tt=seq(from=0, to=tail(entry$timeSeries[,1],1), length.out = 500)

Func=function(pars){
  mod=SoilR::ThreepParallelModel(t=tt,ks=pars[1:3], C0=100*c(pars[4], pars[5], 1-sum(pars[4:5])), In=0, gam1=0, gam2=0)
  Ct=SoilR::getC(mod)
  return(data.frame(time=tt, Ct=rowSums(Ct)))
}

costFunc=function(pars){
  output=Func(pars)
  return(modCost(model=output, obs=Ct, x="time")) 
}

costFunc(inipars)$model # suma de cuadrados residuales

Sfun <- sensFun(costFunc, inipars)  # funciones de sensibilidad
summary(Sfun)
plot(Sfun)
pairs(Sfun)

# cuales parametros pueden ser estimados simultaneamente?
collin(Sfun) # colinearidad

# cual es el set de parametros optimos con el mejor ajuste?
Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,5), upper=c(Inf, Inf, 1,1,1))
bestMod=Func(pars=Fit$par)

# cual es la distribucion de probabilidad de los parametros? 
MCMC <- modMCMC(f = costFunc, p = Fit$par, niter = 10000, var0 = Fit$var_ms_unweighted, wvar0 = 0.1, updatecov = 50, lower=rep(0,5), upper=c(Inf, Inf, 1,1,1))
summary(MCMC)
plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000)

#############################################################################################
# 3 pool series model 
inipars=c(0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25)
Ct=entry$timeSeries[,1:2]
colnames(Ct)=c("time", "Ct")
tt=seq(from=0, to=tail(entry$timeSeries[,1],1), length.out = 500)

Func=function(pars){
  mod=SoilR::ThreepSeriesModel(t=tt,ks=pars[1:3], a21=pars[1]*pars[4], a32=pars[2]*pars[5], C0=100*c(pars[6], pars[7], 1-sum(pars[6:7])), In=0)
  Ct=SoilR::getC(mod)
  return(data.frame(time=tt, Ct=rowSums(Ct)))
}

costFunc=function(pars){
  output=Func(pars)
  return(modCost(model=output, obs=Ct, x="time")) 
}

costFunc(inipars)$model # suma de cuadrados residuales

Sfun <- sensFun(costFunc, inipars)  # funciones de sensibilidad
summary(Sfun)
plot(Sfun)
pairs(Sfun)

# cuales parametros pueden ser estimados simultaneamente?
getOption("max.print")
collin(Sfun) # colinearidad

# cual es el set de parametros optimos con el mejor ajuste?
Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,7), upper=c(Inf, Inf, 1,1,1,1,1))
bestMod=Func(pars=Fit$par)

# cual es la distribucion de probabilidad de los parametros? 
MCMC <- modMCMC(f = costFunc, p = Fit$par, niter = 10000, var0 = Fit$var_ms_unweighted, wvar0 = 0.1, updatecov = 50, lower=rep(0,7), upper=c(Inf, Inf, 1,1,1,1,1))
summary(MCMC)
plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000)