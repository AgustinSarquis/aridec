# Here you can find a procedure to answer the question:
# "What types of model structures can I fit with these data?"
# (modified from Soetaert & Petzoldt, 2010)

# load libraries
library(SoilR)
library(FME)
library(aridec)
library(dplyr)

# load database and select a single entry
db=loadEntries("~/aridec/data/") # specify aridec route in your computer
entry=db$Abera2014 # replace Henry2008 with the specific CitationKey
Ct=entry$timeSeries[,c(1,17)] # select variables Time and V2 (or any variable of your interest)
colnames(Ct)=c("time", "Ct")

# we recommend working with monthly time units
entry$variables$V1$units # check time units
Ct=mutate(Ct, time=time/30) # transform days to months
Ct=mutate(Ct, time=time*12) # transform years to months
Ct=mutate(Ct, time=time/4) # transform weeks to months

###########################################################################################
# Can I fit a 2 pool parallel model?

inipars=c(1, 0.5, 0.5)
tt=seq(from=0, to=tail(Ct[,1],1), length.out = 500)
Func=function(pars){
  mod=SoilR::TwopParallelModel(t=tt,ks=pars[1:2], C0=100*c(pars[3], 1-pars[3]), In=0, gam=0)
  Ct=SoilR::getC(mod)
  return(data.frame(time=tt, Ct=rowSums(Ct)))
}
costFunc=function(pars){
  output=Func(pars)
  return(modCost(model=output, obs=Ct, x="time")) 
}

costFunc(inipars)$model # Get sums of squared residuals

Sfun <- sensFun(costFunc, inipars)  # sensibility functions
summary(Sfun)
plot(Sfun)
pairs(Sfun)

# Get collinearity indexes for each combination of parameters
# collinearity values under 20 are considered acceptable
# you can fit model structures with this data set if collin<20
write.csv(collin(Sfun), '~/Documents/collinearity/Abera2014_2pp.csv') #save it in your computer as .csv file

# What is the optimal set of parameters with the best fit?
Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,3), upper=c(Inf, Inf, 1))
bestMod=Func(pars=Fit$par)

# what is the probability distribution of the parameters? 
MCMC <- modMCMC(f = costFunc, p = Fit$par, niter = 10000, var0 = Fit$var_ms_unweighted, wvar0 = 0.1, updatecov = 50, lower=rep(0,3), upper=c(Inf, Inf, 1))
summary(MCMC)
plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000)

###########################################################################################
# Can I fit a 2 pool series model?
inipars=c(0.5, 0.25, 0.25, 0.15)
tt=seq(from=0, to=tail(Ct[,1],1), length.out = 500)
Func=function(pars){
  mod=SoilR::TwopSeriesModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], C0=100*c(pars[4], 1-pars[4]), In=0)
  Ct=SoilR::getC(mod)
  return(data.frame(time=tt, Ct=rowSums(Ct)))
}
costFunc=function(pars){
  output=Func(pars)
  return(modCost(model=output, obs=Ct, x="time")) 
}

costFunc(inipars)$model # sums of squared residuals

Sfun <- sensFun(costFunc, inipars)  # sensibility functions
summary(Sfun)
plot(Sfun)
pairs(Sfun)

# Get collinearity indexes for each combination of parameters
# collinearity values under 20 are considered acceptable
# you can fit model structures with this data set if collin<20
write.csv(collin(Sfun), '~/Documents/collinearity/Abera2014_2ps.csv')

# which set of parameters has the best fit?
Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,4), upper=c(Inf, Inf, 1,1))
bestMod=Func(pars=Fit$par)

# what is the probability distribution of the parameters? 
MCMC <- modMCMC(f = costFunc, p = Fit$par, niter = 10000, var0 = Fit$var_ms_unweighted, wvar0 = 0.1, updatecov = 50, lower=rep(0,4), upper=c(Inf, Inf, 1,1))
summary(MCMC)
plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000)

###########################################################################################
# Can I fit a 2 pool model with feedback?
inipars=c(1, 0.5, 0.5, 0.5, 0.3)
tt=seq(from=0, to=tail(Ct[,1],1), length.out = 500)
Func=function(pars){
  mod=SoilR::TwopFeedbackModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], a12=pars[2]*pars[4],C0=100*c(pars[5], 1-pars[5]), In=0)
  Ct=SoilR::getC(mod)
  return(data.frame(time=tt, Ct=rowSums(Ct)))
}
costFunc=function(pars){
  output=Func(pars)
  return(modCost(model=output, obs=Ct, x="time")) 
}

costFunc(inipars)$model # sums of squared residuals

Sfun <- sensFun(costFunc, inipars)  # sensibility functions
summary(Sfun)
plot(Sfun)
pairs(Sfun)

# Get collinearity indexes for each combination of parameters
# collinearity values under 20 are considered acceptable
# you can fit model structures with this data set if collin<20
collin(Sfun) 

# which set of parameters has the best fit?
Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,5), upper=c(Inf, Inf, 1,1,1))
bestMod=Func(pars=Fit$par)

# what is the probability distribution of the parameters? 
MCMC <- modMCMC(f = costFunc, p = Fit$par, niter = 10000, var0 = Fit$var_ms_unweighted, wvar0 = 0.1, updatecov = 50, lower=rep(0,5), upper=c(Inf, Inf, 1,1,1))
summary(MCMC)
plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000)

#############################################################################################
# Can I fit a 3 pool parallel model?
inipars=c(0.25, 0.125, 0.125, 0.125, 0.125)
tt=seq(from=0, to=tail(Ct[,1],1), length.out = 500)
Func=function(pars){
  mod=SoilR::ThreepParallelModel(t=tt,ks=pars[1:3], C0=100*c(pars[4], pars[5], 1-sum(pars[4:5])), In=0, gam1=0, gam2=0)
  Ct=SoilR::getC(mod)
  return(data.frame(time=tt, Ct=rowSums(Ct)))
}
costFunc=function(pars){
  output=Func(pars)
  return(modCost(model=output, obs=Ct, x="time")) 
}

costFunc(inipars)$model # sums of squared residuals

Sfun <- sensFun(costFunc, inipars)  # sensibility functions
summary(Sfun)
plot(Sfun)
pairs(Sfun)

# Get collinearity indexes for each combination of parameters
# collinearity values under 20 are considered acceptable
# you can fit model structures with this data set if collin<20
collin(Sfun)

# what is the set of parameters with the best fit?
Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,5), upper=c(Inf, Inf, 1,1,1))
bestMod=Func(pars=Fit$par)

# what is the probability distribution of the parameters?
MCMC <- modMCMC(f = costFunc, p = Fit$par, niter = 10000, var0 = Fit$var_ms_unweighted, wvar0 = 0.1, updatecov = 50, lower=rep(0,5), upper=c(Inf, Inf, 1,1,1))
summary(MCMC)
plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000)

#############################################################################################
# Can I fit a 3 pool series model?
inipars=c(0.5, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25)
tt=seq(from=0, to=tail(Ct[,1],1), length.out = 500)
Func=function(pars){
  mod=SoilR::ThreepSeriesModel(t=tt,ks=pars[1:3], a21=pars[1]*pars[4], a32=pars[2]*pars[5], C0=100*c(pars[6], pars[7], 1-sum(pars[6:7])), In=0)
  Ct=SoilR::getC(mod)
  return(data.frame(time=tt, Ct=rowSums(Ct)))
}
costFunc=function(pars){
  output=Func(pars)
  return(modCost(model=output, obs=Ct, x="time")) 
}

costFunc(inipars)$model # sums of squared residuals

Sfun <- sensFun(costFunc, inipars)  # sensibility functions
summary(Sfun)
plot(Sfun)
pairs(Sfun)

# Get collinearity indexes for each combination of parameters
# collinearity values under 20 are considered acceptable
# you can fit model structures with this data set if collin<20
collin(Sfun) 

# what is the set of parameters with the best fit?
Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,7), upper=c(Inf, Inf, 1,1,1,1,1))
bestMod=Func(pars=Fit$par)

# what is the probability distribution of the parameters? 
MCMC <- modMCMC(f = costFunc, p = Fit$par, niter = 10000, var0 = Fit$var_ms_unweighted, wvar0 = 0.1, updatecov = 50, lower=rep(0,7), upper=c(Inf, Inf, 1,1,1,1,1))
summary(MCMC)
plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000)