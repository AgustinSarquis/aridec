# load libraries
library(SoilR)
library(FME)
library(aridec)
library(dplyr)

# load single entry
db=loadEntries()
entry=db[[8]]
Ct=entry$timeSeries[,c(1,2)] 
colnames(Ct)=c("time", "Ct")

Ct=mutate(Ct, time=time*30) # transformo meses a dias
Ct=mutate(Ct, time=time*365) # anos a dias

# 2 pool parallel model 
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

costFunc(inipars)$model # suma de cuadrados residuales

Sfun <- sensFun(costFunc, inipars)  # funciones de sensibilidad
summary(Sfun)
plot(Sfun)
pairs(Sfun)

# cuales parametros pueden ser estimados simultaneamente?
collin(Sfun) # colinearidad

# cual es el set de parametros optimos con el mejor ajuste?
Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,7), upper=c(Inf, Inf, 1,1,1,1,1))
bestMod=Func(pars=Fit$par)

# cual es la distribucion de probabilidad de los parametros? 
MCMC <- modMCMC(f = costFunc, p = Fit$par, niter = 10000, var0 = Fit$var_ms_unweighted, wvar0 = 0.1, updatecov = 50, lower=rep(0,7), upper=c(Inf, Inf, 1,1,1,1,1))
summary(MCMC)
plot(MCMC, Full = TRUE)
pairs(MCMC, nsample = 1000)

###########################################################################################
# graficos
library(ggplot2)
df=read.csv("~/coldf.csv")
summary(df)
df$X1=as.factor(df$X1)
df$X2=as.factor(df$X2)
df$X3=as.factor(df$X3)
df$X4=as.factor(df$X4)
df$X5=as.factor(df$X5)
df$X6=as.factor(df$X6)
df$X7=as.factor(df$X7)
df$N=as.factor(df$N)
df$length=as.factor(df$length)
# modifique el dataframe para agrupar valores >8
gtwopp=ggplot(subset(df, model %in% "twopp"), aes(N, log.col, color=length))
gtwopp+geom_point(position="jitter", size=3) +
 geom_abline(intercept= log10(20), slope=0) + ylim(0,3.5) + theme_bw() +
  labs(y="log(Collinearity)", x=("Number of parameters in combination")) + 
  labs(colour="Number of time points") +   ggtitle("Two Pool Parallel Model") +
  geom_vline(xintercept = c("2", "3"), linetype=2)

gtwops=ggplot(subset(df, model %in% "twops"), aes(N, log.col, color=length))
gtwops+geom_point(position="jitter", size=3) +
  geom_abline(intercept= log10(20), slope=0) + #ylim(0,3.5) + 
  theme_bw() +
  labs(y="log(Collinearity)", x=("Number of parameters in combination")) + 
  labs(colour="Number of time points") +   ggtitle("Two Pool Series Model") +
  geom_vline(xintercept = c("2", "3", "4"), linetype=2)

gtwopf=ggplot(subset(df, model %in% "twopf"), aes(N, log.col, color=length))
gtwopf+geom_point(position="jitter", size=3) +
  geom_abline(intercept= log10(20), slope=0) + #ylim(0,3.5) + 
  theme_bw() +
  labs(y="log(Collinearity)", x=("Number of parameters in combination")) + 
  labs(colour="Number of time points") +   ggtitle("Two Pool Model with Feedback") +
  geom_vline(xintercept = c("2", "3", "4","5"), linetype=2)

gthreepp=ggplot(subset(df, model %in% "threepp"), aes(N, log.col, color=length))
gthreepp+geom_point(position="jitter", size=3) +
  geom_abline(intercept= log10(20), slope=0) +# ylim(0,10) + 
  theme_bw() +
  labs(y="log(Collinearity)", x=("Number of parameters in combination")) + 
  labs(colour="Number of time points") +   ggtitle("Three Pool Parallel Model") +
  geom_vline(xintercept = c("2", "3", "4","5"), linetype=2)
# use el dataframe sin modiificar
gthreeps=ggplot(subset(df, model %in% "threeps"), aes(N, log.col, color=length))
gthreeps+geom_point(position="jitter", size=3) +
  geom_abline(intercept= log10(20), slope=0) +# ylim(0,10) + 
  theme_bw() +
  labs(y="log(Collinearity)", x=("Number of parameters in combination")) + 
  labs(colour="Number of time points") +   ggtitle("Three Pool Series Model") +
  geom_vline(xintercept = c("2", "3", "4","5","6","7"), linetype=2)

###################################################################################################################
# porcentaje de entradas identificables por modelo
twoppsub = subset(df, model == "twopp", select = c("X1", "X2", "X3", "N", "log.col")) 
all_pars_2pp= filter(twoppsub, N %in% "3")
identifiable=c(all_pars_2pp$log.col <= log10(20))
all_pars_2pp=cbind(all_pars_2pp,identifiable)
count(all_pars_2pp$identifiable)
ks_2pp=filter(twoppsub, X1 %in% "1", X2 %in% "1", X3 %in% "0")
identifiable=c(ks_2pp$log.col <= log10(20))
ks_2pp=cbind(ks_2pp,identifiable)
count(ks_2pp$identifiable)

twopssub = subset(df, model == "twops", select = c("X1", "X2", "X3", "X4","N", "log.col")) 
all_2ps=filter(twopssub, N%in% "4")
identifiable=c(all_2ps$log.col <= log10(20))
all_2ps=cbind(all_2ps,identifiable)
count(all_2ps$identifiable)/20
ks_a_2ps=filter(twopssub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "0")
identifiable=c(ks_a_2ps$log.col <= log10(20))
ks_a_2ps=cbind(ks_a_2ps,identifiable)
count(ks_a_2ps$identifiable)/20

twopfsub = subset (df, model == "twopf", select = c("X1", "X2", "X3", "X4", "X5","N", "log.col"))
all_2pf=filter(twopfsub, N%in% "5")
all_2pf=cbind(all_2pf,identifiable=c(all_2pf$log.col <= log10(20)))
count(all_2pf$identifiable)
ks_as_2pf=filter(twopfsub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "1", X5 %in% "0")
ks_as_2pf=cbind(ks_as_2pf,identifiable=c(ks_as_2pf$log.col <= log10(20)))

threeppsub =subset(df, model == "threepp", select = c("X1", "X2", "X3", "X4", "X5","N", "log.col"))
ks_c1_3pp=filter(threeppsub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "1", X5 %in% "0")
ks_c1_3pp=cbind(ks_c1_3pp, identifiable=c(ks_c1_3pp$log.col <= log10(20)))
count(ks_c1_3pp$identifiable)
ks_c2_3pp=filter(threeppsub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "0", X5 %in% "1")
ks_c2_3pp=cbind(ks_c2_3pp, identifiable=c(ks_c2_3pp$log.col <= log10(20)))
count(ks_c2_3pp$identifiable)
ks_3pp=filter(threeppsub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "0", X5 %in% "0")
ks_3pp=cbind(ks_3pp, identifiable=c(ks_3pp$log.col <= log10(20)))
count(ks_3pp$identifiable)

threepssub= subset(df, model =="threeps", select= c("X1", "X2", "X3", "X4", "X5","X6","X7","N", "log.col"))
ks_as_3ps=filter(threepssub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "1", X5 %in% "1", X6 %in% "0", X7 %in% "0")
ks_as_3ps=cbind(ks_as_3ps,identifiable=c(ks_as_3ps$log.col <= log10(20)))
count(ks_as_3ps$identifiable)                
