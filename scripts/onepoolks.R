# load libraries
library(SoilR)
library(FME)
library(aridec)
library(dplyr)

# load entire database
db=loadEntries("/Users/agustin/Documents/GitHub/aridec/data/")

# load One pool function for a single time series
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
  
  inipars=c(-1) # i.e.: -1*C0/complete[1,2] = -100/100 = -1
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower= -Inf, upper=0)
  bestMod=Func(pars=Fit$par)
 # plot(complete, ylim=c(0,1.2*max(complete[,2])))
 # lines(bestMod)
  summary(Fit)
  SoilRmodel=SoilR::OnepModel(t=tt,k=Fit$par[1], C0=C0, In=0)
  return(list(FME=Fit, k=Fit$par, summary=summary(Fit)))
}

# load loop function for an entire entry
onePentry= function(Ct, time) {lapply(Ct, onepFit, time)}

# load a single entry
entry=db[["Strojan1987"]]
df=entry$timeSeries

######################################################################
# prueba con strojan para ver que pasa cuando cambio la escala temporal
df=mutate(df, Time=Time*7) # weeks to days
df=mutate(df, Time=Time/4) # weeks to months
df=mutate(df, Time=Time/48) # weeks to years
df=mutate(df, Time=Time/30) # days to months

strojankm=data.frame(unlist(onePentry(Ct=df[-1], time = df[,1])))
strojanfull=onePentry(Ct=df[-1], time = df[,1])

plot=data.frame(timeScale=c("days", "weeks", "months", "years"), k=c(ksdays$Larrea$k, ksweeks$Larrea$k, ksmonths$Larrea$k, ksyears$Larrea$k))
perday=c(plot[1,2], plot[2,2]/7, plot[3,2]/30, plot[4,2]/365)
plot=cbind(plot, perday)

plot(complete, ylim=c(0,1.2*max(complete[,2])))
lines(bestMod)
############################################################### eliminar luego
if (entry$variables$V1$units == "years") { 
  Ct=df[-1] 
  time=df[,1] 
  ks=onePentry(Ct, time) 
} else if (entry$variables$V1$units == "months") { 
  df=mutate(df, Time=Time/12) 
  Ct=df[-1] 
  time=df[,1] 
  ks=onePentry(Ct, time)
} else if (entry$variables$V1$units == "weeks") { 
  df=mutate(df, Time=Time/48) 
  Ct=df[-1] 
  time=df[,1] 
  ks=onePentry(Ct, time)
} else { 
  df=mutate(df, Time=Time/365) 
  Ct=df[-1] 
  time=df[,1] 
  ks=onePentry(Ct, time) 
}

write.csv(ks,"~/kxyears/k184.csv") # save k values for the entry as .csv file
