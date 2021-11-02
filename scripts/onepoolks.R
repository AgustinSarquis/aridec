# load libraries
library(SoilR)
library(FME)
library(aridec)
library(dplyr)

# load entire database
db=loadEntries()

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
  
  inipars=c(-1*C0/complete[1,2])
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower= -Inf, upper=0)
  bestMod=Func(pars=Fit$par)
  SoilRmodel=SoilR::OnepModel(t=tt,k=Fit$par[1], C0=C0, In=0)
  return(list(k=Fit$par))
}

# load loop function for an entire entry
onePloop= function(Ct, time) {lapply(Ct, onepFit, time)}

# load a single entry
entry=db[[184]]
df=entry$timeSeries
print(entry$variables$V1$units) # check time units (if units=days, leave unchanged; if not...)

df=mutate(df, Time=Time*7) # transform weeks to days
df=mutate(df, Time=Time*30) # transform months to days
df=mutate(df, Time=Time*365) # transform years to days

Ct=df[-1] # create data frame with mass loss data
time=df[,1] # create time vector
write.csv(onePloop(Ct, time),"~/k182.csv") # save k values for the entry as .csv file
