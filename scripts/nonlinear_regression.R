library(aridec)
library(dplyr)

db=loadEntries("~/aridec/data/")

# funcion para calcular k para una serie de tiempo
nonlineark=function(Mt, t){
  nls=nls(Mt ~ 100*exp(-k*t), start = list(k = 0.1), na.action=na.omit)
  return(coef=summary(nls)$coefficients)
}

kxentry= function(Mt, t) {sapply(Mt, nonlineark, t)}

entry=db[[184]]
df=entry$timeSeries
Mt=df[-1]
if (entry$variables$V1$units == "years") { # transform all time units to months
  df=mutate(df, Time=Time*12)
  t=df[,1]
} else if (entry$variables$V1$units == "weeks") { 
  df=mutate(df, Time=Time/4)
  t=df[,1]
} else if (entry$variables$V1$units == "days") { 
  df=mutate(df, Time=Time/30)
  t=df[,1]
} else { 
  t=df[,1]
  }

entry$citationKey

Zheng2010k=as.data.frame(kxentry(Mt, t))
Zheng2010k=as.data.frame(t(Zheng2010k))
colnames(Zheng2010k)= c("k", "SE", "t value", "p value")
write.csv(Zheng2010k,"~/nonlineark/Zheng2010.csv")
