library(aridec)
library(dplyr)

db=loadEntries("/Users/agustin/Documents/GitHub/aridec/data/")

entry=db[[70]]
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

# funcion para calcular k para una serie de tiempo
nonlineark=function(Mt, t){
    nls=nls(Mt ~ 100*exp(-k*t), start = list(k = 0.1), na.action=na.omit)
    return(coef=summary(nls)$coefficients)
    }

kxentry= function(Mt, t) {sapply(Mt, nonlineark, t)}

Hewins2017k=as.data.frame(kxentry(Mt, t))
Hewins2017k=as.data.frame(t(Hewins2017k))
colnames(Hewins2017k)= c("k", "SE", "t value", "p value")
write.csv(Hewins2017k,"~/Documents/nonlineark/Hewins2017.csv")
