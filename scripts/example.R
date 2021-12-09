# Fig4
# aridec model application example using Berenstecher2021
library(SoilR)
library(FME)
library(aridec)
library(ggplot2)
library(gridExtra)
library(dplyr)

db=loadEntries()
entry=db[["Berenstecher2020"]]
points=data.frame(cbind(entry$timeSeries$Time, entry$timeSeries$PFST))

# Fit a single pool model (Fig 4a)
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
  plot(complete, ylim=c(0,1.2*max(complete[,2])))
  lines(bestMod)
  AIC=(2*length(Fit$par))+2*log(Fit$ms) # formula from SIdb
  AICc=log(Fit$ms)+((n+length(Fit$par))/(n-length(Fit$par)-2)) # for small sample size
  print(paste("AIC = ",AIC))
  SoilRmodel=SoilR::OnepModel(t=tt,k=Fit$par[1], C0=initialCarbon, In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC, AICc=AICc))
}

M1=onepFit(entry$timeSeries[,c(1,2)], initialCarbon = 100) # use columns 1 (Time) and 2 (first variable)
years1=M1$SoilRmodel@times
Ct1=getC(M1$SoilRmodel)
df1=as.data.frame(cbind(years1,Ct1))
graph1=ggplot(df1, aes(years1,Ct1)) +
 geom_line(color="#66C2A5", show.legend = FALSE) +
 ggtitle("(a)") + theme_bw() + ylim(0, 100) +
 theme(axis.title.x=element_blank(),
       axis.text.x=element_blank(),
       plot.title = element_text(face = 'bold')) +
  ylab("Organic matter remaining (%)") +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
 geom_point(data=points, shape=1, aes(x=X1, y=X2))
graph1

# Two pool parallel model using a known value for parameter 3 (lignin % of 9.34; Fig 4b)
twoppFit=function(timeSeries, initialCarbon, inipars=c(1, 0.7)){
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  if(n < 5) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::TwopParallelModel(t=tt,ks=pars[1:2], C0=initialCarbon*c(1-0.0934, 0.0934), In=0, gam=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,2), upper=c(Inf, Inf))
  bestMod=Func(pars=Fit$par)
  print(paste(c("k1=", "k2="),Fit$par))
  plot(complete, ylim=c(0,1.2*max(complete[,2])))
  lines(bestMod)
  AIC=(2*length(Fit$par))+2*log(Fit$ms) # formula from SIdb
  AICc=log(Fit$ms)+((n+length(Fit$par))/(n-length(Fit$par)-2)) # for small sample size
  print(paste("AIC = ",AIC))
  SoilRmodel=SoilR::TwopParallelModel(t=tt,ks=Fit$par[1:2], C0=initialCarbon*c(1-0.0934, 0.0934), In=0, gam=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC, AICc=AICc))
}


M2=twoppFit(entry$timeSeries[,c(1,2)], initialCarbon = 100) 
years2=as.vector(cbind(M2$SoilRmodel@times, M2$SoilRmodel@times, M2$SoilRmodel@times))
Ct2=as.data.frame(getC(M2$SoilRmodel))
total2=rowSums(Ct2)
Ct2=as.vector(cbind(Ct2$V1,Ct2$V2,total2))
df2=as.data.frame(cbind(years2,Ct2,pool=rep(c("fast","slow", "total"), each = 500)))
i <- c(1, 2) 
df2[ , i] <- apply(df2[ , i], 2, function(x) as.numeric(as.character(x)))
graph2=ggplot(df2, aes(years2,Ct2, color=pool)) +
  geom_line() +
  ggtitle("(b)") + theme_bw() + ylim(0, 100) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        plot.title = element_text(face = 'bold')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  ylab("Organic matter remaining (%)") +
  scale_color_manual(values=c("#8DA0CB","#FC8D62","#66C2A5"))#+
  geom_point(data=points, shape=1, aes(x=X1, y=X2))
graph2

# Fit a two pool series model with fixed value of initial proportion of C in pool 2 (Fig 4c)
twopsFit=function(timeSeries, initialCarbon, inipars=c(1, 0.6, 0.6)){
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  if(n < 5) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::TwopSeriesModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], C0=initialCarbon*c(1-0.0934, 0.0934), In=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,3), upper=c(Inf, Inf, 1))
  bestMod=Func(pars=Fit$par)
  print(paste(c("k1=", "k2=", "a21="),Fit$par))
  plot(complete, ylim=c(0,1.2*max(complete[,2])))
  lines(bestMod)
  AIC=(2*length(Fit$par))+2*log(Fit$ms) # formula from SIdb
  AICc=log(Fit$ms)+((n+length(Fit$par))/(n-length(Fit$par)-2)) # for small sample size
  print(paste("AIC = ",AIC))
  SoilRmodel=SoilR::TwopSeriesModel(t=tt,ks=Fit$par[1:2], a21=Fit$par[1]*Fit$par[3], C0=initialCarbon*c(1-0.0934, 0.0934), In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC, AICc=AICc))
}

M3=twopsFit(entry$timeSeries[,c(1,2)], initialCarbon = 100)
years3=as.vector(cbind(M3$SoilRmodel@times, M3$SoilRmodel@times, M3$SoilRmodel@times))
Ct3=as.data.frame(getC(M3$SoilRmodel))
total3=rowSums(Ct3)
Ct3=as.vector(cbind(Ct3$V1,Ct3$V2,total3))
df3=as.data.frame(cbind(years3,Ct3,pool=rep(c("fast","slow", "total"), each = 500)))
i <- c(1, 2) 
df3[ , i] <- apply(df3[ , i], 2, function(x) as.numeric(as.character(x)))
graph3=ggplot(df3, aes(years3,Ct3, color=pool)) +
  geom_line(position = position_dodge(0.1)) +
  ggtitle("(c)") + theme_bw() + 
  theme(plot.title = element_text(face = 'bold')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  xlab("Years") + ylab("Organic matter remaining (%)") +
  scale_color_manual(values=c("#8DA0CB","#FC8D62","#66C2A5"))#+
  geom_point(data=points, shape=1, aes(x=X1, y=X2))
graph3

windows()
grid.arrange(graph1, graph2, graph3, ncol=1)

# to get C release dynamics from the models
berens1=OnepModel(t=seq(from=0, to=tail(points[,1],1), length.out = 500), k=-0.2, C0=100, In=0)
getC(berens1)

berens2=TwopParallelModel(t=seq(from=0, to=tail(points[,1],1), length.out = 500), ks=c(0.24, 1.47e-11), C0=100*c(1-0.0934, 0.0934), In=0, gam=0)
getC(berens2)

berens3=TwopSeriesModel(t=seq(from=0, to=tail(points[,1],1), length.out = 500), ks=c(3.649, 0.244), a21=3.649*1, C0=100*c(1-0.0934, 0.0934), In=0)
getC(berens3)

##########################################################################################
# Andersen2006
entry=db[["CepedaPizarro1990"]]
points=data.frame(cbind(entry$timeSeries$Time, entry$timeSeries$BSBm))
points=na.omit(mutate(points, X1=X1/365)) # days to years

# One pool model

M1=onepFit(points, initialCarbon = 100) # use columns 1 (Time) and 8 (BSBm)
years1=M1$SoilRmodel@times
Ct1=getC(M1$SoilRmodel) # fitted C loss over time
df1=as.data.frame(cbind(years1,Ct1))
graph1=ggplot(df1, aes(years1,Ct1)) +
  geom_line(color="#66C2A5", show.legend = FALSE) +
  ggtitle("(d)") + theme_bw() + ylim(0, 100) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        plot.title = element_text(face = 'bold')) +
  ylab("Organic matter remaining (%)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_point(data=points, shape=1, aes(x=X1, y=X2))+
  scale_x_continuous(breaks=c(0.0,0.25, 0.5, 0.75, 1.0, 1.25))
graph1

# Two pool parallel model using a known value for parameter 3 (lignin % of 6.9; Fig 4b)
twoppFit=function(timeSeries, initialCarbon, inipars=c(1, 0.1)){ # changing inipar no. 2 from 0.7 to 0.1 makes no difference
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  if(n < 5) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::TwopParallelModel(t=tt,ks=pars[1:2], C0=initialCarbon*c(1-0.069, 0.069), In=0, gam=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,2), upper=c(Inf, Inf))
  bestMod=Func(pars=Fit$par)
  print(paste(c("k1=", "k2="),Fit$par))
  plot(complete, ylim=c(0,1.2*max(complete[,2])))
  lines(bestMod)
  AIC=(2*length(Fit$par))+2*log(Fit$ms) # formula from SIdb
  AICc=log(Fit$ms)+((n+length(Fit$par))/(n-length(Fit$par)-2)) # for small sample size
  print(paste("AIC = ",AIC))
  SoilRmodel=SoilR::TwopParallelModel(t=tt,ks=Fit$par[1:2], C0=initialCarbon*c(1-0.069, 0.069), In=0, gam=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC, AICc=AICc))
}

M2=twoppFit(points, initialCarbon = 100) 
years2=as.vector(cbind(M2$SoilRmodel@times, M2$SoilRmodel@times, M2$SoilRmodel@times))
Ct2=as.data.frame(getC(M2$SoilRmodel))
total2=rowSums(Ct2)
Ct2=as.vector(cbind(Ct2$V1,Ct2$V2,total2))
df2=as.data.frame(cbind(years2,Ct2,pool=rep(c("fast","slow", "total"), each = 500)))
i <- c(1, 2) 
df2[ , i] <- apply(df2[ , i], 2, function(x) as.numeric(as.character(x)))
graph2=ggplot(df2, aes(years2,Ct2, color=pool)) +
  geom_line() +
  ggtitle("(e)") + theme_bw() + ylim(0, 100) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        plot.title = element_text(face = 'bold')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  ylab("Organic matter remaining (%)") +
  scale_color_manual(values=c("#8DA0CB","#FC8D62","#66C2A5"))+
  scale_x_continuous(breaks=c(0.0,0.25, 0.5, 0.75, 1.0, 1.25))#+
geom_point(data=points, shape=1, aes(x=X1, y=X2))
graph2

# Fit a two pool series model with fixed value of initial proportion of C in pool 2 (Fig 4c)

twopsFit=function(timeSeries, initialCarbon, inipars=c(1, 0.6, 0.6)){
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  if(n < 5) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::TwopSeriesModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], C0=initialCarbon*c(1-0.069, 0.069), In=0)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }
  
  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,3), upper=c(Inf, Inf, 1))
  bestMod=Func(pars=Fit$par)
  print(paste(c("k1=", "k2=", "a21="),Fit$par))
  plot(complete, ylim=c(0,1.2*max(complete[,2])))
  lines(bestMod)
  AIC=(2*length(Fit$par))+2*log(Fit$ms) # formula from SIdb
  AICc=log(Fit$ms)+((n+length(Fit$par))/(n-length(Fit$par)-2)) # for small sample size
  print(paste("AIC = ",AIC))
  SoilRmodel=SoilR::TwopSeriesModel(t=tt,ks=Fit$par[1:2], a21=Fit$par[1]*Fit$par[3], C0=initialCarbon*c(1-0.069, 0.069), In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC, AICc=AICc))
}

M3=twopsFit(points, initialCarbon = 100)
years3=as.vector(cbind(M3$SoilRmodel@times, M3$SoilRmodel@times, M3$SoilRmodel@times))
Ct3=as.data.frame(getC(M3$SoilRmodel))
total3=rowSums(Ct3)
Ct3=as.vector(cbind(Ct3$V1,Ct3$V2,total3))
df3=as.data.frame(cbind(years3,Ct3,pool=rep(c("fast","slow", "total"), each = 500)))
i <- c(1, 2) 
df3[ , i] <- apply(df3[ , i], 2, function(x) as.numeric(as.character(x)))
graph3=ggplot(df3, aes(years3,Ct3, color=pool)) +
  geom_line() +#position = position_dodge(0.1)) +
  ggtitle("(f)") + theme_bw() + 
  theme(plot.title = element_text(face = 'bold')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  xlab("Years") + ylab("Organic matter remaining (%)") +
  scale_color_manual(values=c("#8DA0CB","#FC8D62","#66C2A5")) +
  scale_x_continuous(breaks=c(0.0,0.25, 0.5, 0.75, 1.0, 1.25))
  #+
  geom_point(data=points, shape=1, aes(x=X1, y=X2))
graph3

windows()
grid.arrange(graph1, graph2, graph3, ncol=1)
