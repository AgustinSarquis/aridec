# Fig7
# aridec model application example using Day2018
library(SoilR)
library(FME)
library(aridec)
library(ggplot2)
library(gridExtra)
library(dplyr)

db=loadEntries('~/aridec/data/')
entry=db[["Day2018"]]
points=data.frame(cbind(entry$timeSeries$Time, entry$timeSeries$Sc_fs))
points=mutate(points, X1=X1/30)

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
  AIC=(2*(length(Fit$par)+1))+2*log(Fit$ms) # formula from SIdb
  AICc=log(Fit$ms)+((n+length(Fit$par)+1)/(n-length(Fit$par)-2)) # for small sample size
  print(paste("AICc = ",AICc))
  SoilRmodel=SoilR::OnepModel(t=tt,k=Fit$par[1], C0=initialCarbon, In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC, AICc=AICc))
}

M1=onepFit(points, initialCarbon = 100) 
months1=M1$SoilRmodel@times
Ct1=getC(M1$SoilRmodel)
df1=as.data.frame(cbind(months1,Ct1))
graph1=ggplot(df1, aes(months1,Ct1)) +
  geom_line(color="#66C2A5", show.legend = FALSE) +
  ggtitle("(a)") + theme_bw() + ylim(0, 100) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        plot.title = element_text(face = 'bold')) +
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
    mod=SoilR::TwopParallelModel(t=tt,ks=pars[1:2], C0=initialCarbon*c(1-0.09, 0.09), In=0, gam=0)
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
  AIC=(2*(length(Fit$par)+1))+2*log(Fit$ms) # formula from SIdb
  AICc=log(Fit$ms)+((n+length(Fit$par)+1)/(n-length(Fit$par)-2)) # for small sample size
  print(paste("AICc = ",AICc))
  SoilRmodel=SoilR::TwopParallelModel(t=tt,ks=Fit$par[1:2], C0=initialCarbon*c(1-0.09, 0.09), In=0, gam=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC, AICc=AICc))
}


M2=twoppFit(points, initialCarbon = 100) 
months2=as.vector(cbind(M2$SoilRmodel@times, M2$SoilRmodel@times, M2$SoilRmodel@times))
Ct2=as.data.frame(getC(M2$SoilRmodel))
total2=rowSums(Ct2)
Ct2=as.vector(cbind(Ct2$V1,Ct2$V2,total2))
df2=as.data.frame(cbind(months2,Ct2,pool=rep(c("fast","slow", "total"), each = 500)))
i <- c(1, 2) 
df2[ , i] <- apply(df2[ , i], 2, function(x) as.numeric(as.character(x)))
graph2=ggplot(df2, aes(months2,Ct2, color=pool)) +
  geom_line() +
  ggtitle("(b)") + theme_bw() + ylim(0, 100) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        plot.title = element_text(face = 'bold')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  ylab("Organic matter remaining (%)") +
  scale_color_manual(values=c("#8DA0CB","#FC8D62","#66C2A5"))+
  geom_point(data=points, shape=1, aes(x=X1, y=X2), inherit.aes = FALSE )
graph2

# Fit a two pool series model with fixed value of initial proportion of C in pool 2 (Fig 4c)
twopsFit=function(timeSeries, initialCarbon, inipars=c(0.02, 0.015, 0.6)){
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  if(n < 5) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)
  
  Func=function(pars){
    mod=SoilR::TwopSeriesModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], C0=initialCarbon*c(1-0.09, 0.09), In=0)
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
  AIC=(2*(length(Fit$par)+1))+2*log(Fit$ms) # formula from SIdb
  AICc=log(Fit$ms)+((n+length(Fit$par)+1)/(n-length(Fit$par)-2)) # for small sample size
  print(paste("AICc = ",AICc))
  SoilRmodel=SoilR::TwopSeriesModel(t=tt,ks=Fit$par[1:2], a21=Fit$par[1]*Fit$par[3], C0=initialCarbon*c(1-0.09, 0.09), In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC, AICc=AICc))
}

M3=twopsFit(points, initialCarbon = 100)
months3=as.vector(cbind(M3$SoilRmodel@times, M3$SoilRmodel@times, M3$SoilRmodel@times))
Ct3=as.data.frame(getC(M3$SoilRmodel))
total3=rowSums(Ct3)
Ct3=as.vector(cbind(Ct3$V1,Ct3$V2,total3))
df3=as.data.frame(cbind(months3,Ct3,pool=rep(c("fast","slow", "total"), each = 500)))
i <- c(1, 2) 
df3[ , i] <- apply(df3[ , i], 2, function(x) as.numeric(as.character(x)))
graph3=ggplot(df3, aes(months3,Ct3, color=pool)) +
  geom_line(position = position_dodge(0.1)) +
  ggtitle("(c)") + theme_bw() + 
  theme(plot.title = element_text(face = 'bold')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  xlab("Months") +
  scale_color_manual(values=c("#8DA0CB","#FC8D62","#66C2A5"))+
  geom_point(data=points, shape=1, aes(x=X1, y=X2),inherit.aes = FALSE)
graph3

pdf('~/Documents/fig07.pdf')
grid.arrange(graph1, graph2, graph3, ncol=1)
dev.off()
