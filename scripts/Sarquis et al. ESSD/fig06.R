library(SoilR)
library(FME)
library(aridec)
library(dplyr)
library(plyr)
library(cowplot)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)

# First, generate the dataframe with collinearity results for the entire database

db=loadEntries("~/aridec/data/") 

inipars2pp=c(1, 0.5, 0.5)
Func2pp=function(pars){
  mod=SoilR::TwopParallelModel(t=tt,ks=pars[1:2], C0=100*c(pars[3], 1-pars[3]), In=0, gam=0)
  Ct=SoilR::getC(mod)
  return(data.frame(time=tt, Ct=rowSums(Ct)))
}
costFunc2pp=function(pars){
  output=Func2pp(pars)
  return(modCost(model=output, obs=Ct, x="time")) 
}

inipars2ps=c(0.5, 0.25, 0.25, 0.15)
Func2ps=function(pars){
  mod=SoilR::TwopSeriesModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], C0=100*c(pars[4], 1-pars[4]), In=0)
  Ct=SoilR::getC(mod)
  return(data.frame(time=tt, Ct=rowSums(Ct)))
}
costFunc2ps=function(pars){
  output=Func2ps(pars)
  return(modCost(model=output, obs=Ct, x="time")) 
}

sink(file='~/Documents/collinearity.csv')

for (i in 1:184) {
  
  entry=db[[i]]
  Ct=entry$timeSeries[,c(1,2)] 
  colnames(Ct)=c("time", "Ct")
  
  if (entry$variables$V1$units == "years") { 
    Ct=mutate(Ct, time=time*12)
  } else if (entry$variables$V1$units == "weeks") { 
    Ct=mutate(Ct, time=time/4.3)
  } else if (entry$variables$V1$units == "days") { 
    Ct=mutate(Ct, time=time/30)
  } 
  
  # Can I fit a 2 pool parallel model?
  tt=seq(from=0, to=tail(Ct[,1],1), length.out = 500)
  
  print(entry$citationKey)
  costFunc2pp(inipars2pp)$model
  print(collin(sensFun(costFunc2pp, inipars2pp)))
  
  # Can I fit a 2 pool series model?
  print(entry$citationKey)
  costFunc2ps(inipars2ps)$model
  print(collin(sensFun(costFunc2ps, inipars2ps)))
  
}

sink() 

############################################################################################
# Now get the number of points in time per data entry and add it to the previously generated dataframe

timepoints=lapply(db, function(x){length(na.omit(x$timeSeries[,2]))})

timepoints=t(as.data.frame(totaltimepoints))

rep.row<-function(x,n){
  data.frame(rep(x,each=n))
}

reptimepoints=rep.row(timepoints,15)

write.csv(reptimepoints,'~/Documents/timepointsrep.csv') # add this as a column to the dataframe

###############################################################################################

df=read.csv("") # load the dataframe
summary(df)
df$X1=as.factor(df$X1)
df$X2=as.factor(df$X2)
df$X3=as.factor(df$X3)
df$X4=as.factor(df$X4)
df$X5=as.factor(df$X5)
df$X6=as.factor(df$X6)
df$X7=as.factor(df$X7)
df$N=as.factor(df$N)
df$length=as.factor(df$Timepoints)

myColors <- brewer.pal(6,"Set2")
names(myColors) <- levels(df$Timepoints)
colScale <- scale_colour_manual(name = "Timepoints",values = myColors)

gtwopp=ggplot(filter(df, Model == "2pp" & X1!= "0" & X2!= "0"), aes(N, log.col, color=Timepoints))
plot1=gtwopp+geom_point(position="jitter", size=3, show.legend = FALSE) +
  geom_abline(intercept= log10(20), slope=0) + ylim(0,NA) + theme_bw() +
  labs(y="log gamma", x=("Parameter combination")) +
  ggtitle("Two Pool Parallel Model") +  geom_vline(xintercept = c("2", "3"), linetype=2) +
  scale_colour_manual(colScale)

gtwops=ggplot(filter(df, Model == "2ps" & X1!= "0" & X2!= "0" & X3!= "0"), aes(N, log.col, color=Timepoints))
plot2=gtwops+geom_point(position="jitter", size=3, show.legend = FALSE) +
  geom_abline(intercept= log10(20), slope=0) + 
  theme_bw() +
  labs(y="log gamma", x=("Number of parameters in combination")) + 
  ggtitle("Two Pool Series Model") +  geom_vline(xintercept = c("3", "4"), linetype=2) +
  scale_colour_manual(colScale)

