# load libraries
library(SoilR)
library(FME)
library(aridec)
library(dplyr)

# load single entry
db=loadEntries("~/aridec/data/")
entry=db$Henry2008
Ct=entry$timeSeries[,c(1,2)] 
colnames(Ct)=c("time", "Ct")

Ct=mutate(Ct, time=time/30) # transform days to months
Ct=mutate(Ct, time=time*12) # transform years to months
Ct=mutate(Ct, time=time/4) # transform weeks to months
Ct=mutate(Ct, time=time/12) # transfrom months to years

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

###########################################################################################
# graficos
library(ggplot2)
library(cowplot)
library(gridExtra)
library(RColorBrewer)
df=read.csv("C:/Users/musi/OneDrive - Facultad de Agronom?a - Universidad de Buenos Aires/aridec manuscript/coldf.csv")
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
df$length.1=as.factor(df$length.1)

myColors <- brewer.pal(6,"Set2")
names(myColors) <- levels(df$length.1)
colScale <- scale_colour_manual(name = "length.1",values = myColors)

gtwopp=ggplot(subset(df, model %in% "twopp"), aes(N, log.col, color=length.1))
plot1=gtwopp+geom_point(position="jitter", size=3, show.legend = FALSE) +
 geom_abline(intercept= log10(20), slope=0) + ylim(0,NA) + theme_bw() +
  labs(y="log ??", x=("Number of parameters in combination")) +
  ggtitle("Two Pool Parallel Model") +  geom_vline(xintercept = c("2", "3"), linetype=2) +
  scale_colour_manual(name= "length.1", values=myColors)

gtwops=ggplot(subset(df, model %in% "twops"), aes(N, log.col, color=length.1))
plot2=gtwops+geom_point(position="jitter", size=3, show.legend = FALSE) +
  geom_abline(intercept= log10(20), slope=0) + #ylim(0,3.5) + 
  theme_bw() +
  labs(y="log ??", x=("Number of parameters in combination")) + 
  ggtitle("Two Pool Series Model") +  geom_vline(xintercept = c("2", "3", "4"), linetype=2) +
  scale_colour_manual(name= "length.1", values=myColors)

gtwopf=ggplot(subset(df, model %in% "twopf"), aes(N, log.col, color=length.1))
plot3=gtwopf+geom_point(position="jitter", size=3, show.legend = FALSE) +
  geom_abline(intercept= log10(20), slope=0) + ylim(0,8) + 
  theme_bw() + labs(y="log ??", x=("Number of parameters in combination")) +
  ggtitle("Two Pool Model with Feedback") +  geom_vline(xintercept = c("2", "3", "4","5"), linetype=2) +
  scale_colour_manual(name= "length.1", values=myColors)

gthreepp=ggplot(subset(df, model %in% "threepp"), aes(N, log.col, color=length.1))
plot4=gthreepp+geom_point(position="jitter", size=3, show.legend = FALSE) +
  geom_abline(intercept= log10(20), slope=0) +# ylim(0,10) + 
  theme_bw() +  labs(y="log ??", x=("Number of parameters in combination")) + 
  ggtitle("Three Pool Parallel Model") +  geom_vline(xintercept = c("2", "3", "4","5"), linetype=2) +
  scale_colour_manual(name= "length.1", values=myColors)

gthreeps=ggplot(subset(df, model %in% "threeps"), aes(N, log.col, color=length.1))
plot5=gthreeps+geom_point(position="jitter", size=3, show.legend = FALSE) +
  geom_abline(intercept= log10(20), slope=0) +# ylim(0,10) + 
  theme_bw() +  labs(y="log ??", x=("Number of parameters in combination")) + 
  ggtitle("Three Pool Series Model") +  geom_vline(xintercept = c("2", "3", "4","5","6","7"), linetype=2) +
  scale_colour_manual(name= "length.1", values=myColors)

# I only made this plot to extract the legend
plotL=gtwopp+geom_point(position="jitter", size=3) +
  geom_abline(intercept= log10(20), slope=0) + ylim(0,NA) + theme_bw() +
  labs(y="log ??", x=("Number of parameters in combination")) +
  ggtitle("Two Pool Parallel Model") +  geom_vline(xintercept = c("2", "3"), linetype=2) +
  scale_colour_manual(name= "Number of time points", values=myColors) +
  theme(legend.position="bottom") 

legend <- cowplot::get_legend(plotL)
windows()
grid.arrange(plot1, plot2, plot3, plot4, plot5, legend, ncol=2)

###################################################################################################################
# porcentaje de entradas identificables por modelo
twoppsub = subset(df, model == "twopp", select = c("X1", "X2", "X3", "N", "log.col")) 
all_pars_2pp= filter(twoppsub, N %in% "3")
identifiable=c(all_pars_2pp$log.col <= log10(20))
all_pars_2pp=cbind(all_pars_2pp,identifiable)
length(all_pars_2pp$identifiable[all_pars_2pp$identifiable==T])
18/29
rest_2pp=filter(twoppsub, N %in% "2")
identifiable=c(rest_2pp$log.col <= log10(20))
rest_2pp=cbind(rest_2pp,identifiable)
length(rest_2pp$identifiable[rest_2pp$identifiable==T])
86/90
ks_2pp=filter(twoppsub, X1 %in% "1", X2 %in% "1", X3 %in% "0")
identifiable=c(ks_2pp$log.col <= log10(20))
ks_2pp=cbind(ks_2pp,identifiable)
length(ks_2pp$identifiable[ks_2pp$identifiable==T])

twopssub = subset(df, model == "twops", select = c("X1", "X2", "X3", "X4","N", "log.col")) 
all_2ps=filter(twopssub, N%in% "4")
identifiable=c(all_2ps$log.col <= log10(20))
all_2ps=cbind(all_2ps,identifiable)
length(all_2ps$identifiable[all_2ps$identifiable==T])
12/25
rest2_2ps=filter(twopssub, N%in% "2")
identifiable=c(rest2_2ps$log.col <= log10(20))
rest2_2ps=cbind(rest2_2ps,identifiable)
length(rest2_2ps$identifiable[rest2_2ps$identifiable==T])
142/150
rest3_2ps=filter(twopssub, N%in% "3")
identifiable=c(rest3_2ps$log.col <= log10(20))
rest3_2ps=cbind(rest3_2ps,identifiable)
length(rest3_2ps$identifiable[rest3_2ps$identifiable==T])
69/101
ks_a_2ps=filter(twopssub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "0")
identifiable=c(ks_a_2ps$log.col <= log10(20))
ks_a_2ps=cbind(ks_a_2ps,identifiable)
length(ks_a_2ps$identifiable[ks_a_2ps$identifiable==T])
22/26

twopfsub = subset (df, model == "twopf", select = c("X1", "X2", "X3", "X4", "X5","N", "log.col"))
all_2pf=filter(twopfsub, N%in% "5")
all_2pf=cbind(all_2pf,identifiable=c(all_2pf$log.col <= log10(20)))
length(all_2pf$identifiable[all_2pf$identifiable==T])
rest2_2pf=filter(twopfsub, N%in% "2")
identifiable=c(rest2_2pf$log.col <= log10(20))
rest2_2pf=cbind(rest2_2pf,identifiable)
length(rest2_2pf$identifiable[rest2_2pf$identifiable==T])
86/150
rest3_2pf=filter(twopfsub, N%in% "3")
identifiable=c(rest3_2pf$log.col <= log10(20))
rest3_2pf=cbind(rest3_2pf,identifiable)
length(rest3_2pf$identifiable[rest3_2pf$identifiable==T])
7/150
rest4_2pf=filter(twopfsub, N%in% "4")
identifiable=c(rest4_2pf$log.col <= log10(20))
rest4_2pf=cbind(rest4_2pf,identifiable)
length(rest4_2pf$identifiable[rest4_2pf$identifiable==T])
0
ks_as_2pf=filter(twopfsub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "1", X5 %in% "0")
ks_as_2pf=cbind(ks_as_2pf,identifiable=c(ks_as_2pf$log.col <= log10(20)))

threeppsub =subset(df, model == "threepp", select = c("X1", "X2", "X3", "X4", "X5","N", "log.col"))
rest4_3pp=filter(threeppsub, N%in% "4")
identifiable=c(rest4_3pp$log.col <= log10(20))
rest4_3pp=cbind(rest4_3pp,identifiable)
length(rest4_3pp$identifiable[rest4_3pp$identifiable==T])
26/75
rest3_3pp=filter(threeppsub, N%in% "3")
identifiable=c(rest3_3pp$log.col <= log10(20))
rest3_3pp=cbind(rest3_3pp,identifiable)
length(rest3_3pp$identifiable[rest3_3pp$identifiable==T])
102/150
rest2_3pp=filter(threeppsub, N%in% "2")
identifiable=c(rest2_3pp$log.col <= log10(20))
rest2_3pp=cbind(rest2_3pp,identifiable)
length(rest2_3pp$identifiable[rest2_3pp$identifiable==T])
134/150
ks_c1_3pp=filter(threeppsub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "1", X5 %in% "0")
ks_c1_3pp=cbind(ks_c1_3pp, identifiable=c(ks_c1_3pp$log.col <= log10(20)))
ks_c2_3pp=filter(threeppsub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "0", X5 %in% "1")
ks_c2_3pp=cbind(ks_c2_3pp, identifiable=c(ks_c2_3pp$log.col <= log10(20)))
ks_3pp=filter(threeppsub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "0", X5 %in% "0")
ks_3pp=cbind(ks_3pp, identifiable=c(ks_3pp$log.col <= log10(20)))
count(ks_3pp$identifiable)

threepssub= subset(df, model =="threeps", select= c("X1", "X2", "X3", "X4", "X5","X6","X7","N", "log.col"))
rest6_3ps=filter(threepssub, N%in% "6")
identifiable=c(rest6_3ps$log.col <= log10(20))
rest6_3ps=cbind(rest6_3ps,identifiable)
length(rest6_3ps$identifiable[rest6_3ps$identifiable==T])
2/35
rest5_3ps=filter(threepssub, N%in% "5")
identifiable=c(rest5_3ps$log.col <= log10(20))
rest5_3ps=cbind(rest5_3ps,identifiable)
length(rest5_3ps$identifiable[rest5_3ps$identifiable==T])
11/105
rest4_3ps=filter(threepssub, N%in% "4")
identifiable=c(rest4_3ps$log.col <= log10(20))
rest4_3ps=cbind(rest4_3ps,identifiable)
length(rest4_3ps$identifiable[rest4_3ps$identifiable==T])
49/175
rest3_3ps=filter(threepssub, N%in% "3")
identifiable=c(rest3_3ps$log.col <= log10(20))
rest3_3ps=cbind(rest3_3ps,identifiable)
length(rest3_3ps$identifiable[rest3_3ps$identifiable==T])
96/175
rest2_3ps=filter(threepssub, N%in% "2")
identifiable=c(rest2_3ps$log.col <= log10(20))
rest2_3ps=cbind(rest2_3ps,identifiable)
length(rest2_3ps$identifiable[rest2_3ps$identifiable==T])
89/105
ks_as_c1_3ps=filter(threepssub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "1", X5 %in% "1", X6 %in% "1", X7 %in% "0")
ks_as_c1_3ps=cbind(ks_as_c1_3ps,identifiable=c(ks_as_c1_3ps$log.col <= log10(20)))
ks_as_c2_3ps=filter(threepssub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "1", X5 %in% "1", X6 %in% "0", X7 %in% "1")
ks_as_c2_3ps=cbind(ks_as_c2_3ps,identifiable=c(ks_as_c2_3ps$log.col <= log10(20)))
ks_as_3ps=filter(threepssub, X1 %in% "1", X2 %in% "1", X3 %in% "1", X4 %in% "1", X5 %in% "1", X6 %in% "0", X7 %in% "0")
ks_as_3ps=cbind(ks_as_3ps,identifiable=c(ks_as_3ps$log.col <= log10(20)))
