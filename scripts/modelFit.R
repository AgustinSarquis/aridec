library(SoilR)
library(FME)
library(aridec)

db=loadEntries()
deco=db[[x]] # replace x with database entry

# after testing for collinearity procede to fit suitable models

# Run a One Pool Model
M1=onepFit(deco$timeSeries[,c(1,2)], initialCarbon = 100) # use columns 1 (Time) and 2 (first variable)
months=M1$SoilRmodel@times
Ct=getC(M1$SoilRmodel)
matplot(months,Ct, type="l",lty=1, col=1, ylab="Carbon remaining (%)")

SA1=1/FMEmodel$par[1] # in this case mean system age and transit time are equal

# Run a Two Pool Parallel Model
M2=twoppFit(deco$timeSeries[,c(1,2)], initialCarbon = 100) 
months=M2$SoilRmodel@times
Ct=getC(M2$SoilRmodel)
matplot(months,Ct, type="l",lty=1, col=2:3, ylim=c(0,100), ylab="Carbon remaining (%)")
legend("topright", c("Slow","Fast"), lty=1, col=c(2:3), bty="n")
points(deco$timeSeries[,c(1,2)])
lines(months, rowSums(Ct), lwd=2)

k=M2$FMEmodel$par[c(1,2)]
A2=diag(-k)
u2=matrix(c(M2$FMEmodel$par[3], 1-M2$FMEmodel$par[3]), ncol=1)
SA2=systemAge(A=A2, u=u2)
TT2=transitTime(A=A2, u=u2)

# Run a two Pool Series Model
M3=twopsFit(deco$timeSeries[,c(1,2)], initialCarbon = 100)
months=M3$SoilRmodel@times
Ct=getC(M3$SoilRmodel)
matplot(months,Ct, type="l",lty=1, col=c(2,3), ylim=c(0,100), ylab="Carbon remaining (%)")
legend("topright", c("Total","Slow","Fast"), lty=1, col=c(1,3,2), bty="n")
points(deco$timeSeries[,c(1,2)])
lines(months, rowSums(Ct), lwd=2)

k=M3$FMEmodel$par[c(1,2)]
A3=diag(-k)
A3[2,1]=M3$FMEmodel$par[3]
u3=matrix(c(M3$FMEmodel$par[4], 1-M3$FMEmodel$par[4]), ncol=1)
SA3=systemAge(A=A3, u=u3)
TT3=transitTime(A=A3, u=u3)

# Run a Two Pool Model with feedback
M4=twopfFit(deco$timeSeries[,c(1,2)], initialCarbon = 100)
months=M4$SoilRmodel@times
Ct=getC(M4$SoilRmodel)
matplot(months,Ct, type="l",lty=1, col=c(2,3), ylim=c(0,100), ylab="Carbon remaining (%)")
legend("topright", c("Total","Fast","Slow"), lty=1, col=c(1,2,3), bty="n")
points(deco$timeSeries[,c(1,2)])
lines(months, rowSums(Ct), lwd=2)

k=M4$FMEmodel$par[c(1,2)]
A4=diag(-k)
A4[2,1]=M4$FMEmodel$par[3]
A4[1,2]=M4$FMEmodel$par[4]
u4=matrix(c(M4$FMEmodel$par[5], 1-M4$FMEmodel$par[5]), ncol=1)
SA4=systemAge(A=A4, u=u4)
TT4=transitTime(A=A4, u=u4)

# Run a Three Pool Parallel Model
M5=threeppFit(deco$timeSeries[,c(1,2)], initialCarbon = 100)
months=M5$SoilRmodel@times
Ct=getC(M5$SoilRmodel)
matplot(months,Ct, type="l",lty=1, col=c(2,3,4), ylim=c(0,100), ylab="Carbon remaining (%)")
legend("topright", c("Total","Pool 1","Pool 2","Pool 3"), lty=1, col=c(1,2,3,4), bty="n")
points(deco$timeSeries[,c(1,2)])
lines(months, rowSums(Ct), lwd=2)

k=M5$FMEmodel$par[c(1,2,3)]
A5=diag(-k)
u5=matrix(c(M5$FMEmodel$par[4], M5$FMEmodel$par[5], 1-M5$FMEmodel$par[4]-M5$FMEmodel$par[5]), ncol=1)
SA5=systemAge(A=A5, u=u5)
TT5=transitTime(A=A5, u=u5)

# Run a Three Pool Series Model
M6=threepsFit(deco$timeSeries[,c(1,2)], initialCarbon = 100)
months=M6$SoilRmodel@times
Ct=getC(M6$SoilRmodel)
matplot(months,Ct, type="l",lty=1, col=c(2,3,4), ylim=c(0,100), ylab="Carbon remaining (%)")
legend("topright", c("Total","Pool 1","Pool 2","Pool 3"), lty=1, col=c(1,2,3,4), bty="n")
points(deco$timeSeries[,c(1,2)])
lines(months, rowSums(Ct), lwd=2)

k=M6$FMEmodel$par[c(1,2,3)]
A6=diag(-k)
A6[2,1]=M6$FMEmodel$par[4]
A6[3,2]=M6$FMEmodel$par[5]
u6=matrix(c(M6$FMEmodel$par[6], M6$FMEmodel$par[7], 1-M6$FMEmodel$par[6]-M6$FMEmodel$par[7]), ncol=1)
SA6=systemAge(A=A6, u=u6)
TT6=transitTime(A=A6, u=u6)

# Plot all models in one graph
Mlist=list(M1,M2,M3,M4,M5,M6)
Rs=sapply(Mlist, function(x){rowSums(getC(x$SoilRmodel))})
modelNames=c("One-pool", "Two-pool parallel", "Two-pool series", "Two-pool feedback", "Three-pool parallel", "Three-pool series")
plot(deco$timeSeries[,c(1,2)], ylab="C remaining (%)")
matlines(months, Rs, col=2:6, lty=1)
legend("topright", modelNames, lty=1, col=2:6, bty="n")

# Models sum of squared residules (SSR), mean squared residuals (MSR) and AIC
statistics=data.frame(npar=sapply(Mlist, function(x){length(x$FMEmodel$par)}), SSR=sapply(Mlist, function(x){x$FMEmodel$ssr}), MSR=sapply(Mlist, function(x){x$FMEmodel$ms}), AIC=sapply(Mlist, function(x){x$AIC}))
row.names(statistics)<-modelNames

# Plot models in array
Ct=sapply(Mlist, function(x){getC(x$SoilRmodel)})
par(mfrow=c(2,2), mar=c(4,4.5,1,0.5))

matplot(months,Ct[[2]], type="l",lty=1, col=4:2, ylab=expression(paste("C remaining (%)")),
          main=modelNames[2], ylim=c(0,100), bty="n", font.main=1)
points(deco$timeSeries[,c(1,2)], pch=19, cex=0.5)
lines(months, Rs[,2])

matplot(months,Ct[[3]], type="l",lty=1, col=4:2, ylab=expression(paste("C remaining (%)")),
        main=modelNames[3], ylim=c(0,100), bty="n", font.main=1)
points(deco$timeSeries[,c(1,2)], pch=19, cex=0.5)
lines(months, Rs[,3])

matplot(months,Ct[[4]], type="l",lty=1, col=4:2, ylab=expression(paste("C remaining (%)")),
        main=modelNames[4], ylim=c(0,100), bty="n", font.main=1)
points(deco$timeSeries[,c(1,2)], pch=19, cex=0.5)
lines(months, Rs[,4])

matplot(months,Ct[[5]], type="l",lty=1, col=4:2, ylab=expression(paste("C remaining (%)")),
        main=modelNames[5], ylim=c(0,100), bty="n", font.main=1)
points(deco$timeSeries[,c(1,2)], pch=19, cex=0.5)
lines(months, Rs[,5])

legend("topright", c("Total","Pool 1", "Pool 2", "Pool 3"), lty=1, col=c(1,4:2), bty="n")
par(mfrow=c(1,1))
