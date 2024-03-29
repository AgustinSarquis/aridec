#' Fits a two pool model with series structure to a time series
#'
#' @param timeSeries A time series of entry values
#' @param initialCarbon The initial amount of carbon in units that correspond to the time series data
#' @param inipars vector of parameter values for the initial search of the optimization algorithm
#' @return R list with an FME model object, a SoilR model object, and the AIC value
#' @export
#' @import FME
#' @import graphics
#' @importFrom stats complete.cases
#' @examples
#' \donttest{
#' aridec=loadEntries()
#' entry=aridec[["20"]]
#' b=twopsFit(timeSeries = entry$timeSeries[,c(1,2)],
#' initialCarbon=100,
#' inipars=c(0.005, 0.00001, 0.1, 0.01))
#' }
twopsFit=function(timeSeries, initialCarbon, inipars=c(1, 0.5, 0.5, 0.3)){
  #  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Rt=cumsum(timeSeries[complete.cases(timeSeries),2]))
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  if(n < 5) stop("Time series is too short. No degrees of freedom")
  tt=seq(from=0, to=tail(complete[,1],1), length.out = 500)

  Func=function(pars){
    mod=SoilR::TwopSeriesModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], C0=initialCarbon*c(pars[4], 1-pars[4]), In=0)
    #    Rt=SoilR::getAccumulatedRelease(mod)
    Ct=SoilR::getC(mod)
    return(data.frame(time=tt, Ct=rowSums(Ct)))
  }

  costFunc=function(pars){
    output=Func(pars)
    return(modCost(model=output, obs=complete))
  }

  Fit=modFit(f=costFunc, p=inipars, method="Marq", lower=rep(0,4), upper=c(Inf, Inf, 1,1))
  bestMod=Func(pars=Fit$par)
  print(paste(c("k1=", "k2=", "a21=", "Proportion of C0 in pool 1="),Fit$par))
  plot(complete, ylim=c(0,1.2*max(complete[,2])))
  lines(bestMod)
  AIC=(2*length(Fit$par))-2*log(Fit$ms)
  print(paste("AIC = ",AIC))
  SoilRmodel=SoilR::TwopSeriesModel(t=tt,ks=Fit$par[1:2], a21=Fit$par[1]*Fit$par[3], C0=initialCarbon*c(Fit$par[4], 1-Fit$par[4]), In=0)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, AIC=AIC))
}
