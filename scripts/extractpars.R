extractpars <- function(list, par) {
  pars=lapply(list, FUN=function(x){x[par]}) 
  return(data.frame(pars=unlist(pars)))
}


