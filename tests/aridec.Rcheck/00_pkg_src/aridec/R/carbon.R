#' Creates a data frame with the carbon content in litter samples
#'
#' @param database A list with the aridec structure
#' @return A data frame with the carbon content (%) of the litter samples from the database
#' @export
#' @examples
#' \dontrun{
#' aridec=loadEntries(path='/aridec/data/')
#' C=carbon(database=aridec)
#' }
carbon= function(database){
  mean=lapply(database, FUN=function(x){x$initConditions$carbon})
  return(data.frame(carbon=unlist(mean)
  ))
}
