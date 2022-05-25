#' Creates a data frame with the nitrogen content of litter samples
#'
#' @param database A list with the aridec structure
#' @return A data frame with the nitrogen content (%) of the litter samples from the database
#' @export
#' @examples
#' \dontrun{
#' aridec=loadEntries(path='/aridec/data/')
#' N=nitrogen(database=aridec)
#' }
nitrogen = function(database) {
  mean=lapply(database, FUN=function(x){x$initConditions$nitrogen})
  return(data.frame(nitrogen=unlist(mean)
  ))
}
