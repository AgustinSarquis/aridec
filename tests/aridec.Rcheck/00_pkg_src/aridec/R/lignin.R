#' Creates a data frame with the lignin content in litter samples
#'
#' @param database A list with the aridec structure
#' @return A data frame with the lignin content (%) of the litter samples from the database
#' @export
#' @examples
#' \dontrun{
#' aridec=loadEntries(path='/aridec/data/')
#' lignin=lignin(database=aridec)
#' }
lignin= function(database){
  mean=lapply(database, FUN=function(x){x$initConditions$lignin})
  return(data.frame(lignin=unlist(mean)
  ))
}
