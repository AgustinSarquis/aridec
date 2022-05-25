#' Creates a data frame with the list of litter samples' plant parts
#'
#' @param database A list with the aridec structure
#' @return A data frame with the list of the litter samples' plant parts from the database
#' @export
#' @examples
#' \dontrun{
#' aridec=loadEntries(path='/aridec/data/')
#' material=material(database=aridec)
#' }
material= function(database){
  material=lapply(database, FUN=function(x){x$initConditions$material})
  return(data.frame(material=unlist(material)
  ))
}
