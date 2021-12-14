#' Creates a data frame with the species list of litter samples
#'
#' @param database A list with the aridec structure
#' @return A data frame with the list of species of the litter samples from the database
#' @export
#' @examples
#' species=species(database=aridec)
species= function(database){
  species=lapply(database, FUN=function(x){x$initConditions$species})
  return(data.frame(species=unlist(species) 
  ))
}
