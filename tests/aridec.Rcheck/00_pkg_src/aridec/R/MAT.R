#' Creates a data frame with mean annual temperature values of the sites
#'
#' @param database A list with the aridec structure
#' @return A data frame with the mean annual temperature values (Celsius degrees) from the database
#' @export
#' @examples
#' \dontrun{
#' aridec=loadEntries(path='/aridec/data/')
#' MAT=MAT(database)
#' }
MAT <- function(database) {
  MAT=lapply(database, FUN=function(x){x$siteInfo$MAT})
  return(data.frame(MAT=unlist(MAT)))
}
