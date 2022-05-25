#' Creates a data frame with soil orders of the sites
#'
#' @param database A list with the aridec structure
#' @return A data frame with the soil orders from the database
#' @export
#' @examples
#' \dontrun{
#' aridec=loadEntries(path='/aridec/data/')
#' soilorder=soilorder(database=aridec)
#' }
soilorder <- function(database) {
  order=lapply(database, FUN=function(x){x$siteInfo$soilTaxonomy$soilOrder})
  return(data.frame(order=unlist(order)))
}
