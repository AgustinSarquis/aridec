#' Creates a data frame with elevation values of the sites
#'
#' @param database A list with the aridec structure
#' @return A data frame with the elevation values (mm) from the database
#' @export
#' @examples
#' elevation=elevation(database=aridec)
elevation <- function(database) {
  elevation=lapply(database, FUN=function(x){x$siteInfo$elevation})
  return(data.frame(elevation=unlist(elevation)))
}
