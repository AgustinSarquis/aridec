#' Creates a data frame with the coordinates of the sites
#'
#' @param database A list with the aridec structure
#' @return A data frame with the longitude and latitude of sites from the database
#' @export
#' @examples
#' coor=coordinates(database=aridec)
coordinates <- function(database) {
  long=lapply(database, FUN=function(x){x$siteInfo$coordinates$longitude})
  lat=lapply(database, FUN=function(x){x$siteInfo$coordinates$latitude})
  return(data.frame(longitude=unlist(long), latitude=unlist(lat)))
}