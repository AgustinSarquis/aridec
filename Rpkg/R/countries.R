#' Creates a data frame with the countries of the sites
#'
#' @param database A list with the aridec structure
#' @return A data frame with the countries from the database
#' @export
#' @examples
#' countries=countries(database=aridec)
countries <- function(database) {
  country=lapply(database, FUN=function(x){x$siteInfo$country})
  return(data.frame(country=unlist(country)))
}
