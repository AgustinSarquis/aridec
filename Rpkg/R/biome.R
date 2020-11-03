#' Creates a data frame with the ecosystem type of the sites
#'
#' @param database A list with the aridec structure
#' @return A data frame with the ecosystem types from the database
#' @export
#' @examples
#' biome=biome(database=aridec)
biome <- function(database) {
  biome=lapply(database, FUN=function(x){x$siteInfo$landCover})
  return(data.frame(biome=unlist(biome)))
}

