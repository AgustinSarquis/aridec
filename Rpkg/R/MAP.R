#' Creates a data frame with mean annual precipitation values of the sites
#'
#' @param database A list with the aridec structure
#' @return A data frame with the mean annual precipitation values (mm) from the database
#' @export
#' @examples
#' aridec<-loadEntries(path='~/Repos/aridec/data/')
#' MAP=MAP(database=aridec)
MAP <- function(database) {
  MAP=lapply(database, FUN=function(x){x$siteInfo$MAP})
  return(data.frame(MAP=unlist(MAP)))
}

