#' Load all entries of the aridec dataset
#'
#' @param path character string with the path where isdb is stored
#' @return R list with all entries
#' @export
#' @examples
#' aridec=loadEntries()
loadEntries <- function(path="~/Repos/aridec/data/") {
  entryNames=list.dirs(path, full.names=FALSE, recursive=FALSE)

  longList=lapply(entryNames, FUN=readEntry, path=path)

 return(longList)
}


