#' Load all entries of the aridec dataset
#'
#' @param path character string with the path where aridec data is stored
#' @return R list with all entries
#' @export
#' @examples
#' aridec=loadEntries()
loadEntries <- function(path="~/aridec/data/") {
  entryNames=list.dirs(path, full.names=FALSE, recursive=FALSE)

  longList=lapply(entryNames, FUN=readEntry, path=path)
  names(longList)<-entryNames

 return(longList)
}


