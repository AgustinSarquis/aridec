#' Read single entry of the aridec database
#'
#' @param path character string with the path where aridec is stored
#' @param entryName character string with the name of the entry in the database
#' @return R list with the entry
#' @export
#' @import yaml
#' @import utils
#' @examples
#' \dontrun{
#' Adair2017=readEntry(path = '~/aridec/data/', entryName="Adair2017")
#' }
readEntry <- function(path, entryName) {

    entry=yaml::yaml.load_file(input=paste(path,entryName,"/metadata.yaml",sep=""))
    csv=utils::read.csv(file=paste(path,entryName,"/timeSeries.csv",sep=""))
    entry[["timeSeries"]]<-csv
    init=utils::read.csv(file=paste(path,entryName,"/initConditions.csv",sep=""))
    entry[["initConditions"]]<-init
    assign(entryName, entry)

    return(entry)
}

