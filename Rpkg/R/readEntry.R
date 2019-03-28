#' Read single entry of the aridec database
#'
#' @param path character string with the path where aridec is stored
#' @param entryName character string with the name of the entry in the database
#' @return R list with the entry
#' @export
#' @import yaml
#' @examples
#' Adair2017=readEntry(entryName="Adair2017")
readEntry <- function(path="~/Repos/aridec/data/", entryName) {

    entry=yaml::yaml.load_file(input=paste(path,entryName,"/metadata.yaml",sep=""))
    csv=read.csv(file=paste(path,entryName,"/data.csv",sep=""))
    entry[["data"]]<-csv
    assign(entryName, entry)

    return(entry)
}
