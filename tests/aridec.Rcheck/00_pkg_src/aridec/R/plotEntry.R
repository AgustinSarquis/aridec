#' Plot individual entries of the aridec dataset
#'
#' @param entry character string with the name of the entry to be plotted
#' @return A plot
#' @export
#' @examples
#' \dontrun{
#' aridec=loadEntries(path='/aridec/data/')
#' plotEntry(entry=aridec[["Adair2017"]])
#' }
plotEntry=function(entry){
  x=entry$timeSeries$Time
  n=ncol(entry$timeSeries)
  ys=entry$timeSeries[,2:n]

  matplot(x, ys,type="b",lty=1, pch=19, main=entry$citationKey,
          xlab=paste("Time (",entry$variables$V1$units, ")"),
          ylab=paste(entry$variables$V2$units))
  y.names=names(entry$timeSeries)[-1]
  legend("topright",y.names,col=1:(n-1),pch=19,lty=1,bty="n")
}
