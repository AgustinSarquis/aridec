#' Plot individual entries of the aridec dataset
#'
#' @param entry character string with the name of the entry to be plotted
#' @return A plot
#' @export
#' @examples
#' plotEntry(entry=aridec[["Adair2017"]])
plotEntry=function(entry){
  x=entry$data$Time
  n=ncol(entry$data)
  ys=entry$data[,2:n]
  
  matplot(x, ys,type="b",lty=1, pch=19, main=entry$citationKey,
          xlab=paste("Time (",entry$Variables$V1$Units, ")"),
          ylab=paste(entry$Variables$V2$Units))
  y.names=names(entry$data)[-1]
  legend("topright",y.names,col=1:(n-1),pch=19,lty=1,bty="n")
}
