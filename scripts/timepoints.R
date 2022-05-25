# if you just want to know the number of time points in a specific mass loss variable
timepoints=lapply(db, function(x){length(na.omit(x$timeSeries[,2]))})

# to Create a data frame with the number of points in time per entry
timepoints <- function(dataframe){ # creates a function applied to a single dataframe
  if (anyNA(dataframe) == TRUE) {
    print(NA)
  } else {
    nrow(dataframe)
  }
}

timePoints <- function(database) { # creates a function applied to the whole database
  timePoints=lapply(database, FUN=function(x){timepoints(x$timeSeries)})
  return(data.frame(timePoints=unlist(timePoints)))
}

# to Create a data frame with the number of variables per entry
variables <- function(dataframe){ # creates a function applied to a single dataframe
    ncol(dataframe)
  }

Variables <- function(database) { # creates a function applied to the whole database
  Variables=lapply(database, FUN=function(x){variables(x$timeSeries)})
  return(data.frame(Variables=unlist(Variables)))
}

# timePoints=timePoints(database)

timePoints2=na.omit(timePoints) # omits entries with NAs

# to obtain data points per variable...
actualnrow <- function(x){
  NROW(na.omit(x))
}
groupnrow <- function(x){
  lapply(x[,-1], actualnrow)
}
totalnrow <- function(db){
  datapoints=lapply(db, FUN=function(x){groupnrow(x$timeSeries)})
  return(data.frame(datapoints=unlist(datapoints)))
}

# totalpointsaridec=totalnrow(aridec)

# to identify which entries have a certain number of time points (change number)

which(timePoints$timePoints==4)

# maximum value?
max(timePoints$timePoints, na.rm=TRUE)
