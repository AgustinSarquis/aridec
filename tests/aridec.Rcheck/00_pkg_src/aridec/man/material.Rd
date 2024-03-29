% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/material.R
\name{material}
\alias{material}
\title{Creates a data frame with the list of litter samples' plant parts}
\usage{
material(database)
}
\arguments{
\item{database}{A list with the aridec structure}
}
\value{
A data frame with the list of the litter samples' plant parts from the database
}
\description{
Creates a data frame with the list of litter samples' plant parts
}
\examples{
\dontrun{
aridec=loadEntries(path='/aridec/data/')
material=material(database=aridec)
}
}
