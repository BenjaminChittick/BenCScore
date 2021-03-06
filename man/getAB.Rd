% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getAB.R
\name{getAB}
\alias{getAB}
\title{getAB}
\usage{
getAB(screen, role, variable)
}
\arguments{
\item{screen}{A data frame or data table containing screening data with 
columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
Compound_Plate, Value. See documentation for further description of column
values}

\item{role}{character string indicating which well type to subset data must
be one of: 'positive control', 'negative control', or 'compound'}
}
\description{
This function splits screening data by a user supplied factor e.g. a compound
ID and estimates a prior distribution of variance from replicates using the 
estimateAB function
}
\examples{
Load screen
data(mdmx)
ab <- getAB(mdmx, 'compound', 'Compound_ID')[[1]]
}

