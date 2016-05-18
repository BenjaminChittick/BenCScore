"
Pipeline for Analysis of Primary High Throughput Screens
Benjamin Chittick

This script contains functions for scoring compounds from high-throughput 
screens
"

#library dependencies
library(data.table)

#' addModels
#' 
#' Add z-score, RVM t-test, Bayesian RVM, and custom Bayesian models to HTS
#' screening data
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' @param alternative String indicating direction of comparison. Must be either
#' 'less' or 'greater'
#' @param ROPE Region of Practical Equivalence a practical cutoff at which a
#' compound could be active but so minimally so we do not care
#' @param c Numeric value for the inactive popluation in RVM t-test. Use same
#' value as ROPE.
#' values
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' mdmx <- applyModels(mdmx)

addModels <- function(screen, alternative, ROPE, c=0){
  
  #check screen
  checkScreen(screen)
  checkAlternative(alternative)
  
  screen <- zScore(screen, alternative)
  screen <- RVMtTest(screen, alternative, c)
  screen <- tBayes(screen, alternative, ROPE)
  
  return(screen)
}
