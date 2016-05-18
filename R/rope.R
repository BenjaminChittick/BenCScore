"
Pipeline for Analysis of Primary High Throughput Screens
Benjamin Chittick

This script contains functions for choosing a region of practical equivalence

References:
http://www.indiana.edu/~kruschke/BEST/BEST.pdf
"

#library dependencies
library(data.table)

#' autoRope
#' 
#' From a provided maximum detectable EC50 and by assuming a one site effective 
#' concentration 50 curve this method estimates the minimum signal cutoff for an 
#' active hit at the given screening concentration.
#' 
#' @param assay.conc numeric value of assay concentration (same units as max.con)
#' @param max.conc numeric value of maximum detectable EC50 (same units as assay.conc)
#' @param modifier numeric modifier of signal
#' @export
#' @examples
#' autoRope(25, 100, 0.8)

autoRope <- function(assay.conc, max.conc, modifier=1, hillslope=1,
                     top=100, bottom=0){
  
  #assumes inhibition curve:
  #y = bottom + (top - bottom)/(1+10^((logIC50 - x)*hillslope))
  
  logIC50 <- log10(max.conc)
  x <- log10(assay.conc)
  
  pct.effect <- bottom + (top - bottom)/(1 + 10^((logIC50 - x)*hillslope))
  pct.effect <- pct.effect/modifier
  
  return(pct.effect)
}