"
Pipeline for Analysis of Primary High Throughput Screens
Benjamin Chittick

This script contains miscellaneous functions for scoring compounds from 
high-throughput screens
"

#library dependencies
library(data.table)

#' plateStats
#' 
#' Add number of observations, plate median, and plate MAD
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' mdmx <- plateStats(mdmx)
#' head(mdmx)

plateStats <- function(screen){
  
  #check screen
  checkScreen(screen)
  
  #plate statistics
  #number of observations
  n.obs <- tapply(screen$Value, screen$Compound_ID, length)
  screen$Number_Observations <- unsplit(n.obs, screen$Compound_ID)
  
  #subset compounds
  cpd <- screen[screen$Well_Type == 'compound', ]
  cpd <- screen[screen$Compound_ID != '', ]
  
  #median/mad
  plate.median <- tapply(cpd$Value, cpd$Plate_ID, 
                         function(x) median(x, na.rm=TRUE))
  plate.mad <- tapply(cpd$Value, cpd$Plate_ID, 
                      function(x) mad(x, na.rm=TRUE))
  screen$Plate_Median <- unsplit(plate.median, screen$Plate_ID)
  screen$Plate_MAD <- unsplit(plate.mad, screen$Plate_ID)
  
  return(screen)
}

#' wellPlatePair
#' 
#' Concatenate Compound_Plate and Well coordinates to create new Well_Plate_Pair
#' column
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' mdmx <- wellPlatePair(mdmx)

wellPlatePair <- function(screen){
  #create a Well_Plate_Pair variable from the screening data
  
  screen$Well <- paste(sprintf('%02d', screen$Row),
                       sprintf('%02d', screen$Column), sep='-')
  
  screen$Well_Plate_Pair <- paste(screen$Well, screen$Compound_Plate, sep='-')
  
  return(screen)
}

#' percentEffect
#' 
#' Calculate percent effect relative to positive and neutral controls
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' mdmx <- percentEffect(mdmx)

percentEffect <- function(screen){

  pct_effect <- function(plate){
    neutral <- median(plate$Value[plate$Well_Type == 'neutral control'])
    positive <- median(plate$Value[plate$Well_Type == 'positive control'])
    
    plate$Value <- (plate$Value - neutral)/abs(positive - neutral)*100
    return(plate)
  }
  
  screen <- split(screen, screen$Plate_ID)
  screen <- lapply(screen, pct_effect)
  screen <- rbindlist(screen)
}