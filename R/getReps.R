"
Pipeline for Analysis of Primary High Throughput Screens
Benjamin Chittick
"

#library dependencies
library(data.table)

#' repScatter
#' 
#' A function to tie together extracting compound and control replicates for 
#' plotting
#' 
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' replicates <- repScatter(mdmx)

repScatter <- function(screen){
  
  #checks
  checkScreen(screen)
  
  #add Well_Plate_Pair variable to screen
  screen <- wellPlatePair(screen)
   
  positive.replicates <- getReps(screen, 'positive control', 'Well_Plate_Pair')
  neutral.replicates <- getReps(screen, 'neutral control', 'Well_Plate_Pair')
  compound.replicates <- getReps(screen, 'compound', 'Compound_ID')
  
  scatter.data <- rbindlist(list(positive.replicates, 
                                 neutral.replicates, 
                                 compound.replicates), fill=TRUE)
  
  return(scatter.data)
}


#' getReps
#' 
#' Get replicates by factor
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @param role character string indicating which well type to subset data must
#' be one of: 'positive control', 'negative control', or 'compound' 
#' @param factor column name containing factor by which to cast data
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' compound.replicates <- getReps(mdmx, 'compound', 'Compound_ID')
#' 
#' neutral.replicates <- getReps(mdmx, 'neutral control', '')

getReps <- function(screen, role, variable){
  #split, factor, and cast data to generate screening replicates in wide format
  
  #checks
  checkScreen(screen)
  checkRole(role)
  
  #subset to well role
  subgroup <- screen[screen$Well_Type == role, ]
  
  #add factor handel and subset to just factor and values
  subgroup$Variable <- subgroup[, variable, with=FALSE]
  if (any(subgroup$Variable == '')) {
    warning('Blank values found in data factor variable. 
            Rows containing blank factor values were ignored')
    subgroup <- subgroup[subgroup$Variable != '', ]
  }
  subgroup <- subgroup[, c('Well_Type', 'Variable', 'Value'), with=FALSE]
  
  #split and apply casting variable
  replicates <- tapply(subgroup$Value, subgroup$Variable, 
                       function(x) paste('Replicate', 1:length(x), sep='_'))
  subgroup$Replicates <- unsplit(replicates, subgroup$Variable)
  
  #cast data
  replicates <- dcast.data.table(subgroup, Well_Type + Variable ~ Replicates,
                                 value.var='Value')
}