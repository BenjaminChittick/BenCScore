#' checkScreen
#' 
#' Check screen data for necessary data labels and values
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' checkScreen(mdmx)

checkScreen <- function(screen){
  
  #check for necessary column names
  if(!all(c('Plate_ID', 'Row', 'Column',
            'Compound_ID', 'Well_Type', 
            'Run_ID', 'Compound_Plate', 'Value') %in% names(screen))){
    stop('Screen data must contain columns: Plate_ID, Row, Column, Compound_ID,
         Well_Type, Run_ID, Compound_Plate, Value. For details see 
         documentation')
  }
  
  #valid well types are: compound, positive control, neutral control
  if(!all(screen$Well_Type %in% 
          c('compound', 'positive control', 'neutral control'))){
    stop('Valid Well Types are: compound, positive control, neutral control.
         Change Well Types so they are only: compound, positive control, 
         or neutral control')
  }
}

#' checkAlternative
#' 
#' Check alternative flag is one of 'less' or 'greater'
#' 
#' @param alternative Character string specifying the alternative hypothesis
#' must be one of: 'less', 'greater'
#' @export
#' @examples
#' checkAlternative('less')

checkAlternative <- function(alternative){
  #check alternative hypothesis
  if(!(alternative %in% list('greater', 'less'))){
    stop('alternative must be one of: "greater", "less"')
  }
}

#' checkRole
#' 
#' Check role flag is one of 'positive control', 'neutral control', 'compound'
#' 
#' @param role Character string specifying well role must be one of: 'less', 
#' 'greater'
#' @export
#' @examples
#' Load screen
#' checkRole('positive control')

checkRole <- function(role){
  if (!(role %in% list('positive control', 'neutral control', 'compound'))) {
    stop('Invalid role for Well_Type')
  }
}