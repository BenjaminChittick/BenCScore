"
Pipeline for Analysis of Primary High Throughput Screens
Benjamin Chittick
"

#library dependencies

#' platePattern
#' 
#' Create a heat map of plate wells by run
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @export
#' @examples
#' library(RColorBrewer)
#' library(ggplot)
#' Load screen
#' data(mdmx)
#' pattern <- platePattern(mdmx)
#' myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
#' ggplot(pattern, aes(factor(Column), factor(Row))) +
#'  facet_wrap('Run_ID') + 
#'  geom_tile(aes(fill=Background), color='white') +
#'  scale_fill_gradientn(colours = myPalette(100), limits=c(-150, 100)) +
#'  ggtitle('Plate Well Heatmap') + 
#'  xlab('Column') + 
#'  ylab('Row') +
#'  theme(axis.text=element_text(size=16), 
#'        axis.title=element_text(size=16), 
#'        legend.text=element_text(size=16),
#'        title=element_text(size=16),
#'        legend.title=element_blank())


platePattern <- function(screen){
  
  #check screen
  checkScreen(screen)
  
  #zero pad row and column index
  screen$Row <- sprintf('%02d', screen$Row)
  screen$Column <- sprintf('%02d', screen$Column)
  
  screen$Well <- paste(screen$Row, screen$Column, sep='-')
  
  pattern <- data.frame(Row=rep(1:max(screen$Row), each=max(screen$Column)),
                        Column=rep(1:max(screen$Column), times=max(screen$Row)))
  
  #zero pad row and column
  pattern$Row <- sprintf('%02d', pattern$Row)
  pattern$Column <- sprintf('%02d', pattern$Column)
  pattern$Well <- paste(pattern$Row, pattern$Column, sep='-')
  
  #split data by run and calculate well median
  runData <- split(screen, screen$Run_ID)
  runPattern <- lapply(runData, 
                       function(x) tapply(x$Value, x$Well, 
                                          function(a) median(a, na.rm=TRUE)))
  runPattern <- lapply(runPattern, 
                       function(x) data.frame('Well'=names(x),'Value'=x))
  
  #merge run pattern with default pattern to catch masked wells
  runPattern <- lapply(runPattern, 
                       function(x) merge(pattern, x, by='Well', all.x=TRUE))
  for (i in 1:length(runPattern)){
    runPattern[[i]]$Run_ID <- names(runPattern[i])
  }
  
  runPattern <- rbindlist(runPattern)
  setnames(runPattern, 'Value', 'Background')
  
  return(runPattern)
}

#' fixPattern
#' 
#' fix patterns on plates by run
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @export
#' @examples
#' library(RColorBrewer)
#' library(ggplot)
#' Load screen
#' data(mdmx)
#' pattern <- platePattern(mdmx)
#' myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
#' ggplot(pattern, aes(factor(Column), factor(Row))) +
#'  facet_wrap('Run_ID') + 
#'  geom_tile(aes(fill=Background), color='white') +
#'  scale_fill_gradientn(colours = myPalette(100), limits=c(-150, 150)) +
#'  ggtitle('Plate Well Heatmap') + 
#'  xlab('Column') + 
#'  ylab('Row') +
#'  theme(axis.text=element_text(size=16), 
#'        axis.title=element_text(size=16), 
#'        legend.text=element_text(size=16),
#'        title=element_text(size=16),
#'        legend.title=element_blank())
#' 
#' mdmx <- fixPattern(mdmx)
#' pattern <- platePattern(mdmx)
#' ggplot(pattern, aes(factor(Column), factor(Row))) +
#'  facet_wrap('Run_ID') + 
#'  geom_tile(aes(fill=Background), color='white') +
#'  scale_fill_gradientn(colours = myPalette(100), limits=c(-150, 150)) +
#'  ggtitle('Plate Well Heatmap') + 
#'  xlab('Column') + 
#'  ylab('Row') +
#'  theme(axis.text=element_text(size=16), 
#'        axis.title=element_text(size=16), 
#'        legend.text=element_text(size=16),
#'        title=element_text(size=16),
#'        legend.title=element_blank())


fixPattern <- function(screen){
  
  #checks
  checkScreen(screen)
  
  #This function cannot be run twice
  if ('Background' %in% names(screen)) {
    stop('Background variable already exists. 
         Has pattern correction already been run?
         To run fixPattern remove or change Background column name.')
  }
  
  byRun <- function(run){

    pattern <- platePattern(run)
    pattern <- pattern[, c('Well', 'Run_ID', 'Background'), with=FALSE]
      
    #create well id 
    run$Well <- paste(sprintf('%02d', run$Row),
                      sprintf('%02d', run$Column), sep='-')
      
    run <- merge(run, pattern, by=c('Well', 'Run_ID'))
      
    positive.median <- median(run$Value[run$Well_Type == 'positive control'], 
                              na.rm=TRUE)
 
    #check if run has multiple plates
    if (length(unique(run$Plate_ID)) > 1){
      #remove well position background
      run$Value <- run$Value - run$Background
      
      #add back positive control median
      run$Value[run$Well_Type == 'positive control'] <- 
        run$Value[run$Well_Type == 'positive control'] + positive.median
    }
    
    return(run)
  }

  run.data <- split(screen, screen$Run_ID)
  test <- lapply(run.data, byRun)
  screen <- rbindlist(lapply(run.data, byRun))

  return(screen)
}