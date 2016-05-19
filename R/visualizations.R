"
Pipeline for Analysis of Primary High Throughput Screens
Benjamin Chittick

This script contains visualizing functions for HTS analysis 
"

#library dependencies
library(ggplot2)

#' columnPlot
#' 
#' Make a column-wise plot of a high-throughput screen
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' columnPlot(mdmx)

columnPlot <- function(screen){
  
  #check screen
  checkScreen(screen)
  
  #guess limits
  upper <- quantile(screen$Value, 0.99, na.rm=TRUE)
  lower <- quantile(screen$Value, 0.01, na.rm=TRUE)
  
  ggplot(screen, aes(x=Column, y=Value, color=Well_Type)) + 
    scale_color_brewer(palette='Set1') +
    geom_point() + 
    facet_wrap(~Run_ID) + 
    ylim(lower, upper) + 
    ggtitle('Column-Wise Plot') + 
    xlab('Column') + 
    ylab('Value') +
    theme(axis.text=element_text(size=16), 
          axis.title=element_text(size=16), 
          legend.text=element_text(size=16),
          title=element_text(size=16),
          legend.title=element_blank())
}

#' scatterPlot
#' 
#' Make a replicate scatter plot of a primary high-throughput screen
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' scatterPlot(mdmx)

scatterPlot <- function(screen){
  
  #check screen
  checkScreen(screen)
  
  #the screening campaign must be subset to the primary
  if (any(!(screen$Run_Intent == 'PRIMARYSCREEN'))) {
    warning('Scatter plot only applies to primary screens. 
             Subsetting plot to primary screen')
    screen <- screen[screen$Run_Intent == 'PRIMARYSCREEN', ]
  }
  
  #get replicates
  if (all(c('neutral control', 'compound', 'positive control') %in% 
          unique(screen$Well_Type))) {
    replicates <- repScatter(screen)
  } else {
    #add Well_Plate_Pair variable to screen
    screen <- wellPlatePair(screen)
    for (well_type in unique(screen$Well_Type)) {
      if (well_type == 'compound') variable = 'Compound_ID'
      else variable = 'well_Plate_Pair'
      replicates <- getReps(screen, well_type, 'Well_Plate_Pair')
    }
  }

  #guess limits
  upper <- quantile(c(replicates$Replicate_1, replicates$Replicate_2), 0.999,
                    na.rm=TRUE)
  lower <- quantile(c(replicates$Replicate_1, replicates$Replicate_2), 0.0001,
                    na.rm=TRUE)
  
  #plot
  ggplot(replicates, aes(x=Replicate_1, y=Replicate_2, color=Well_Type)) + 
    scale_color_brewer(palette='Set1') +
    geom_point() + 
    xlab('Replicate 1') + 
    ylab('Replicate 2') + 
    xlim(lower, upper) + 
    ylim(lower, upper) +
    ggtitle('Replicate Scatter') + 
    theme(axis.text=element_text(size=16), 
          axis.title=element_text(size=16), 
          legend.text=element_text(size=16),
          title=element_text(size=16),
          legend.title=element_blank())
}

#' wellHeatMap
#' 
#' Make a heat map of well backgrounds
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' wellHeatMap(mdmx)

wellHeatMap <- function(screen){
  
  #check screen
  checkScreen(screen)
  
  #positional effects
  pattern <- platePattern(screen)
  
  #guess limits
  upper <- quantile(screen$Value, 0.99, na.rm=TRUE)
  lower <- quantile(screen$Value, 0.01, na.rm=TRUE)
  
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  ggplot(pattern, aes(factor(Column), factor(Row))) +
    facet_wrap('Run_ID') + 
    geom_tile(aes(fill=Background), color='white') +
    scale_fill_gradientn(colours = myPalette(100), limits=c(lower, upper)) +
    ggtitle('Plate Well Heatmap') + 
    xlab('Column') + 
    ylab('Row') +
    theme(axis.text=element_text(size=16), 
          axis.title=element_text(size=16), 
          legend.text=element_text(size=16),
          title=element_text(size=16),
          legend.title=element_blank())
}

#' logisticFit
#' 
#' 4 parameter logistic binding function
#' 
#' @param x concentration
#' @param IC50 log(IC50)
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' wellHeatMap(mdmx)

logisticFit <- function(x, top, bottom, IC50, hillslope){
  y.hat <- top + (bottom - top)/(1 + 10^((log(IC50) - x)*hillslope))
  return(y.hat)
}

#' plotCurve
#' 
#' Plot fitted binding curves
#' 
#' @param compound data.frame row with Compound_ID, IC50, IC50_Error, Hillslope,
#' Top, Bottom, Value, and Concentration.
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' wellHeatMap(mdmx)

plotCurve <- function(screen){
  
  IC50 <- compound$IC50[1]
  Hillslope <- compound$Hillslope[1]
  compound.id <- compound$Compound_ID[1]
  Top <- compound$Top[1]
  Bottom <- compound$Bottom[1]
  IC50_Error <- compound$IC50_Error[1]
  
  p <- ggplot(compound, aes(x=log(Concentration), y=Value)) + 
    geom_point() + 
    stat_function(fun=logisticFit, args=list(Top, Bottom, IC50, Hillslope)) +
    ylab('% Probe Bound') +
    ggtitle(compound.id) + 
    geom_text(aes(x=3.5, y=20, label = paste('IC50', round(IC50, 0), '+/-', 
                                             round(IC50_Error, 1))), size=3) + 
    geom_text(aes(x=3.5, y=10, label=paste('slope', round(Hillslope, 1))), 
              size=3) + 
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=10), 
          legend.text=element_text(size=10),
          title=element_text(size=10))
  
  print(p)
}

#' screenPlot
#' 
#' Plot fitted binding curves for a concentration-response screen
#' 
#' @param screen data.frame row with Compound_ID, Value, and Concentration.
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' screenPlot(confirmation)

screenPlot <- function(screen){
 
  confirmation <- screen[screen$Run_Intent == 'CONFIRMATION_IN_DOSE', ]
  confirmation <- screenFit(confirmation)
  
  cpd <- confirmation[confirmation$Well_Type == 'compound', ]
  cpd <- cpd[cpd$Compound_ID != '', ]
  
  x.min <- min(cpd$Concentration)
  x.max <- max(cpd$Concentration)
  y.min <- min(cpd$Value)
  y.max <- max(cpd$Value)
  
  cpd <- cpd[order(cpd$IC50), ]
  id <- unique(cpd$Compound_ID)
  cpd <- split(cpd, cpd$Compound_ID)
  
  curves <- lapply(id, function(x) plotCurves(cpd[[x]], x.min, x.max, y.min, y.max))
}