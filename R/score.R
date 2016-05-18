"
Pipeline for Analysis of Primary High Throughput Screens
Benjamin Chittick

This script contains functions for scoring compounds from high-throughput 
screens
"

#library dependencies
library(data.table)

#' zScore
#' 
#' This function calcluates the compounds z-score
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' mdmx <- zScore(mdmx)
#' head(mdmx)

zScore <- function(screen, alternative){
  
  #check screen
  checkScreen(screen)
  checkAlternative(alternative)
  
  #check for plate stat
  if(!all(c('Number_Observations', 'Plate_Median', 'Plate_MAD') 
          %in% names(screen))){
    screen <- plateStats(screen)
  }
  
  #add robust z-score
  screen$Z <- (screen$Value - screen$Plate_Median)/screen$Plate_MAD/1.4826
  #take mean of compound z-scores across plates
  Z <- tapply(screen$Z, screen$Compound_ID, 
                    function(x) mean(x, na.rm=TRUE))
  screen$Z <- unsplit(Z, screen$Compound_ID)
  
  #p-value
  if (alternative == 'less'){
    screen$P_Value_Z <- dnorm(screen$Z, 0, 1)
    screen$P_Value_Z[screen$Z > 0] <- 0.5
  }
  if (alternative == 'greater'){
    screen$P_Value_Z <- dnorm(screen$Z, 0, 1)
    screen$P_Value_Z[screen$Z < 0] <- 0.5
  }
  
  return(screen)
}

#' RVMtTest
#' 
#' This function calcluates the compounds t-statistic from a Random Variance
#' Model t-test
#' 
#' References:
#' Malo, N. et al. Experimental design and statistical methods for improved hit 
#' detection in high-throughput screening. J. Biomol. Screen. 15, 990–1000 
#' (2010).
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' mdmx <- RVMtTest(mdmx)
#' head(mdmx)

RVMtTest <-function(screen, alternative, c=0){

  #check screen
  checkScreen(screen)
  checkAlternative(alternative)
  
  #check for plate stat
  if(!all(c('Number_Observations', 'Plate_Median', 'Plate_MAD') 
          %in% names(screen))){
    screen <- plateStats(screen)
  }
  
  #get alpha and beta
  cpd.ab <- getAB(screen, 'compound', 'Compound_ID')[[1]]
  
  #calculate median signal for compounds
  #if there are only 2 replicates the mean is returned
  Compound_Mean <- tapply(screen$Value, screen$Compound_ID, 
                     function(x) median(x, na.rm=TRUE))
  screen$Compound_Mean <- unsplit(Compound_Mean, screen$Compound_ID)
  
  #calculate compound variance
  Compound_Variance <- tapply(screen$Value, screen$Compound_ID, 
                    function(x) var(x, na.rm=TRUE))
  screen$Compound_Variance <- unsplit(Compound_Variance, screen$Compound_ID)
  #catach NA variance
  screen$Compound_Variance[is.na(screen$Compound_Variance)] <- 0 
  
  #2*cpd.ab[2] is used because EstimateAB uses 1/b parameter
  screen$Pooled_SD <- sqrt(((screen$Number_Observations - 1)*
                             screen$Compound_Variance + 2*cpd.ab[2])/
                            (screen$Number_Observations - 1 + 2*cpd.ab[1]))
  
  #default constant of c=0 is chosen for the RVM one-sample t-test
  screen$RVM_t <- sqrt(screen$Number_Observations)*(screen$Compound_Mean-c)/
    screen$Pooled_SD
  
  #p-value
  if (alternative == 'less'){
    screen$P_Value_RVM_t <- dt(screen$RVM_t, 
                              df=(screen$Number_Observations - 1 + 2*cpd.ab[1]))
    screen$P_Value_RVM_t[screen$RVM_t > 0] <- 0.5
  }
  if (alternative == 'greater'){
    screen$P_Value_RVM_t <- dt(screen$RVM_t, 
                              df=(screen$Number_Observations - 1 + 2*cpd.ab[1]))
    screen$P_Value_RVM_t[screen$RVM_t < 0] <- 0.5
  }
  
  return(screen)
}

#' nBayes
#' 
#' Deprecated: Estimate the probability of a compound coming from the neutral
#' control distribution assuming the compounds and neutral controls are normally 
#' distributed with unknown mu and tau.
#' 
#' @param screen A data frame or data table containing screening data with 
#' variables: Plate_ID, Row, Column, Compound_ID, Well_Type, Plate_Mask, 
#' Well_Mask, Value. For an explination of variables see read me.
#' @param alternative A character string specifying the alternet hypothesis. 
#' Must be one of "two.sided", "greater", or "less". 
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' mdmx <- nBayes(mdmx, alternative='less')
#' head(mdmx)

nBayes <- function(screen, alternative, ROPE, mu0=0.0, n0=0.01, size=20000){

  #checks
  checkScreen(screen)
  checkAlternative(alternative)
  
  #add Well_Plate_Pairs
  screen <- wellPlatePair(screen)
  
  #build null.post
  #the null priors are built by calling replicates across plates by well
  #position similar to how compounds are replicated by stamping two different
  #assay plates

  neg.ab <- getAB(screen, 'neutral control', 'Well_Plate_Pair')[[1]]
  neutral.controls <- screen[screen$Well_Type == 'neutral control', ]
  null.post <- nSample(neutral.controls$Value, a=neg.ab[1], b=neg.ab[2],
                       mu0=mu0, n0=n0, size=size)
  
  #contrast compounds to null posterior
  compound.replicates <- getReps(screen, 'compound', 'Compound_ID')
  replicates <- grep('Replicate_[0-9]+', names(compound.replicates))
  cpd.ab <- getAB(screen, 'compound', 'Compound_ID')[[1]]
  
  cpd.prob <- apply(compound.replicates[, replicates, with=FALSE], 1, 
                    function(x) nContrast(x, a=cpd.ab[1], b=cpd.ab[2], mu0=mu0,
                                          n0=n0, null.post=null.post, 
                                          alternative, size, ROPE))
  probs <- rbindlist(cpd.prob)
  probs$Compound_ID <- compound.replicates$Variable
  screen <- merge(screen, probs, by='Compound_ID', all.x=TRUE)
}

#' tBayes
#' 
#' Estimate the probability of a compound coming from the neutral control 
#' distribution assuming the compounds and neutral controls are t distributed 
#' with unknown mean, and variance. The compound means are assumed to come from
#' t-distribution fit to all compounds in the screen.
#' 
#' @param screen A data frame or data table containing screening data with 
#' variables: Plate_ID, Row, Column, Compound_ID, Well_Type, Plate_Mask, 
#' Well_Mask, Value. For an explination of variables see read me.
#' @param alternative Character string specifying the alternet hypothesis. 
#' Must be one of "two.sided", "greater", or "less". 
#' @param runjags.method Character string specifying the method by which jags is 
#' run, default is set to 'simple' to avoid conflict with snow library
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' mdmx <- tBayes(mdmx, alternative='less')
#' head(mdmx)

tBayes <- function(screen, alternative, ROPE){

  #checks
  checkScreen(screen)
  checkAlternative(alternative)
  
  #add Well_Plate_Pairs
  #find wellPlatePair and move to plateStats, wrap it in an if statment
  screen <- wellPlatePair(screen)
  
  #the null priors are built by calling replicates across plates
  neutral.controls <- getReps(screen, 'neutral control', 'Well_Plate_Pair')
  idx.neutral.reps <- grep('Replicate_[0-9]+', names(neutral.controls))
  neutral.replicates <- neutral.controls[, idx.neutral.reps, with=FALSE]
  neg.data <- getAB(screen, 'neutral control', 'Well_Plate_Pair')
  neg.ab <- neg.data[[1]]
  
  mu.prior <- median(unlist(neutral.replicates), na.rm=TRUE)
  var.prior <- var(apply(neutral.replicates, 1,
                         function(x) mean(x, na.rm=TRUE)), na.rm=TRUE)
  nu.prior <- median(neg.data[[3]])

  #sample neutral control posterior
  if (require(snow)) detach(package:snow, unload=TRUE)
  null.post <- tSample(unlist(neutral.replicates), mu.prior, var.prior, 
                       nu.prior, a=neg.ab[1], b=neg.ab[2])
  
  #contrast compounds to null posterior
  compounds <- getReps(screen, 'compound', 'Compound_ID')
  idx.cpd.reps <- grep('Replicate_[0-9]+', names(compounds))
  compound.replicates <- compounds[, idx.cpd.reps, with=FALSE]
  cpd.data <- getAB(screen, 'compound', 'Compound_ID')
  cpd.ab <- cpd.data[[1]]
  
  mu.prior <- median(unlist(compound.replicates), na.rm=TRUE)
  var.prior <- var(apply(compound.replicates, 1,
                         function(x) mean(x, na.rm=TRUE)), na.rm=TRUE)
  nu.prior <- median(cpd.data[[3]])
  
  cpd.prob <- tParallel(data=compound.replicates,
                        mu.prior=mu.prior,
                        var.prior=var.prior,
                        nu.prior=nu.prior,
                        a=cpd.ab[1], b=cpd.ab[2], null.post=null.post,
                        ROPE=ROPE, alternative=alternative)
    
  probs <- rbindlist(cpd.prob)
  probs$Compound_ID <- compounds$Variable
  screen <- merge(screen, probs, by='Compound_ID', all.x=TRUE)
  
  return(screen)
}