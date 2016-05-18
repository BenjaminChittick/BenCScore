"
Pipeline for Analysis of Primary High Throughput Screens
Benjamin Chittick

This script contains functions for estimating the prior distribution of variance

References:
Estimating the distribution of variance:
Wright, G. W., & Simon, R. M. (2003). A random variance model for detection of
differential gene expression in small microarray experiments. Bioinformatics, 
19(18), 2448–2455. doi:10.1093/bioinformatics/btg345

Rational for priors:
http://www.indiana.edu/~kruschke/BEST/BEST.pdf
"

#library dependencies

#' getAB
#' 
#' This function splits screening data by a user supplied factor e.g. a compound
#' ID and estimates a prior distribution of variance from replicates using the 
#' estimateAB function
#' 
#' @param screen A data frame or data table containing screening data with 
#' columns: Plate_ID, Row, Column, Compound_ID, Well_Type, Run_ID, 
#' Compound_Plate, Value. See documentation for further description of column
#' values
#' @param role character string indicating which well type to subset data must
#' be one of: 'positive control', 'negative control', or 'compound' 
#' @export
#' @examples
#' Load screen
#' data(mdmx)
#' ab <- getAB(mdmx, 'compound', 'Compound_ID')[[1]]

getAB <- function(screen, role, variable){
  
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
  
  #variance and degrees of freedom
  variable.var <- tapply(subgroup$Value, subgroup$Variable, 
                       function(x) var(x, na.rm=TRUE))
  variable.degf <- tapply(subgroup$Value, subgroup$Variable, 
                        function(x) length(x) - 1)
  
  variable.degf <- variable.degf[!is.na(variable.var) & variable.var != 0]
  variable.var <- variable.var[!is.na(variable.var) & variable.var != 0]
  
  #use EstimateAB to find compound alpha and beta paramters
  ab <- estimateAB(variable.var, variable.degf)
  
  return(list(ab, variable.var, variable.degf))
}

#' estimateAB
#' 
#' This function fits the alpha (a) and beta (b) parameters for the inverse 
#' gamma distribution of variance by maximum likelihood estimation. The log 
#' likelihood is calculated by the flik function.
#' 
#' Reference:
#' Wright, G. W., & Simon, R. M. (2003). A random variance model for detection 
#' of differential gene expression in small microarray experiments. 
#' Bioinformatics, 19(18), 2448–2455. doi:10.1093/bioinformatics/btg345
#' 
#' @param variance A vector of observed variances
#' @param n A vector of degrees of fredom for each observed variance 
#' (number of observations - 1)
#' @export
#' @examples
#' x <- data.frame(true.var = 1 / rgamma(100000, 30, 3100))
#' x$obs.1 <- sapply(x$true.var, function(x) rnorm(1, 0, sqrt(x)))
#' x$obs.2 <- sapply(x$true.var, function(x) rnorm(1, 0, sqrt(x)))
#' x$obs.var <- apply(x[, c('obs.1', 'obs.2')], 1, var)
#' x$degf <- 1
#' x.ab <- estimateAB(x$obs.var, x$degf)
#' hist(x$true.var, freq=FALSE)
#' library(MCMCpack)
#' prediction <- dinvgamma(50:200, x.ab[1], x.ab[2])
#' lines(x=50:200, y=prediction)

estimateAB <- function(variance, n){
  
  strt <- c(2, 1)
  a <- nlminb(strt, flik, variance=variance, m=n, lower=c(0, 0))
  param <- a$par
  
  #Check if estimation of a and b was successful
  if (param[1] == strt[1] & param[2] == strt[2]) {
    stop('Initial values for a and b were returned.
         Estimation failed.')
  }
  
  #JAGS gamma distribution parameters for beta
  param[2] <- 1 / param[2]
  
  return (param)
  }

#' flik
#' 
#' This function calculates the log likelihood that the ratio of 
#' variance and the expected value of the fitted distribution come from an F
#' distribution with m and 2*a degrees of freedom. Where m is number of 
#' observations - 1, and a is the alpha parameter for the inverse gamma 
#' distribution.
#' 
#' @param param A list containing the values for a and b.
#' @param variance A vector of the observed variances
#' @param m A vector of degrees of fredom for each observed variance 
#' (number of observations - 1)
#' @export

flik <- function(param, variance, m){
  
  a <- param[1]
  b <- param[2]
  x <- variance*(a*b) # This is the "ratio" of the obs. and expected variance
  n <- 2*a
  out <- log(df(x, m, n))+log(a*b) #log(a*b) comes from the change in density
  # that comes from the change of variables from sig to a*b*sig
  
  return(sum(-out))
}