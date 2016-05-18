"
Pipeline for Analysis of Primary High Throughput Screens
Benjamin Chittick

This script contains functions for logistic curve fitting
"

#library dependencies
library(data.table)

#' curveFit
#' 
#' When the only tool you have is a hammer, every problem looks like a nail. 
#' curveFit, fits a 4 parameter logistic inhibition curve according to
#' y[i] ~ N(mu[i], tau)
#' mu <- top + (bottom - top)/(1 + 10^((log(IC50) - log(x))*hillslope))
#' 
#' IC50 <- exp(1/100)
#' hillslope ~ dunif(0.1, 3)
#' top ~ dnorm(0, 0.01)
#' bottom ~ dnorm(-100, 0.01)
#' tau <- pow(sigma, -2)
#' sigma ~ dunif(0, 100)
#' 
#' @param x Numeric of compound concentrations
#' @param y Numeric measured values aligned with x
#' @export
#' @examples
#' x <- c(2, 4, 8, 16, 32, 64, 128, 256)
#' y <- -100/(1 + 10^((log(27) - log(x))*1.23))
#' y.error <- y + rnorm(length(y), 0, 5)
#' plot(y.error, x)
#' estimate <- curveFit(x, y.error)
#' estimate

curveFit <- function(x, y, top.mu, bottom.mu,
                     runjags.method='rjparallel'){
  
  #Make sure library dependencies are loaded
  library(coda)
  library(rjags)
  library(runjags)
  
  # The model for sampling the data posterior mean and variance
  model.str <- 'model {
    for (i in 1:N) {
      y[i] ~ dnorm(mu[i], tau)
      mu[i] <- top + (bottom - top)/(1+10^((log(IC) - log(x[i]))*hillslope))
    }
    top ~ dnorm(top.mu, 0.001)
    bottom ~ dnorm(bottom.mu, 0.001)
    IC ~ dexp(1/1000)
    hillslope ~ dunif(0.1, 3)
    sigma ~ dunif(0, 100)
    tau <- pow(sigma, -2)
  }'
  
  # JAGS uses precision instead of standard deviation.
  # Precision: tau = 1 / variance
  inits1 <- list(.RNG.name='base::Super-Duper', .RNG.seed=7)
  inits2 <- list(.RNG.name='base::Wichmann-Hill', .RNG.seed=13)
  data.model <- try(autorun.jags(model=model.str,
                             data=list(y=y,
                                       x=x,
                                       N=length(y),
                                       top.mu=top.mu,
                                       bottom.mu=bottom.mu),
                             monitor=c('IC', 'hillslope', 'top', 'bottom'),
                             inits=list(inits1,
                                        inits2),
                             method=runjags.method))
  
  chain <- data.frame(as.mcmc(data.model))
  IC50 <- mean(chain$IC)
  IC50.error <- sd(chain$IC)
  hillslope <- mean(chain$hillslope)
  hillslope.error <- sd(chain$hillslope)
  top <- mean(chain$top)
  bottom <- mean(chain$bottom)
  
  #catch errors in data.model
  if (class(data.model) == 'try-error') {
    IC50 = NA
    IC50.error = NA
    hillslope = NA
    hillslope.error = NA
    top <- NA
    bottom <- NA
  }
  
  estimate <- data.frame('IC50'=IC50,
                         'IC50_Error'=IC50.error,
                         'Hillslope'=hillslope,
                         'Hillslope_Error'=hillslope.error,
                         'Top'=top,
                         'Bottom'=bottom)
  
  return(estimate)
}

#' 
#' 
#' Take HTS confirmation in dose and fit curves using curveFit
#' 
#' @param screen data.frame or data.table of confirmation in dose screen
#' @export

screenFit <- function(screen){
   
  if(!all(c('Concentration','Compound_ID', 'Value') %in% names(screen))){
    stop('Screen data must contain columns: Compound_ID, Concentration, and
          Value')
  }
  
  #make sure snow is not attached
  if (require(snow)) detach("package:snow", unload=TRUE)
  
  bottom <- screen[screen$Well_Type == 'positive control', ]
  bottom.mu <- median(bottom$Value)
  
  top <- screen[screen$Well_Type == 'neutral control', ]
  top.mu <- median(top$Value)
  
  cpd <- screen[screen$Well_Type == 'compound', ]
  cpd <- cpd[cpd$Compound_ID != '', ]

  by.compounds <- split(cpd, cpd$Compound_ID)
  curve.fits <- lapply(by.compounds, 
                       function(x) curveFit(x$Concentration, 
                                            x$Value, 
                                            top.mu=top.mu,
                                            bottom.mu=bottom.mu))
  
  compound.id <- names(curve.fits)
  curve.fits <- rbindlist(curve.fits)
  curve.fits$Compound_ID <- compound.id
  
  screen <- merge(screen, curve.fits, by='Compound_ID', all.x=TRUE)
  
  #reattach snow
  require(snow)
  
  return(screen)
}  
  