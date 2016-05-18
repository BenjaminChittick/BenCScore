"
Pipeline for Analysis of Primary High Throughput Screens
Benjamin Chittick

This script contains functions for sampling Bayesian models of compound activity
in high-throughput screens

References:
  http://www.indiana.edu/~kruschke/BEST/BEST.pdf
"

#library dependencies
library(coda)
library(rjags)
library(runjags)

#' tSample
#' 
#' Samples the posterior distribution of a compound's activity according to: \cr
#' \cr
#' data ~ dt(mu, tau, nu) \cr
#' mu ~ dt(mu.prior, tau.prior, nu.prior) \cr
#' tau ~ dgamma(a, b) \cr
#' nu ~ dexp(lambda) \cr 
#' \cr
#' The sampling is done by MCMC via autorun.jags.
#' 
#' @param data list of compound activity data
#' @param mu.prior location parameter for prior distribution of mu, typically the average of 
#' all compound activity.
#' @param var.prior spread parameter for prior distribution of mu, typically the square 
#' standard error of the mean for all compounds. 
#' @param a alpha parameter for prior distribution of compound variance. See EstimateAB.
#' @param b beta parameter for prior distribution of compound variance. See EstimateAB.
#' @param lambda paramater for normality of t-distribution. Defaults to 1/29. 
#' Based on the Bayesian Estimation Superseeds T-test http://www.indiana.edu/~kruschke/BEST/BEST.pdf.
#' @param runjags.method method by which jags is run, default is rjparallel
#' @export
#' @examples
#' compound <- rnorm(2, 100, 10)
#' mu.prior <- 0
#' var.prior <- 100
#' nu.prior <- 1
#' a <- 30
#' b <- 3100
#' posterior <- tSample(compound, mu.prior, var.prior, nu.prior, a, b)
#' hist(posterior$mu, 100)
#' hist(posterior$tau, 100)

tSample <- function(data, mu.prior, var.prior, nu.prior, a, b, lambda=1/29,
                    runjags.method='rjparallel'){
  #Make sure library dependencies are loaded
  library(coda)
  library(rjags)
  library(runjags)
  
  # The model for sampling the data posterior mean and variance
  model.str <- 'model {
    for (i in 1:N) {
      dat[i] ~ dt(mu, tau, nu)
    }
    mu ~ dt(mu.prior, tau.prior, nu.prior)
    tau ~ dgamma(a, b)
    nuPrime ~ dexp(lambda)
    nu <- nuPrime + 1
  }'
  
  # JAGS uses precision instead of standard deviation.
  # Precision: tau = 1 / variance
  inits1 <- list(.RNG.name='base::Super-Duper', .RNG.seed=7)
  inits2 <- list(.RNG.name='base::Wichmann-Hill', .RNG.seed=13)
  data.model <- autorun.jags(model=model.str,
                             data=list(dat=data,
                                       N=length(data),
                                       mu.prior=mu.prior,
                                       tau.prior=1/var.prior,
                                       nu.prior=nu.prior,
                                       lambda=lambda,
                                       a=a, 
                                       b=b),
                             monitor=c('mu', 'tau'),
                             inits=list(inits1,
                                        inits2),
                             method=runjags.method)
  
  chain <- data.frame(as.mcmc(data.model))

  return (chain)
}

#' tContrast
#' 
#' Wraps the tSample function and compares its markov chain to a provided null
#' markov chain. Returns estimated mean, standard deviation, and probability of 
#' belonging to the null posterior.
#' 
#' @param data activity data for a compound
#' @param mu.prior location parameter for prior distribution of mu, typically the average of 
#' all compound activity.
#' @param var.prior spread parameter for prior distribution of mu, typically the square 
#' standard error of the mean for all compounds. 
#' @param a alpha parameter for prior distribution of compound variance. See EstimateAB.
#' @param b beta parameter for prior distribution of compound variance. See EstimateAB.
#' @param null.post a MCMC sampling of a null posterior with column name mu.
#' @param alternative a character string specifying the relationship to the null posterior, 
#' must be one of "two.sided", "greater", or "less". 
#' @param rope Region of Practical Equivalence. This value specifies how close 
#' the sampled posterior needs to be to the null posterior to consider them 
#' equivalent, default value is 0.
#' @param runjags.method method by which jags is run, default is rjparallel
#' @export
#' @examples
#' compound <- rnorm(2, 100, 10)
#' null.reference <- rnorm(320, 0, 10)
#' pooled <- c(compound, null.reference)
#' mu.prior <- median(pooled)
#' var.prior <- var(pooled)
#' nu.prior <- 1
#' a <- 30
#' b <- 3100
#' ROPE <- 15
#' null.post <- tSample(null.reference, mu.prior, var.prior, nu.prior, a, b)
#' posterior <- tContrast(compound, mu.prior, var.prior, nu.prior, a, b, 
#'                        null.post, ROPE, alternative='greater')
#' posterior

tContrast <- function(data, mu.prior, var.prior, nu.prior, a, b, null.post, 
                      ROPE, alternative, runjags.method='rjparallel'){
  
  posterior <-  tSample(data=data,
                        mu.prior=mu.prior,
                        var.prior=var.prior,
                        nu.prior=nu.prior,
                        a=a, b=b,
                        runjags.method=runjags.method)
  
  estimate <- data.frame('Mu'=mean(posterior$mu),
                         'Tau'=mean(posterior$tau))

  delta <- posterior$mu - null.post$mu
  if (alternative == 'greater') {
    estimate$Alternate_Probability_t <- mean(delta > ROPE) - 1/dim(posterior)[1]
  }
  if (alternative == 'less'){
    estimate$Alternate_Probability_t <- mean(delta < ROPE) - 1/dim(posterior)[1]
  }
  
  if(estimate$Alternate_Probability_t < 1/dim(posterior)[1]){
    estimate$Alternate_Probability_t <- 1/dim(posterior)[1]
  } 
  
  return(estimate)
}

#' tParallel
#' 
#' This function runs tContrast sampling in parallel
#' 
#' @param data compound activity data with individual compounds in rows and 
#' replicates in columns
#' @param mu.prior location parameter for prior distribution of mu. Typically the average of 
#' all compound activity.
#' @param var.prior spread parameter for prior distribution of mu. Typically the square 
#' standard error of the mean for all compounds. 
#' @param a alpha parameter for prior distribution of compound variance. See EstimateAB.
#' @param b beta parameter for prior distribution of compound variance. See EstimateAB.
#' @param null.post a MCMC sampling of a null posterior with column name mu.
#' @param alternative a character string specifying the relationship to the null posterior, 
#' must be one of "two.sided", "greater", or "less". 
#' @param rope Region of Practical Equivalence. This value specifies how close 
#' the sampled posterior needs to be to the null posterior to consider them 
#' equivalent, default value is 0.
#' @param runjags.method method by which jags is run, default is set to 'simple' to 
#' avoid conflict with snow library
#' @export
#' @examples
#' set.seed(7)
#' compounds <- matrix(rnorm(20, seq(0, 100, 10), 3), nrow=10)
#' null.reference <- rnorm(32, 0, 3)
#' pooled <- c(null.reference, sapply(compounds, function(x) x))
#' mu.prior <- median(pooled)
#' var.prior <- var(null.reference)
#' nu.prior <- 1
#' a <- 30
#' b <- 3100
#' ROPE <- 15
#' 
#' null.post <- tSample(null.reference, mu.prior, var.prior, nu.prior, a, b)
#' posteriors <- tParallel(compounds, mu.prior, var.prior, nu.prior, a, b,
#'                         null.post, ROPE, alternative='greater')
#' rbindlist(posteriors)
#' compounds

tParallel <- function(data, mu.prior, var.prior, nu.prior, a, b, null.post, 
                      ROPE, alternative, runjags.method='simple'){
  
  # snow does not send the variables to the workers unless they are referenced
  # before hand
  mu.prior; var.prior; nu.prior; a; b; null.post; alternative; ROPE
  runjags.method
  null.post <- null.post
  
  #runjags.method='rjags'
  #rjags is a better method but sometimes the
  #Gelman-Rubin statistic will not converge
  
  library(snow) # snow conflicts with runjags method rjparallel
  n.cores <- parallel::detectCores()
  
  clus <- makeCluster(n.cores)
  
  clusterExport(clus, c('mu.prior', 
                        'var.prior',
                        'nu.prior',
                        'a', 
                        'b', 
                        'null.post', 
                        'alternative',
                        'ROPE',
                        'runjags.method'),
                envir=environment())
  
  clusterExport(clus, 'tContrast')
  clusterExport(clus, 'tSample')
  clusterExport(clus, 'autorun.jags')
  clusterExport(clus, 'as.mcmc')
  
  estimate <- parApply(clus, data, 1,
                       function(x) tContrast(data=x, mu.prior=mu.prior,
                                             var.prior=var.prior,
                                             nu.prior=nu.prior,
                                             a=a, b=b, null.post=null.post,
                                             alternative=alternative,
                                             ROPE=ROPE,
                                             runjags.method=runjags.method))
  
  detach(package:snow, unload=TRUE)
  
  return(estimate)
}

#' nSample
#' 
#' Depricated: sample a posterior distribution of the form
#'   mu|x,tau ~ N(mu, tau)
#' where mu and tau are unknown
#' 
#' @param x an array of data
#' @param a alpha paramter
#' @param b beta parameter
#' @param mu0 mu ~ N(mu0, n0)
#' @param n0 mu ~ N(mu0, n0)
#' @export
#' @examples
#' compound <- rnorm(2, 100, 10)
#' a <- 30 
#' b <- 3100
#' mu0 <- 0.0
#' n0 <- 0.0001
#' size <- 20000
#' 
#' posterior <- nSample(compound, a, b, mu0, n0, size)
#' hist(posterior$mu)

nSample <- function(x, a, b, mu0, n0, size){
  #sample a posterior distribution of the form
  # mu|x,tau ~ N(mu, tau)
  #where mu and tau are unknown
  
  n <- sum(!is.na(x))
  x.bar <- median(x, na.rm=TRUE)
  
  #tau posterior parameters
  a.prime <- a + n/2
  b.prime <- b + 1/2*sum((x-x.bar)^2, na.rm=TRUE) + n*n0/(2*(n + n0))*(x.bar - mu0)^2
  
  #sample tau|x
  tau <- rgamma(size, a.prime, b.prime)
  
  #mu posterior parameters
  mu.prime <- n*tau/(n*tau + n0*tau)*x.bar + n0*tau/(n*tau + n0*tau)*mu0
  tau.prime <- n*tau + n0*tau
  
  #sample mu|x,tau
  mu <- rnorm(size, mu.prime, sqrt(1/tau.prime))
  
  chain <- data.frame('mu'=mu, 'tau'=tau)
  
  return(chain)
}

#' nContrast
#' 
#' Depricated: estimate the probability a compound is from the neutral
#' distribution
#' 
#' @param x an array of data
#' @param a alpha paramter
#' @param b beta parameter
#' @param mu0 mu ~ N(mu0, n0)
#' @param n0 mu ~ N(mu0, n0)
#' @param null.post data.frame or data.table of a null posterior with column mu
#' @export
#' @examples
#' set.seed(7)
#' compound <- rnorm(2, 100, 10)
#' null.reference <- rnorm(32, 0, 10)
#' a <- 30
#' b <- 3100
#' mu0 <- 0.0
#' n0 <- 0.0001
#' size <- 20000
#'
#' null.post <- nSample(null.reference, a, b, mu0, n0, size)
#' posterior <- nContrast(compound, a, b, mu0, n0, null.post, 
#'                        alternative='greater', size)
#' posterior

nContrast <- function(x, a, b, mu0, n0, null.post, alternative, size, ROPE){
  
  posterior <- nSample(x, a, b, mu0, n0, size)
  
  estimate <- data.frame('Mu_n'=mean(posterior$mu),
                         'Tau_n'=mean(posterior$tau))
  
  delta <- posterior$mu - null.post$mu
  if (alternative == 'greater') {
    estimate$Alternate_Probability_nDist <- mean(delta > ROPE)
  }
  if (alternative == 'less'){
    estimate$Alternate_Probability_nDist <- mean(delta < ROPE)
  }
  
  return(estimate)
}