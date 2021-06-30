library(Rcpp)       # For c++ intergration
#library(ptmc)

# ptmc solver
model <- list(
  
  namesOfParameters =  c("mu"),
  
  # Generate the initial step of Markov chain
  samplePriorDistributions = function() {
    s1 <- runif(1, min=0, max=1)
    s1
  },
  
  # Evaluate the log prior
  evaluateLogPrior = function(params) {
    if (params < 0)
      return(log(0))
    else if (params > 1)
      return(log(0))
    else 
      p1 <- dunif(params, min=0, max=1, log=TRUE)
    
    p1
  },
  
  # Evaluate the log likelihood
  evaluateLogLikelihood = function(params, covar) {
    nAA <- 50
    nAa <- 21
    naa <- 29
    ll <- (2*nAA)*log(params) + nAa*log(2*params*(1-params)) + 2*naa*log(1-params)
    ll
  }
)

settingsPT <-  list(
  numberChainRuns = 4,
  numberTempChains = 5,
  iterations = 1000, 
  burninPosterior = 1,
  thin = 10,
  consoleUpdates = 1,
  numberFittedPar = 1,        
  onAdaptiveCov = FALSE,
  updatesAdaptiveCov = 1,
  burninAdaptiveCov = 100,
  onAdaptiveTemp = TRUE,
  updatesAdaptiveTemp = 1,
  onDebug = FALSE,
  lowerParBounds = c(-100000.0),
  upperParBounds = c(100000.0)
)

data <- list(
  obsdata=0,
  poly=0,
  dem=0,
  agegroup=0
)

#outPT <- ptmc_func(model, settingsPT)

