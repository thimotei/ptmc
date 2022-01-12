#' @useDynLib ptmc
#' @importFrom Rcpp sourceCpp
#' @import coda
#' @import parallel
#' @import tidyr
#' @import dplyr
#' @import foreach
#' @import doParallel
#' @importFrom magrittr %>% %<>%
NULL

#' Create a model
#'
#' @param model The model to run ptmc.
#' @param data Data used in the calibration process
#' @param settings settings
#' @return Returns a list with the fist element being the mcmc samples formatting for analysis and plottig with the CODA package. The second is the log posterior value at each time step.
#'
#' @export
ptmc_discrete_func <- function(model, data_list, settings, par = NULL) {
  if (length(par) == 0) {
    par <- rep(list(list(type = "None")), settings[["numberChainRuns"]])
    output <- get_discrete_outputB(model, data_list, settings, FALSE, par)
  } else {
    output <- get_discrete_outputB(model, data_list, settings, TRUE, par)
  }
  output
}

get_discrete_outputB <- function(model, data_list, settings, update_ind, par) {
  nCores <- detectCores() - 1
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  outPTpost <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTdiscrete <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTlp <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTtemp <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTacc <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTpar <- vector(mode = "list", length = settings[["numberChainRuns"]])

  # Run the chains in parallel
  out_raw <- list()
  if (settings[["runParallel"]]) {
    out_raw <- foreach(i = 1:settings[["numberChainRuns"]], .packages = c('ptmc','coda')) %dopar% {
      run_ptmc_discrete(model, data_list, settings, update_ind, par[[i]])
    }
  } else {
    for (i in 1:settings[["numberChainRuns"]]) {
      out_raw[[i]] <- run_ptmc_discrete(model, data_list, settings, update_ind, par[[i]])
    }
  }

  for(i in 1:settings[["numberChainRuns"]]) {
    out_post <- out_raw[[i]][["output"]][, 1:settings$numberFittedPar]
    outPTpar[[i]] <- out_raw[[i]][["PTMCpar"]]
    if (settings$numberFittedPar > 1){
        colnames(out_post) <- model[["namesOfParameters"]]
    }
    outPTpost[[i]] <- mcmc(out_post)
    outPTdiscrete[[i]] <- out_raw[[i]][["discrete"]]
  #  cat(out_raw[[i]][["discrete"]])
    outPTlp[[i]] <- out_raw[[i]][["output"]][, settings$numberFittedPar + 1]
    outPTtemp[[i]] <- out_raw[[i]][["output"]][, settings$numberFittedPar + 2]
    outPTacc[[i]] <- out_raw[[i]][["output"]][, settings$numberFittedPar + 3]
  }

  outlpv <- data.frame(matrix(unlist(outPTlp), nrow = length(outPTlp[[1]])))
  colnames(outlpv) <- c(1:settings[["numberChainRuns"]])
  outlpv <- outlpv %>% gather(colnames(outlpv), key="chain_no",value="lpost")
  outlpv$sample_no <-rep(1:length(outPTlp[[1]]), settings[["numberChainRuns"]])

 # outdiscretev <- data.frame(matrix(unlist(outPTdiscrete), nrow = length(outPTdiscrete[[1]])))
 # colnames(outdiscretev) <- c(1:settings[["numberChainRuns"]])
  #outdiscretev <- outdiscretev %>% gather(colnames(outdiscretev), key="chain_no",value="lpost")
  #outdiscretev$sample_no <- rep(1:length(outdiscretev[[1]]), settings[["numberChainRuns"]])

  outltempv <- data.frame(matrix(unlist(outPTtemp), nrow=length(outPTtemp[[1]])))
  colnames(outltempv) <- c(1:settings[["numberChainRuns"]])
  outltempv <- outltempv %>% gather(colnames(outltempv), key="chain_no", value="temperature")
  outltempv$sample_no <- rep(1:length(outPTtemp[[1]]), settings[["numberChainRuns"]])
  
  outlaccv <- data.frame(matrix(unlist(outPTacc), nrow=length(outPTacc[[1]])))
  colnames(outlaccv) <- c(1:settings[["numberChainRuns"]])
  outlaccv <- outlaccv %>% gather(colnames(outlaccv), key="chain_no", value="acceptance rate")
  outlaccv$sample_no <- rep(1:length(outPTacc[[1]]), settings[["numberChainRuns"]])

  output <- list(
    mcmc = as.mcmc.list(outPTpost),
    discrete = outPTdiscrete,
    lpost = outlpv,
    temp = outltempv,
    acc = outlaccv,
    outPTpar = outPTpar
  )
  output
}

get_outputA <- function(model, settings, update_ind, par) {
  outPTpost <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTlp <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTtemp <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTacc <- vector(mode = "list", length = settings[["numberChainRuns"]])

  for (i in 1:settings[["numberChainRuns"]]) {
    out_raw <- run_ptmc(model, settings, update_ind, par[[i]])
    out_post <- out_raw[, 1:settings$numberFittedPar]
    if (settings$numberFittedPar > 1){
        colnames(out_post) <- model[["namesOfParameters"]]
    }
    outPTpost[[i]] <- mcmc(out_post)
    
    outPTlp[[i]] <- out_raw[, settings$numberFittedPar + 1]
    outPTtemp[[i]] <- out_raw[, settings$numberFittedPar + 2]
    outPTacc[[i]] <- out_raw[, settings$numberFittedPar + 3]
  }

  outlpv <- data.frame(matrix(unlist(outPTlp), nrow=length(outPTlp[[1]])))
  colnames(outlpv) <- c(1:settings[["numberChainRuns"]])
  outlpv <- outlpv %>% gather(colnames(outlpv), key="chain_no",value="lpost")
  outlpv$sample_no <-rep(1:length(outPTlp[[1]]), settings[["numberChainRuns"]])

  outltempv <- data.frame(matrix(unlist(outPTtemp), nrow=length(outPTtemp[[1]])))
  colnames(outltempv) <- c(1:settings[["numberChainRuns"]])
  outltempv <- outltempv %>% gather(colnames(outltempv), key="chain_no", value="temperature")
  outltempv$sample_no <- rep(1:length(outPTtemp[[1]]), settings[["numberChainRuns"]])
  
  outlaccv <- data.frame(matrix(unlist(outPTacc), nrow=length(outPTacc[[1]])))
  colnames(outlaccv) <- c(1:settings[["numberChainRuns"]])
  outlaccv <- outlaccv %>% gather(colnames(outlaccv), key="chain_no", value="acceptance rate")
  outlaccv$sample_no <- rep(1:length(outPTacc[[1]]), settings[["numberChainRuns"]])

  output <- list(
    mcmc = as.mcmc.list(outPTpost),
    lpost = outlpv,
    temp = outltempv,
    acc = outlaccv
  )
  output
}
