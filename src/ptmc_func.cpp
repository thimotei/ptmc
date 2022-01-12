#include <Rcpp.h>
#include <RcppEigen.h>

#include "./headers/mvn.hpp"
#include "./headers/ptmc.hpp"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]

using RPTMC = ptmc::PTMC;

void init_samplePriorDistributions(RPTMC* model, Rcpp::Function samplePriorDistributions) {
  auto func = [samplePriorDistributions]() {
    PutRNGstate();
    auto rData = samplePriorDistributions();
    GetRNGstate();
    return Rcpp::as<VectorXd>(rData);
  };
  model->samplePriorDistributions = func;
}

void init_evaluateLogPrior(RPTMC* model, Rcpp::Function evaluateLogPrior) {
  auto func = [evaluateLogPrior](VectorXd params) {
    PutRNGstate();
    auto rData = evaluateLogPrior(params);
    GetRNGstate();
    return Rcpp::as<double>(rData);
  };
  model->evaluateLogPrior = func;
}

void init_evaluateLogLikelihood(RPTMC* model, Rcpp::Function evaluateLogLikelihood) {
  auto func = [evaluateLogLikelihood](VectorXd params, MatrixXd covariance, List dataList) {
    PutRNGstate();
    auto rData = evaluateLogLikelihood(params, covariance, dataList);
    GetRNGstate();
    return Rcpp::as<double>(rData);
  };
  model->evaluateLogLikelihood = func;
}


// [[Rcpp::export]]
List run_ptmc(Rcpp::List model, Rcpp::List dataList, Rcpp::List settings, bool update_ind, Rcpp::List PTMCpar)
{
  RPTMC PTMC; MatrixXd output;
  init_samplePriorDistributions(&PTMC, model["samplePriorDistributions"]);
  init_evaluateLogPrior(&PTMC, model["evaluateLogPrior"]);
  init_evaluateLogLikelihood(&PTMC, model["evaluateLogLikelihood"]);

  if (update_ind) {
    PTMC.updateClass(settings, dataList, PTMCpar);
  } else {  
    PTMC.initialiseClass(settings, dataList);
  }

  output = PTMC.runPTMCC();
  PTMCpar = PTMC.savePTMCpar();
  return Rcpp::List::create(_["output"] = output, _["PTMCpar"] = PTMCpar);

}
