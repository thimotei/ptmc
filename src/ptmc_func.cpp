#include <Rcpp.h>
#include <RcppEigen.h>

#include "./headers/mvn.hpp"
#include "./headers/ptmc.hpp"

#include "./headers/ptmc_discrete.hpp"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]

////////////////////////////////////////
////////////////////////////////////////
//////////// CONTINUOUS PTMC //////////
///////////////////////////////////
///////////////////////////////////

using RPTMC = ptmc::PTMC;

using RPTMC_D = ptmc_discrete::PTMC_D;

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



////////////////////////////////////////
////////////////////////////////////////
//////////// DISCRETE PTMC //////////
///////////////////////////////////
///////////////////////////////////



void init_samplePriorDistributions_discrete(RPTMC_D* model, Rcpp::Function samplePriorDistributions) {
  auto func = [samplePriorDistributions]() {
    PutRNGstate();
    auto rData = samplePriorDistributions();
    GetRNGstate();
    return Rcpp::as<VectorXd>(rData);
  };
  model->samplePriorDistributions = func;
}

void init_evaluateLogPrior_discrete(RPTMC_D* model, Rcpp::Function evaluateLogPrior) {
  auto func = [evaluateLogPrior](VectorXd params, VectorXi discrete) {
    PutRNGstate();
    auto rData = evaluateLogPrior(params, discrete);
    GetRNGstate();
    return Rcpp::as<double>(rData);
  };
  model->evaluateLogPrior = func;
}

void init_evaluateLogLikelihood_discrete(RPTMC_D* model, Rcpp::Function evaluateLogLikelihood) {
  auto func = [evaluateLogLikelihood](VectorXd params, VectorXi discrete, MatrixXd covariance, List dataList) {
    PutRNGstate();
    auto rData = evaluateLogLikelihood(params, discrete, covariance, dataList);
    GetRNGstate();
    return Rcpp::as<double>(rData);
  };
  model->evaluateLogLikelihood = func;
}

void init_initialiseDiscrete_discrete(RPTMC_D* model, Rcpp::Function initialiseDiscrete) {
  auto func = [initialiseDiscrete]() {
    PutRNGstate();
    auto rData = initialiseDiscrete();
    GetRNGstate();
    return Rcpp::as<VectorXi>(rData);
  };
  model->initialiseDiscrete = func;
}

void init_discreteSampling_discrete(RPTMC_D* model, Rcpp::Function discreteSampling) {
  auto func = [discreteSampling](VectorXi discrete) {
    PutRNGstate();
    auto rData = discreteSampling(discrete);
    GetRNGstate();
    return Rcpp::as<VectorXi>(rData);
  };
  model->discreteSampling = func;
}


// [[Rcpp::export]]
List run_ptmc_discrete(Rcpp::List model, Rcpp::List dataList, Rcpp::List settings, bool update_ind, Rcpp::List PTMCpar)
{
  RPTMC_D PTMC; List output_full;
  MatrixXd output;
  MatrixXi discrete;
  init_samplePriorDistributions_discrete(&PTMC, model["samplePriorDistributions"]);
  init_evaluateLogPrior_discrete(&PTMC, model["evaluateLogPrior"]);
  init_evaluateLogLikelihood_discrete(&PTMC, model["evaluateLogLikelihood"]);
  init_initialiseDiscrete_discrete(&PTMC, model["initialiseDiscrete"]);
  init_discreteSampling_discrete(&PTMC, model["discreteSampling"]);

  if (update_ind) {
    PTMC.updateClass(settings, dataList, PTMCpar);
  } else {  
    PTMC.initialiseClass(settings, dataList);
  }


  output_full = PTMC.runPTMCC();
  PTMCpar = PTMC.savePTMCpar();
  output = output_full[0];
  discrete = output_full[1];
  return Rcpp::List::create(_["output"] = output, _["discrete"] = discrete, _["PTMCpar"] = PTMCpar);

}
