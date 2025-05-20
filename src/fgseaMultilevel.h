#pragma once

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::export]]
DataFrame fgseaMultilevelCpp(const NumericVector& enrichmentScores,
                             const NumericVector& ranks,
                             int pathwaySize, int sampleSize,
                             int seed, double eps, bool sign,
                             double moveScale = 1.0, bool logStatus = false);

