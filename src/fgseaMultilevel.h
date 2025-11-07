#pragma once

#include <Rcpp.h>
using namespace Rcpp;

//' Calculates low GSEA p-values for a given gene set size using the multilevel split Monte Carlo approach.
//'
//' @param enrichmentScores A vector of enrichment scores, for which p-values should be calculated
//' @param ranks An integer vector with the gene-level statistics
//' @param pathwaySize A scalar with the size of the gene set
//' @param seed Random seed
//' @param eps P-values below eps aren't calculated
//' @param sign Controls whether ES^+ or ES score is used
//' @param moveScale Controls the number of MCMC iterations on each level
//' @param logStatus Controls whether debug output should be shown
//'
//' @keyword internal
//'
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::export]]
DataFrame fgseaMultilevelCpp(const NumericVector& enrichmentScores,
                             const SEXP& ranks,
                             int pathwaySize, int sampleSize,
                             int seed, double eps, bool sign,
                             double moveScale = 1.0, bool logStatus = false);

