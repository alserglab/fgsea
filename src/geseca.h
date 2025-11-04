#ifndef GESECA_H
#define GESECA_H


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::export]]
Rcpp::List gesecaCpp(const NumericMatrix & E, const NumericVector & inpScores,
                    unsigned genesetSize, unsigned sampleSize, int seed, double eps);

#endif // GESECA_H
