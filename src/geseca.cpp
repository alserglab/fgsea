#include "geseca.h"
#include "ScoreRuler.h"


Rcpp::List gesecaCpp(const NumericMatrix & E, const NumericVector & inpScores,
                        unsigned genesetSize, unsigned sampleSize, int seed, double eps){
    std::vector<double> scores = as<std::vector<double> > (inpScores);
    std::vector<std::vector<float> > expressionMatrix(E.nrow());
    for (unsigned i = 0; i < E.nrow(); i++){
        NumericVector currentRow = E(i, _);
        expressionMatrix[i] = as<std::vector<float> >(currentRow);
    }

    ScoreRuler ruler(expressionMatrix, sampleSize, genesetSize);
    std::vector<double> pvals;

    double maxScore = *std::max_element(scores.begin(), scores.end());
    ruler.extend(maxScore, seed, eps);

    Rcpp::List res;

    for (auto score : scores){
        std::pair<double, double> pval = ruler.getPvalue(score, eps);
        res.push_back(Rcpp::List::create(Rcpp::Named("pval") = pval.first,
                                         Rcpp::Named("log2err") = pval.second));
    }

    return res;
}
