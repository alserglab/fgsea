#pragma once
#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
using namespace std;

// [[Rcpp::export]]
List calcGseaStatCumulativeBatch(
        NumericVector const& stats,
        double gseaParam,
        NumericVector const& pathwayScores,
        IntegerVector const& pathwaysSizes,
        int iterations,
        int seed,
        std::string scoreType);

//' Calculates GSEA statistic values for all the prefixes of a gene set
//'
//' Takes \emph{O(k^\{3/2\})} time, where \emph{k} is a size of `selectedSize`.
//' @param stats Named numeric vector with gene-level statistics
//'     sorted in decreasing order (order is not checked)
//' @param selectedStats indexes of selected genes in a `stats` array
//' @param gseaParam GSEA weight parameter (0 is unweighted, suggested value is 1)
//' @return Numeric vector of GSEA statistics for all prefixes of selectedStats.
//' @export
//' @examples
//' data(exampleRanks)
//' data(examplePathways)
//' ranks <- sort(exampleRanks, decreasing=TRUE)
//' es <- calcGseaStatCumulative(ranks, na.omit(match(examplePathways[[1]], names(ranks))), 1)
// [[Rcpp::export]]
NumericVector calcGseaStatCumulative(
        NumericVector const& stats,
        IntegerVector const& selectedStats, // Indexes start from one!
        double gseaParam,
        std::string scoreType = "std"
);

//' Calculates GSEA statistic valus for all gene sets in `selectedStats` list.
//'
//' Takes \emph{O(n + mKlogK)} time, where n is the number of genes, m is the number of gene sets,
//' and k is the mean gene set size.
//' @param stats Numeric vector of gene-level statistics sorted in decreasing order
//' @param selectedGenes List of integer vector with integer gene IDs (from 1 to n)
//' @param geneRanks Integer vector of gene ranks
//' @return Numeric vector of GSEA statistics of the same length as `selectedGenes` list
// [[Rcpp::export]]
NumericVector calcGseaStatBatchCpp(
        NumericVector const & stats,
        List const & selectedGenes,
        IntegerVector const & geneRanks // Ranks are 1-based!
);
