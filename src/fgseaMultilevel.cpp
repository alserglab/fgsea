#include "fgseaMultilevel.h"
#include "fgseaMultilevelSupplement.h"
using namespace std;

DataFrame fgseaMultilevelCpp(const NumericVector& enrichmentScores,
                             const NumericVector& ranks, int pathwaySize,
                             int sampleSize, int seed,double eps,
                             bool sign, double moveScale, bool logStatus)
{
    vector<double> posRanks = as<std::vector<double> >(ranks);
    for (int i = 0; i < posRanks.size(); i++) {
        posRanks[i] = abs(posRanks[i]);
    }
    vector<double> negRanks = posRanks;
    reverse(negRanks.begin(), negRanks.end());

    const vector<double> esVector = as<std::vector<double> >(enrichmentScores);

    EsRuler esRulerPos(posRanks, sampleSize, pathwaySize, moveScale, logStatus);
    EsRuler esRulerNeg(negRanks, sampleSize, pathwaySize, moveScale, logStatus);

    double maxES = *max_element(begin(esVector), end(esVector));
    double minES = *min_element(begin(esVector), end(esVector));
    if (maxES >= 0){
        esRulerPos.extend(abs(maxES), seed, eps);
    }
    if (minES < 0){
        esRulerNeg.extend(abs(minES), seed, eps);
    }

    vector<double> pvalRes;
    vector<bool> isCpGeHalf;
    vector<double> log2err;

    int nrow = esVector.size();
    for (int i = 0; i < nrow; i++){
        // pair<double, bool> resPair;
        tuple <double, bool, double> res;
        double currentES = esVector[i];
        res = (currentES >= 0) ? esRulerPos.getPvalue(abs(currentES), eps, sign) : esRulerNeg.getPvalue(abs(currentES), eps, sign);

        pvalRes.push_back(get<0>(res));
        isCpGeHalf.push_back(get<1>(res));
        log2err.push_back(get<2>(res));
    }

    // return vector with pvalues and vector with conditional probability result
    return DataFrame::create(Named("cppMPval") = pvalRes,
                             Named("cppIsCpGeHalf") = isCpGeHalf,
                             Named("log2err") = log2err);
}
