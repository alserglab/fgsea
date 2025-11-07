#ifndef FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H
#define FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H

#include <vector>
#include <random>
#include <set>
#include <cmath>
#include <algorithm>
#include <functional>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>

#include "esCalculation.h"
#include "util.h"
#include "Rcpp.h"

using namespace std;

class EsRuler {
private:
    using hash_t = uint64_t;
    using gsea_t = pair<score_t, hash_t>;

    struct Level {
        //  score, is_positive
        vector<pair<gsea_t, bool>> lowScores, highScores;
        gsea_t bound;
    };

    struct PerturbateResult {
        int moves;
        int iters;
    };

    const bool logStatus;
    // const vector<double> doubleRanks;
    const vector<int64_t> ranks;
    vector<hash_t> geneHashes;
    const unsigned int sampleSize;
    const unsigned int pathwaySize;
    const double movesScale;
    bool incorrectRuler = false;

    vector<vector<int>> currentSamples;
    int oldSamplesStart = 0;
    vector<Level> levels;
    // TODO: remove commented parts
    // vector<unsigned int> probCorrector;

    // void duplicateSamples();
    bool resampleGenesets(random_engine_t &rng);

    vector<int> chunkLastElement;
    int chunksNumber;

    struct SampleChunks {
        vector<int64_t> chunkSum;
        vector<vector<int>> chunks;
        SampleChunks(int);
    };

    PerturbateResult perturbate(vector<int64_t> const& ranks, int k, SampleChunks& sampleChunks,
               gsea_t bound, random_engine_t &rng);

    PerturbateResult perturbate_iters(vector<int64_t> const& ranks, int k, SampleChunks& sampleChunks,
               gsea_t bound, random_engine_t &rng, int iters);

    PerturbateResult perturbate_success(vector<int64_t> const& ranks, int k, SampleChunks& sampleChunks,
               gsea_t bound, random_engine_t &rng, int successes);

    PerturbateResult perturbate_until(vector<int64_t> const& ranks, int k, SampleChunks& sampleChunks, gsea_t bound, random_engine_t &rng, std::function<bool(int, int)> const& f);

    int chunkLen(int ind);

    hash_t calcHash(const vector<int>& curSample);

    vector<int64_t> scaleRanks(vector<double> const& ranks);

public:

    EsRuler(const vector<int64_t> &inpRanks,
            unsigned int inpSampleSize,
            unsigned int inpPathwaySize,
            double inpMovesScale,
            bool inpLog);

    ~EsRuler();

    void extend(double ES, int seed, double eps);

    // pair<double, bool> getPvalue(double ES, double eps, bool sign);
    tuple <double, bool, double> getPvalue(double ES, double eps, bool sign);
};

double betaMeanLog(unsigned long a, unsigned long b);
// double getErrorPerLevel(unsigned long k, unsigned long n);
double getVarPerLevel(unsigned long k, unsigned long n);

// TODO: do we still need it?
pair<double, bool> calcLogCorrection(const vector<unsigned int> &probCorrector,
                                     long probCorrIndx, unsigned int sampleSize);


#endif //FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H
