#ifndef FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H
#define FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H

#include <vector>
#include <random>
#include <set>
#include <cmath>
#include <algorithm>

#include "util.h"


using namespace std;

class EsRuler {
private:

    const vector<double> &ranks;
    const unsigned int sampleSize;
    const unsigned int pathwaySize;


    vector<double> enrichmentScores;
    vector<vector<int> > currentSamples;

    vector<unsigned int> probCorrector;

    void duplicateSamples();

    vector<int> chunkLastElement;
    int chunksNumber;

    struct SampleChunks {
        vector<double> chunkSum;
        vector<vector<int>> chunks;
        SampleChunks(int);
    };

    int perturbate(const vector<double> &ranks, int k, SampleChunks &cusSampleChunks,
               double bound, random_engine_t& rng);

    int chunkLen(int ind);

public:

    EsRuler(const vector<double> &inpRanks, unsigned int inpSampleSize, unsigned int inpPathwaySize);

    ~EsRuler();

    void extend(double ES, int seed, double eps);

    pair<double, bool> getPvalue(double ES, double eps, bool sign);
};

pair<double, bool> calcLogCorrection(const vector<unsigned int> &probCorrector,
                                     long probCorrIndx, unsigned int sampleSize);


#endif //FGSEAMULTILEVELCPP_FGSEAMULTILEVELSUPPLEMENT_H
