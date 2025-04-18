#ifndef SCORE_RULER_H
#define SCORE_RULER_H

#include <vector>
#include <random>
#include <algorithm>

#include "util.h"
#include "ScoreCalculation.h"



class ScoreRuler {
private:
    const unsigned n;
    const unsigned m;
    std::vector<float> expressionMatrix;
    const unsigned sampleSize;
    const unsigned genesetSize;

    const double upPrmtr = 0.2;
    const unsigned itersPerStep;

    std::vector<double> scores;
    std::vector<std::vector<unsigned> > currentSample;
    std::vector<std::vector<float> > currentProfiles;


    void duplicateSampleElements();
    int updateElement(std::vector<unsigned> & element, std::vector<float> & profile,
                      double threshold, random_engine_t & mtGen);
public:
    ScoreRuler(const std::vector<std::vector<float> > & inpE,
               unsigned inpSampleSize, unsigned inpGenesetSize);
    ~ScoreRuler();

    void extend(double inpScore, int seed, double eps);
    std::pair<double, double> getPvalue(double inpScore, double eps);

    // pair<double, bool> getPvalue(double ES, double eps, bool sign);
};

#endif //SCORE_RULER_H
