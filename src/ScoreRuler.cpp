#include "ScoreRuler.h"
#include "util.h"
#include "Rcpp.h"

using namespace Rcpp;

ScoreRuler::ScoreRuler(const std::vector<std::vector<float> > & inpE,
                       unsigned inpSampleSize, unsigned inpGenesetSize):
        n(inpE.size()), m(inpE[0].size()),
        sampleSize(inpSampleSize),
        genesetSize(inpGenesetSize),
        itersPerStep(std::max(unsigned(1), unsigned(inpGenesetSize * upPrmtr))) {
    currentSample.resize(inpSampleSize);
    currentProfiles.resize(inpSampleSize);

    expressionMatrix = std::vector<float>(n*m);
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j) {
            expressionMatrix[i*m + j] = inpE[i][j];
        }
    }
}



ScoreRuler::~ScoreRuler() = default;


void ScoreRuler::duplicateSampleElements(){
    /*
     * Replaces sample elements that are less than median with elements
     * that are greater that the median
     */

    std::vector<std::pair<double, unsigned> > scoreAndIndex(sampleSize);

    for (unsigned elemIndex = 0; elemIndex < sampleSize; elemIndex++) {
        double elemScore = getScore(currentProfiles[elemIndex]);
        scoreAndIndex[elemIndex] = std::make_pair(elemScore, elemIndex);
    }
    std::sort(scoreAndIndex.begin(), scoreAndIndex.end());

    for (unsigned elemIndex = 0; 2 * elemIndex < sampleSize; elemIndex++){
        scores.push_back(scoreAndIndex[elemIndex].first);
    }

    std::vector<std::vector<unsigned> > tempSample;
    std::vector<std::vector<float> > tempProfiles;
    for (unsigned elemIndex = 0; 2 * elemIndex < sampleSize - 2; elemIndex++) {
        for (unsigned rep = 0; rep < 2; rep++) {
            tempSample.push_back(currentSample[scoreAndIndex[sampleSize - 1 - elemIndex].second]);
            tempProfiles.push_back(currentProfiles[scoreAndIndex[sampleSize - 1 - elemIndex].second]);
        }
    }
    tempSample.push_back(currentSample[scoreAndIndex[sampleSize >> 1].second]);
    tempProfiles.push_back(currentProfiles[scoreAndIndex[sampleSize >> 1].second]);
    std::swap(currentSample, tempSample);
    std::swap(currentProfiles, tempProfiles);
}

void ScoreRuler::extend(double inpScore, int seed, double eps) {
    random_engine_t mtGen(static_cast<uint32_t>(seed));

    // fill currentSample
    for (unsigned elemIndex = 0; elemIndex < sampleSize; elemIndex++) {
        std::vector<int> comb = combination(0, n - 1, genesetSize, mtGen);
        currentSample[elemIndex] = std::vector<unsigned>(comb.begin(), comb.end());
        currentProfiles[elemIndex] = getProfile(expressionMatrix, currentSample[elemIndex], m);
    }

    duplicateSampleElements();

    while (scores.back() <= inpScore - 1e-10){
        int moves = 0;
        int totalIters = 0;
        for (moves = 0; moves < sampleSize * genesetSize;) {
            for (unsigned elemIndex = 0; elemIndex < sampleSize; elemIndex++) {
                moves += updateElement(currentSample[elemIndex], currentProfiles[elemIndex],
                                       scores.back(), mtGen);
                totalIters += itersPerStep;
            }

            if (double(moves)/totalIters < 0.01) {
                break;
            }
        }

        if (double(moves)/totalIters < 0.01) {
            break;
        }

        duplicateSampleElements();
        if (eps != 0){
            unsigned long k = scores.size() / ((sampleSize + 1) / 2);
            if (k > - log2(0.5 * eps)) {
                break;
            }
        }
    }
}

std::pair<double, double> ScoreRuler::getPvalue(double inpScore, double eps){
    unsigned long halfSize = (sampleSize + 1) / 2;

    auto it = scores.begin();
    if (inpScore >= scores.back()){
        it = scores.end() - 1;
    }
    else{
        it = lower_bound(scores.begin(), scores.end(), inpScore);
    }

    unsigned long indx = 0;
    (it - scores.begin()) > 0 ? (indx = (it - scores.begin())) : indx = 0;

    unsigned long k = (indx) / halfSize;
    unsigned long remainder = sampleSize -  (indx % halfSize);

    double adjLog = betaMeanLog(halfSize, sampleSize);
    double adjLogPval = k * adjLog + betaMeanLog(remainder + 1, sampleSize);

    double pval = std::max(0.0, std::min(1.0, exp(adjLogPval)));
    double log2err = multilevelError(k+1, sampleSize);

    if (inpScore > scores.back()) {
        log2err = std::numeric_limits<double>::infinity();
    }

    return std::make_pair(pval, log2err);
}


int ScoreRuler::updateElement(std::vector<unsigned> & element,
                              std::vector<float> & profile,
                              double threshold,
                              random_engine_t &mtGen){
    // unsigned n = expressionMatrix.size();

    uid_wrapper uid_n(0, n - 1, mtGen);
    uid_wrapper uid_k(0, genesetSize - 1, mtGen);

    std::vector<bool> used(n);
    for (auto el: element) {
        used[el] = 1;
    }

    unsigned niters = itersPerStep;
    int moves = 0;
    std::vector<float> newProfile(profile.size());
    for (unsigned i = 0; i < niters; i++){
        unsigned toDrop = uid_k();
        unsigned indxOld = element[toDrop];
        unsigned indxNew = uid_n();

        if (used[indxNew]) {
            continue;
        }

        adjustProfile(expressionMatrix, profile, newProfile, indxNew, indxOld, m);
        double newScore = getScore(newProfile);

        if (newScore >= threshold){
            used[element[toDrop]] = 0;
            used[indxNew] = 1;
            element[toDrop] = indxNew;
            std::swap(profile, newProfile);
            moves++;
        }

    }
    return moves;
}
