#include "fgseaMultilevelSupplement.h"
#include "esCalculation.h"
#include "util.h"

pair<double, bool> calcLogCorrection(const vector<unsigned int> &probCorrector,
                                     long probCorrIndx, unsigned int sampleSize){
    double result = 0.0;

    unsigned long halfSize = (sampleSize + 1) / 2;
    unsigned long remainder = sampleSize - probCorrIndx % (halfSize);

    double condProb = betaMeanLog(probCorrector[probCorrIndx] + 1, remainder);
    result += condProb;

    if (exp(condProb) >= 0.5){
        return make_pair(result, true);
    }
    else{
        return make_pair(result, false);
    }
}

EsRuler::EsRuler(const vector<double> &inpRanks, unsigned int inpSampleSize, unsigned int inpPathwaySize) :
    ranks(inpRanks), sampleSize(inpSampleSize), pathwaySize(inpPathwaySize) {
    currentSamples.resize(inpSampleSize);
}

EsRuler::~EsRuler() = default;

void EsRuler::duplicateSamples() {
    /*
     * Removes samples with an enrichment score less than the median value and
     * replaces them with samples with an enrichment score greater than the median
     * value
    */
    vector<pair<double, int> > stats(sampleSize);
    vector<int> posEsIndxs;
    int totalPosEsCount = 0;

    for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
        double sampleEsPos = calcPositiveES(ranks, currentSamples[sampleId]);
        double sampleEs = calcES(ranks, currentSamples[sampleId]);
        if (sampleEs > 0) {
            totalPosEsCount++;
            posEsIndxs.push_back(sampleId);
        }
        stats[sampleId] = make_pair(sampleEsPos, sampleId);
    }
    sort(stats.begin(), stats.end());
    for (int sampleId = 0; 2 * sampleId < sampleSize; sampleId++) {
        enrichmentScores.push_back(stats[sampleId].first);
        if (find(posEsIndxs.begin(), posEsIndxs.end(), stats[sampleId].second) != posEsIndxs.end()) {
            totalPosEsCount--;
        }
        probCorrector.push_back(totalPosEsCount);
    }

    vector<vector<int> > new_sets;
    for (int sampleId = 0; 2 * sampleId < sampleSize - 2; sampleId++) {
        for (int rep = 0; rep < 2; rep++) {
            new_sets.push_back(currentSamples[stats[sampleSize - 1 - sampleId].second]);
        }
    }
    new_sets.push_back(currentSamples[stats[sampleSize >> 1].second]);
    swap(currentSamples, new_sets);
}

EsRuler::SampleChunks::SampleChunks(int chunksNumber) : chunkSum(chunksNumber), chunks(chunksNumber) {}

void EsRuler::extend(double ES, int seed, double eps) {
    random_engine_t gen(seed);

    for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
        currentSamples[sampleId] = combination(0, ranks.size() - 1, pathwaySize, gen);
        sort(currentSamples[sampleId].begin(), currentSamples[sampleId].end());

        double currentES = calcES(ranks, currentSamples[sampleId]);
    }

    chunksNumber = max(1, (int) sqrt(pathwaySize));
    chunkLastElement = vector<int>(chunksNumber);
    chunkLastElement[chunksNumber - 1] = ranks.size();
    vector<int> tmp(sampleSize);
    vector<SampleChunks> samplesChunks(sampleSize, SampleChunks(chunksNumber));

    duplicateSamples();
    while (enrichmentScores.back() <= ES - 1e-10){
        for (int i = 0, pos = 0; i < chunksNumber - 1; ++i) {
            pos += (pathwaySize + i) / chunksNumber;
            for (int j = 0; j < sampleSize; ++j) {
                tmp[j] = currentSamples[j][pos];
            }
            nth_element(tmp.begin(), tmp.begin() + sampleSize / 2, tmp.end());
            chunkLastElement[i] = tmp[sampleSize / 2];
        }

        for (int i = 0; i < sampleSize; ++i) {
            fill(samplesChunks[i].chunkSum.begin(), samplesChunks[i].chunkSum.end(), 0.0);
            for (int j = 0; j < chunksNumber; ++j) {
                samplesChunks[i].chunks[j].clear();
            }
            int cnt = 0;
            for (int pos : currentSamples[i]) {
                while (chunkLastElement[cnt] <= pos) {
                    ++cnt;
                }
                samplesChunks[i].chunks[cnt].push_back(pos);
                samplesChunks[i].chunkSum[cnt] += ranks[pos];
            }
        }

        for (int moves = 0; moves < sampleSize * pathwaySize;) {
            for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
                moves += perturbate(ranks, pathwaySize, samplesChunks[sampleId], enrichmentScores.back(), gen);
            }
        }

        for (int i = 0; i < sampleSize; ++i) {
            currentSamples[i].clear();
            for (int j = 0; j < chunksNumber; ++j) {
                for (int pos : samplesChunks[i].chunks[j]) {
                    currentSamples[i].push_back(pos);
                }
            }
        }

        double prevTopScore = enrichmentScores.back();
        duplicateSamples();
        if (enrichmentScores.back() <= prevTopScore) {
            // couldn't advance the median
            // :TODO: should be more accurately addressed
            break;
        }
        if (eps != 0){
            unsigned long k = enrichmentScores.size() / ((sampleSize + 1) / 2);
            if (k > - log2(0.5 * eps)) {
                break;
            }
        }
    }
}


pair<double, bool> EsRuler::getPvalue(double ES, double eps, bool sign) {
    unsigned long halfSize = (sampleSize + 1) / 2;

    auto it = enrichmentScores.begin();
    bool goodError = true;
    if (ES >= enrichmentScores.back()){
        it = enrichmentScores.end() - 1;
        if (ES > enrichmentScores.back() + 1e-10) {
            goodError = false; // went beyond the ruler
        }
    }
    else{
        it = lower_bound(enrichmentScores.begin(), enrichmentScores.end(), ES);
    }

    unsigned long indx = 0;
    (it - enrichmentScores.begin()) > 0 ? (indx = (it - enrichmentScores.begin())) : indx = 0;

    unsigned long k = (indx) / halfSize;
    unsigned long remainder = sampleSize -  (indx % halfSize);

    double adjLog = betaMeanLog(halfSize, sampleSize);
    double adjLogPval = k * adjLog + betaMeanLog(remainder + 1, sampleSize);

    if (sign) {
        return make_pair(max(0.0, min(1.0, exp(adjLogPval))), goodError);
    } else {
        pair<double, bool> correction = calcLogCorrection(probCorrector, indx, sampleSize);
        double resLog = adjLogPval + correction.first;
        return make_pair(max(0.0, min(1.0, exp(resLog))), goodError && correction.second);
    }
}

int EsRuler::chunkLen(int ind) {
    if (ind == 0) {
        return chunkLastElement[0];
    }
    return chunkLastElement[ind] - chunkLastElement[ind - 1];
}

int EsRuler::perturbate(const vector<double> &ranks, int k, EsRuler::SampleChunks &sampleChunks,
               double bound, random_engine_t &rng) {
    double pertPrmtr = 0.1;
    int n = (int) ranks.size();
    uid_wrapper uid_n(0, n - 1, rng);
    uid_wrapper uid_k(0, k - 1, rng);
    double NS = 0;
    for (int i = 0; i < (int) sampleChunks.chunks.size(); ++i) {
        for (int pos : sampleChunks.chunks[i]) {
            NS += ranks[pos];
        }
    }
    double q1 = 1.0 / (n - k);
    int iters = max(1, (int) (k * pertPrmtr));
    int moves = 0;

    int candVal = -1;
    bool hasCand = false;
    int candX = 0;
    double candY = 0;

    for (int i = 0; i < iters; i++) {
        int oldInd = uid_k();
        int oldChunkInd = 0, oldIndInChunk = 0;
        int oldVal;
        {
            int tmp = oldInd;
            while ((int) sampleChunks.chunks[oldChunkInd].size() <= tmp) {
                tmp -= sampleChunks.chunks[oldChunkInd].size();
                ++oldChunkInd;
            }
            oldIndInChunk = tmp;
            oldVal = sampleChunks.chunks[oldChunkInd][oldIndInChunk];
        }

        int newVal = uid_n();

        int newChunkInd = upper_bound(chunkLastElement.begin(), chunkLastElement.end(), newVal) - chunkLastElement.begin();
        int newIndInChunk = lower_bound(sampleChunks.chunks[newChunkInd].begin(), sampleChunks.chunks[newChunkInd].end(), newVal) - sampleChunks.chunks[newChunkInd].begin();

        if (newIndInChunk < (int) sampleChunks.chunks[newChunkInd].size() && sampleChunks.chunks[newChunkInd][newIndInChunk] == newVal) {
            if (newVal == oldVal) {
                ++moves;
            }
            continue;
        }

        sampleChunks.chunks[oldChunkInd].erase(sampleChunks.chunks[oldChunkInd].begin() + oldIndInChunk);
        sampleChunks.chunks[newChunkInd].insert(
            sampleChunks.chunks[newChunkInd].begin() + newIndInChunk - (oldChunkInd == newChunkInd && oldIndInChunk < newIndInChunk ? 1 : 0),
            newVal);

        NS = NS - ranks[oldVal] + ranks[newVal];

        sampleChunks.chunkSum[oldChunkInd] -= ranks[oldVal];

        sampleChunks.chunkSum[newChunkInd] += ranks[newVal];

        if (hasCand) {
            if (oldVal == candVal) {
                hasCand = false;
            }
        }
        if (hasCand) {
            if (oldVal < candVal) {
                candX++;
                candY -= ranks[oldVal];
            }
            if (newVal < candVal) {
                candX--;
                candY += ranks[newVal];
            }
        }

        double q2 = 1.0 / NS;

        if (hasCand && -q1 * candX + q2 * candY > bound) {
            ++moves;
            continue;
        }

        int curX = 0;
        double curY = 0;
        bool ok = false;
        int last = -1;

        bool fl = false;
        for (int i = 0; i < (int) sampleChunks.chunks.size(); ++i) {
            if (q2 * (curY + sampleChunks.chunkSum[i]) - q1 * curX < bound) {
                curY += sampleChunks.chunkSum[i];
                curX += chunkLastElement[i] - last - 1 - (int) sampleChunks.chunks[i].size();
                last = chunkLastElement[i] - 1;
            } else {
                for (int pos : sampleChunks.chunks[i]) {
                    curY += ranks[pos];
                    curX += pos - last - 1;
                    if (q2 * curY - q1 * curX > bound) {
                        ok = true;
                        hasCand = true;
                        candX = curX;
                        candY = curY;
                        candVal = pos;
                        break;
                    }
                    last = pos;
                }
                if (ok) {
                    break;
                }
                curX += chunkLastElement[i] - 1 - last;
                last = chunkLastElement[i] - 1;
            }
        }

        if (!ok) {
        	NS = NS - ranks[newVal] + ranks[oldVal];

            sampleChunks.chunkSum[oldChunkInd] += ranks[oldVal];

            sampleChunks.chunkSum[newChunkInd] -= ranks[newVal];

            sampleChunks.chunks[newChunkInd].erase(
                sampleChunks.chunks[newChunkInd].begin() + newIndInChunk - (oldChunkInd == newChunkInd && oldIndInChunk < newIndInChunk ? 1 : 0));
            sampleChunks.chunks[oldChunkInd].insert(sampleChunks.chunks[oldChunkInd].begin() + oldIndInChunk, oldVal);

            if (hasCand) {
                if (newVal == candVal) {
                    hasCand = false;
                }
            }
            if (hasCand) {
                if (oldVal < candVal) {
                    candX--;
                    candY += ranks[oldVal];
                }
                if (newVal < candVal) {
                    candX++;
                    candY -= ranks[newVal];
                }
            }
        } else {
            ++moves;
        }
    }
    return moves;
}
