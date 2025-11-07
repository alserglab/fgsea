#include "fgseaMultilevelSupplement.h"
#include <Rcpp.h>

double getVarPerLevel(unsigned long k, unsigned long n){
    return boost::math::trigamma(k) - boost::math::trigamma(n + 1);
}

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

EsRuler::EsRuler(const vector<int64_t> &inpRanks,
                 unsigned int inpSampleSize,
                 unsigned int inpPathwaySize,
                 double inpMovesScale,
                 bool inpLog) :
    logStatus(inpLog),
    // doubleRanks(inpRanks),
    // ranks(scaleRanks(doubleRanks)),
    ranks(inpRanks),
    geneHashes(inpRanks.size()),
    sampleSize(inpSampleSize),
    pathwaySize(inpPathwaySize),
    movesScale(inpMovesScale) {
    currentSamples.resize(inpSampleSize);

}

EsRuler::~EsRuler() = default;

vector<int64_t> EsRuler::scaleRanks(vector<double> const& ranks) {
    int64_t const MAX_NS = score_t::getMaxNS();
    double const curNS = accumulate(ranks.begin(), ranks.end(), double(0));
    vector<int64_t> result(ranks.size());
    // keep input integers integer
    double const scale = floorl(MAX_NS / curNS);
    if (logStatus) {
        Rcpp::Rcout << "Scaling ranks by " << scale << endl;
    }
    for (int i = 0; i < int(ranks.size()); ++i) {
        result[i] = int64_t(ranks[i] * scale);
    }
    return result;
}

bool EsRuler::resampleGenesets(random_engine_t &rng) {
    //  <score, is_positive, ind>
    vector<tuple<gsea_t, int, int>> stats(sampleSize);

    for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
        auto sampleEsPos = calcPositiveES(ranks, currentSamples[sampleId]);
        auto sampleEs = calcES(ranks, currentSamples[sampleId]);
        hash_t sampleHash = calcHash(currentSamples[sampleId]);
        stats[sampleId] = {gsea_t{sampleEsPos, sampleHash},
                           (sampleEs.getNumerator() >= 0),
                           sampleId};
    }
    sort(stats.begin(), stats.end());

    int startFrom = 0;
    auto centralValue = get<0>(stats[sampleSize / 2]);
    for (int sampleId = 0; sampleId < sampleSize; sampleId++){
        if (get<0>(stats[sampleId]) >= centralValue) {
            startFrom = sampleId;
            break;
        }
    }

    if (startFrom == 0) {
        while (startFrom < sampleSize && get<0>(stats[startFrom]) == get<0>(stats[0])) {
            ++startFrom;
        }
    }

    if (startFrom == sampleSize) {
        if (logStatus) {
            Rcpp::Rcout << "Got all equal values. Ending multilevel process\n";
        }
        return true;
    }

    levels.emplace_back();
    levels.back().bound = get<0>(stats[startFrom - 1]);   //  greater
    for (int i = 0; i < startFrom; ++i) {
        levels.back().lowScores.emplace_back(get<0>(stats[i]), get<1>(stats[i]));
    }
    for (int i = startFrom; i < sampleSize; ++i) {
        levels.back().highScores.emplace_back(get<0>(stats[i]), get<1>(stats[i]));
    }

    uid_wrapper uid(0, sampleSize - startFrom - 1, rng);

    auto gen_new_sample = [&] {
        int ind = uid() + startFrom;
        return currentSamples[get<2>(stats[ind])];
    };

    vector<vector<int> > new_sets;
    for (int i = 0; i < startFrom; i++){
        new_sets.push_back(gen_new_sample());
    }
    for (int i = startFrom; i < sampleSize; ++i) {
        new_sets.push_back(currentSamples[get<2>(stats[i])]);
    }

    oldSamplesStart = startFrom;
    swap(currentSamples, new_sets);
    return true;
}

EsRuler::SampleChunks::SampleChunks(int chunksNumber) : chunkSum(chunksNumber), chunks(chunksNumber) {}

EsRuler::hash_t EsRuler::calcHash(const vector<int>& curSample) {
    hash_t res = 0;
    for (int i : curSample) {
        res ^= geneHashes[i];
    }
    return res;
}

void EsRuler::extend(double ES_double, int seed, double eps) {
    random_engine_t gen(static_cast<uint32_t>(seed));
    int const n = (int) ranks.size();
    int const k = pathwaySize;

    for (int i = 0; i < n; ++i) {
        geneHashes[i] = gen();
    }

    for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
        currentSamples[sampleId] = combination(0, ranks.size() - 1, pathwaySize, gen);
        sort(currentSamples[sampleId].begin(), currentSamples[sampleId].end());
    }

    if (!resampleGenesets(gen)) {
        if (logStatus) {
            Rcpp::Rcout << "Could not advance in the start" << endl;
        }
        incorrectRuler = true;
        return;
    }

    chunksNumber = max(1, (int) sqrt(pathwaySize));
    chunkLastElement = vector<int>(chunksNumber);
    chunkLastElement[chunksNumber - 1] = ranks.size();
    vector<int> tmp(sampleSize);
    vector<SampleChunks> samplesChunks(sampleSize, SampleChunks(chunksNumber));

    // flooring ES an exact score, but 1.0 should stay 1.0
    // works best if getMaxNS is a power of two
    score_t NEED_ES{score_t::getMaxNS(), int64_t(score_t::getMaxNS() * ES_double), 1, 0};

    double adjLogPval = 0;
    for (int levelNum = 1; levels.back().bound.first < NEED_ES; ++levelNum) {
        adjLogPval += betaMeanLog(int(levels.back().highScores.size() + 1), sampleSize);
        if (eps != 0 && adjLogPval < log(eps)) {
            break;
        }

        if (logStatus) {
            Rcpp::Rcout << std::setprecision(15) << std::fixed << "Iteration " << levelNum << ": score=" << levels.back().bound.first.getDouble() << ", hash=" << levels.back().bound.second << endl;
        }

        for (int i = 0, pos = 0; i < chunksNumber - 1; ++i) {
            pos += (pathwaySize + i) / chunksNumber;
            for (int j = 0; j < sampleSize; ++j) {
                tmp[j] = currentSamples[j][pos];
            }
            nth_element(tmp.begin(), tmp.begin() + sampleSize / 2, tmp.end());
            chunkLastElement[i] = tmp[sampleSize / 2];
        }

        for (int i = 0; i < sampleSize; ++i) {
            fill(samplesChunks[i].chunkSum.begin(), samplesChunks[i].chunkSum.end(), int64_t(0));
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

        int nIterations = 0;
        int nAccepted = 0;
        int needAccepted = movesScale * sampleSize * pathwaySize / 2;
        for (; nAccepted < needAccepted; nIterations++) {
            for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
                auto perturbResult = perturbate(ranks, k, samplesChunks[sampleId], levels.back().bound, gen);
                nAccepted += perturbResult.moves;
            }
        }
        for (int i = 0; i < nIterations; i++) {
            for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
                perturbate(ranks, k, samplesChunks[sampleId], levels.back().bound, gen);
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

        auto lastSize = levels.size();
        if (!resampleGenesets(gen)) {
            incorrectRuler = true;
            if (logStatus) {
                Rcpp::Rcout << "Could not advance after level " << levelNum << endl;
            }
        }
        if (lastSize == levels.size()) {
            break;
        }
    }
}

// p-value, true, error
// TODO: remove bool
tuple <double, bool, double> EsRuler::getPvalue(double ES_double, double eps, bool sign){
    if (incorrectRuler){
        return make_tuple(std::nan("1"), true, std::nan("1"));
    }

    int const n = (int) ranks.size();
    int const k = pathwaySize;
    // flooring ES an exact score, but 1.0 should stay 1.0
    // works best if getMaxNS is a power of two
    score_t ES_score{score_t::getMaxNS(), int64_t(score_t::getMaxNS() * ES_double), 1, 0};
    gsea_t ES{ES_score, 0};
    //  Calculate Pr[<score, hash> >= ES]

    double adjLogPval = 0;
    double lvlsVar = 0;

    for (auto& lvl : levels) {
        if (ES <= lvl.bound) {
            int cntLast = 0;
            int cntPositive = 0;
            //  highScores > lvl.bound >= ES
            for (auto[x, isPositive] : lvl.highScores) {
                cntLast += 1;
                cntPositive += isPositive;
            }
            for (auto[x, isPositive] : lvl.lowScores) {
                if (x >= ES) {
                    cntLast += 1;
                    cntPositive += isPositive;
                }
            }

            int numerator = (sign ? cntLast : cntPositive);
            if (numerator == 0) {
                adjLogPval += betaMeanLog(1, sampleSize);
                return make_tuple(max(0.0, min(1.0, exp(adjLogPval))), true, std::nan("1"));
            }

            adjLogPval += betaMeanLog(numerator, sampleSize);
            lvlsVar += getVarPerLevel(numerator, sampleSize);

            double log2err = sqrt(lvlsVar) / log(2);
            return make_tuple(max(0.0, min(1.0, exp(adjLogPval))), true, log2err);
        }

        int nhigh = int(lvl.highScores.size());
        nhigh += 1;
        adjLogPval += betaMeanLog(nhigh, sampleSize);
        lvlsVar += getVarPerLevel(nhigh, sampleSize);
    }

    auto& lastLevel = levels.back();
    //  ES > lastLevel.bound
    //  highScores -> uniform [> lastLevel.bound]
    int cntLast = 0;
    int cntPositive = 0;
    for (auto[x, isPositive] : lastLevel.highScores) {
        if (x >= ES) {
            cntLast += 1;
            cntPositive += isPositive;
        }
    }

    int numerator = (sign ? cntLast : cntPositive);

    if (numerator == 0) {
        adjLogPval += betaMeanLog(1, int(lastLevel.highScores.size()));
        return make_tuple(max(0.0, min(1.0, exp(adjLogPval))), true, std::nan("1"));
    }

    adjLogPval += betaMeanLog(numerator, int(lastLevel.highScores.size()));
    lvlsVar += getVarPerLevel(numerator, int(lastLevel.highScores.size()));

    double log2err = sqrt(lvlsVar) / log(2);
    return make_tuple(max(0.0, min(1.0, exp(adjLogPval))), true, log2err);
}

int EsRuler::chunkLen(int ind) {
    if (ind == 0) {
        return chunkLastElement[0];
    }
    return chunkLastElement[ind] - chunkLastElement[ind - 1];
}

EsRuler::PerturbateResult EsRuler::perturbate(vector<int64_t> const& ranks, int k, EsRuler::SampleChunks& sampleChunks, gsea_t bound, random_engine_t &rng) {
    double pertPrmtr = 0.1;
    int iters = max(1, (int) (k * pertPrmtr));
    return perturbate_iters(ranks, k, sampleChunks, bound, rng, iters);
}

EsRuler::PerturbateResult EsRuler::perturbate_iters(vector<int64_t> const& ranks, int k, EsRuler::SampleChunks& sampleChunks, gsea_t bound, random_engine_t &rng, int need_iters) {
    return perturbate_until(ranks, k, sampleChunks, bound, rng, [need_iters](int moves, int iters) {
        return iters >= need_iters;
    });
}

EsRuler::PerturbateResult EsRuler::perturbate_success(vector<int64_t> const& ranks, int k, EsRuler::SampleChunks& sampleChunks, gsea_t bound, random_engine_t &rng, int need_successes) {
    return perturbate_until(ranks, k, sampleChunks, bound, rng, [need_successes](int moves, int iters) {
        return moves >= need_successes;
    });
}

EsRuler::PerturbateResult EsRuler::perturbate_until(vector<int64_t> const& ranks, int k, EsRuler::SampleChunks& sampleChunks, gsea_t bound, random_engine_t &rng, std::function<bool(int, int)> const& f) {
    int n = int(ranks.size());
    uid_wrapper uid_n(0, n - 1, rng);
    uid_wrapper uid_k(0, k - 1, rng);

    int64_t NS = 0;
    hash_t curHash = 0;
    for (int i = 0; i < (int) sampleChunks.chunks.size(); ++i) {
        for (int pos : sampleChunks.chunks[i]) {
            NS += ranks[pos];
            curHash ^= geneHashes[pos];
        }
    }
    int candVal = -1;
    bool hasCand = false;
    int candX = 0;
    int64_t candY = 0;

    int moves = 0;
    int iters = 0;
    while (!f(moves, iters)) {
        iters += 1;
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

        // TODO: might be beneficial to replace with linear search for small sizes
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
        curHash ^= geneHashes[oldVal] ^ geneHashes[newVal];
        sampleChunks.chunkSum[oldChunkInd] -= ranks[oldVal];
        sampleChunks.chunkSum[newChunkInd] += ranks[newVal];

        bool strictly = (curHash <= bound.second);

        auto check = [&](score_t const& score) {
            return strictly ? score > bound.first : score >= bound.first;
        };

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

        if (hasCand && check(score_t{NS, candY, n - k, candX})) {
            ++moves;
            continue;
        }

        int curX = 0;
        int64_t curY = 0;
        bool ok = false;
        int last = -1;

        bool fl = false;
        for (int i = 0; i < (int) sampleChunks.chunks.size(); ++i) {
            if (!check(score_t{NS, curY + sampleChunks.chunkSum[i], n - k, curX})) {
                curY += sampleChunks.chunkSum[i];
                curX += chunkLastElement[i] - last - 1 - (int) sampleChunks.chunks[i].size();
                last = chunkLastElement[i] - 1;
            } else {
                for (int pos : sampleChunks.chunks[i]) {
                    curY += ranks[pos];
                    curX += pos - last - 1;
                    if (check(score_t{NS, curY, n - k, curX})) {
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
            curHash ^= geneHashes[newVal] ^ geneHashes[oldVal];

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
    return {moves, iters};
}
