#include "fgseaMultilevelSupplement.h"
#include <Rcpp.h>

void my_assert(bool x, string msg) {
    if (!x) {
        std::cerr << "FAIL" << " " << msg << std::endl;
        throw;
    }
}

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

EsRuler::EsRuler(const vector<double> &inpRanks,
                 unsigned int inpSampleSize,
                 unsigned int inpPathwaySize,
                 double inpMovesScale,
                 bool inpLog) :
    logStatus(inpLog),
    doubleRanks(inpRanks),
    ranks(scaleRanks(doubleRanks)),
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
    double const scale = MAX_NS / curNS;
    if (logStatus) {
        Rcpp::Rcout << "Scaling ranks by " << scale << endl;
    }
    for (int i = 0; i < int(ranks.size()); ++i) {
        result[i] = int64_t(ranks[i] * scale);
    }
    return result;
}

bool EsRuler::resampleGenesets(mt19937 &rng) {
    vector<pair<gsea_t, int>> stats(sampleSize);
    vector<int> posEsIndxs;
    int totalPosEsCount = 0;
    
    for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
        auto sampleEsPos = calcPositiveES(ranks, currentSamples[sampleId]);
        auto sampleEs = calcES(ranks, currentSamples[sampleId]);
        if (sampleEs.getNumerator() > 0) {
            totalPosEsCount++;
            posEsIndxs.push_back(sampleId);
        }
        hash_t sampleHash = calcHash(currentSamples[sampleId]);
        stats[sampleId] = {gsea_t{sampleEsPos, sampleHash}, sampleId};
    }
    sort(stats.begin(), stats.end());
    
    int startFrom = 0;

    auto choose_median = [&] {
        auto centralValue = stats[sampleSize / 2].first;
        for (int sampleId = 0; sampleId < sampleSize; sampleId++){
            if (stats[sampleId].first >= centralValue) {
                startFrom = sampleId;
                break;
            }
        }

        if (startFrom == 0) {
            while (startFrom < sampleSize && stats[startFrom].first == stats[0].first) {
                ++startFrom;
            }
        }
    };
    auto choose_min = [&] {
        while (startFrom < sampleSize && stats[startFrom].first == stats[0].first) {
            ++startFrom;
        }
    };

    choose_median();
    // choose_min();
    
    if (startFrom == sampleSize) {
        if (logStatus) {
            Rcpp::Rcout << "all values are equal!\n";
        }
        return true;
    }
    
    levels.emplace_back();
    levels.back().bound = stats[startFrom - 1].first;   //  greater
    for (int i = 0; i < startFrom; ++i) {
        levels.back().lowScores.push_back(stats[i].first);
    }
    for (int i = startFrom; i < sampleSize; ++i) {
        levels.back().highScores.push_back(stats[i].first);
    }

    uniform_int_distribution<> uid(0, sampleSize - startFrom - 1);

    auto gen_new_sample = [&] {
        int ind = uid(rng) + startFrom;
        return currentSamples[stats[ind].second];
    };

    vector<vector<int> > new_sets;
    for (int i = 0; i < startFrom; i++){
        new_sets.push_back(gen_new_sample());
    }
    for (int i = startFrom; i < sampleSize; ++i) {
        new_sets.push_back(currentSamples[stats[i].second]);
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

void EsRuler::extend(double ES, int seed, double eps) {
    mt19937 gen(seed);
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

    int evMovesCnt = movesScale * pathwaySize;
    poisson_distribution<> distr(evMovesCnt);

    //  maybe here -1 is needed
    score_t NEED_ES{score_t::getMaxNS(), int64_t(roundl(score_t::getMaxNS() * ES)) - 1, n - k, 0};
    
    for (int levelNum = 1; levels.back().bound.first < NEED_ES; ++levelNum) {
        if (logStatus) {
            Rcpp::Rcout << std::setprecision(15) << std::fixed << "Iteration " << levelNum << ": " << levels.back().bound.first.getDouble() << ' ' << levels.back().bound.second << ' ' << ES << endl;
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
                my_assert(cnt < chunksNumber, "chunks constructed incorrectly");
                samplesChunks[i].chunks[cnt].push_back(pos);
                samplesChunks[i].chunkSum[cnt] += ranks[pos];
            }
        }

        auto run_full = [&] {
            int iterations = 0;
            int needMoves = movesScale * sampleSize * pathwaySize / 2; 
            int moves = 0;
            int iters = 0;
            for (; moves < needMoves; iterations++) {
                for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
                    auto nSuccess = perturbate(ranks, k, samplesChunks[sampleId], levels.back().bound, gen);
                    moves += nSuccess.moves;
                    iters += nSuccess.iters;
                }
            }
            for (int i = 0; i < iterations; i++) {
                for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
                    auto nSuccess = perturbate(ranks, k, samplesChunks[sampleId], levels.back().bound, gen);
                    moves += nSuccess.moves;
                    iters += nSuccess.iters;
                }
            }
            if (logStatus) {
                Rcpp::Rcout << "Acceptance rate: " << 1.0 * moves / iters << std::endl;
            }
        };
        auto run_partial = [&] {
            for (int sampleId = 0; sampleId < oldSamplesStart; sampleId++) {
                int movesCnt = distr(gen);
                auto iters = perturbate_success(ranks, k, samplesChunks[sampleId], levels.back().bound, gen, movesCnt / 2).iters;
                perturbate_iters(ranks, k, samplesChunks[sampleId], levels.back().bound, gen, iters);
            }
        };
        run_full();
        // run_partial();
        
        for (int i = 0; i < sampleSize; ++i) {
            currentSamples[i].clear();
            for (int j = 0; j < chunksNumber; ++j) {
                for (int pos : samplesChunks[i].chunks[j]) {
                    currentSamples[i].push_back(pos);
                }
            }
            my_assert(std::is_sorted(currentSamples[i].begin(), currentSamples[i].end()), "chunks constructed incorrectly");
            if (logStatus) {
                auto gsea_score = calcPositiveES(ranks, currentSamples[i]);
                auto hash_value = calcHash(currentSamples[i]);
                gsea_t score{gsea_score, hash_value};
                if (score <= levels.back().bound) {
                    Rcpp::Rcout << std::setprecision(15) << std::fixed << score.first.getDouble() << " " << score.second << std::endl;
                }
                my_assert(score > levels.back().bound, "Perturbed sets have scores < bound");
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

        /*if (eps != 0){
            unsigned long k = enrichmentScores.size() / ((sampleSize + 1) / 2);
            if (k > - log2(0.5 * eps)) {
                break;
            }
        }*/
    }
    if (logStatus) {
        Rcpp::Rcout << "Done: " << levels.back().bound.first.getDouble() << " " << NEED_ES.getDouble() << endl;
    }
}

tuple <double, bool, double> EsRuler::getPvalue(double ES_double, double eps, bool sign){
    if (incorrectRuler){
        return make_tuple(std::nan("1"), true, std::nan("1"));
    }

    int const n = (int) ranks.size();
    int const k = pathwaySize;
    score_t ES_score{score_t::getMaxNS(), int64_t(roundl(score_t::getMaxNS() * ES_double)), n - k, 0};
    gsea_t ES{ES_score, 0};
    int lastLevel = 0;
    while (lastLevel < int(levels.size()) - 1 && levels[lastLevel].bound < ES) {
        ++lastLevel;
    }
    //  [0..lastLevel) => P(highScore)
    //  lastLevel => P(x >= ES)

    double adjLogPval = 0;
    double lvlsVar = 0;

    for (int i = 0; i < lastLevel; ++i) {
        auto& lvl = levels[i];
        int nhigh = int(lvl.highScores.size());
        nhigh += 1;
        adjLogPval += betaMeanLog(nhigh, sampleSize);
        lvlsVar += getVarPerLevel(nhigh, sampleSize);
    }
    int cntLast = 0;
    auto& lastLvl = levels[lastLevel];
    for (auto x : lastLvl.lowScores) {
        cntLast += (x > ES);
    }
    for (auto x : lastLvl.highScores) {
        cntLast += (x > ES);
    }
    adjLogPval += betaMeanLog(cntLast, sampleSize);
    lvlsVar += getVarPerLevel(cntLast, sampleSize);
    
    double log2err = sqrt(lvlsVar) / log(2);
    return make_tuple(max(0.0, min(1.0, exp(adjLogPval))), true, log2err);
}

int EsRuler::chunkLen(int ind) {
    if (ind == 0) {
        return chunkLastElement[0];
    }
    return chunkLastElement[ind] - chunkLastElement[ind - 1];
}

EsRuler::PerturbateResult EsRuler::perturbate(vector<int64_t> const& ranks, int k, EsRuler::SampleChunks& sampleChunks, gsea_t bound, mt19937 &rng) {
    double pertPrmtr = 0.1;
    int iters = max(1, (int) (k * pertPrmtr));
    return perturbate_iters(ranks, k, sampleChunks, bound, rng, iters);
}

EsRuler::PerturbateResult EsRuler::perturbate_iters(vector<int64_t> const& ranks, int k, EsRuler::SampleChunks& sampleChunks, gsea_t bound, mt19937 &rng, int need_iters) {
    return perturbate_until(ranks, k, sampleChunks, bound, rng, [need_iters](int moves, int iters) {
        return iters >= need_iters;
    });
}

EsRuler::PerturbateResult EsRuler::perturbate_success(vector<int64_t> const& ranks, int k, EsRuler::SampleChunks& sampleChunks, gsea_t bound, mt19937 &rng, int need_successes) {
    return perturbate_until(ranks, k, sampleChunks, bound, rng, [need_successes](int moves, int iters) {
        return moves >= need_successes;
    });
}

EsRuler::PerturbateResult EsRuler::perturbate_until(vector<int64_t> const& ranks, int k, EsRuler::SampleChunks& sampleChunks, gsea_t bound, mt19937 &rng, std::function<bool(int, int)> const& f) {
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