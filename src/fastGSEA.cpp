#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <random>
#include <chrono>

#include "util.h"

using namespace std;

const double eps = 1e-13;

// not actually a segment tree, but a square root heuristic
template <class T>
class SegmentTree {
    private:
        vector<T> t;
        vector<T> b;
        int n;
        int k;
        int k2;
        int logK;
        int blockMask;
    public:
        SegmentTree(int n_) {
            k = 1;
            logK = 0;
            while (k * k < n_) {
                k <<= 1;
                logK++;
            }
            k2 = (n_ - 1) / k + 1;
            n = k * k;

            blockMask = k - 1;

            t = vector<T>(n, 0);
            b = vector<T>(k2, 0);
        }

        // works in O(sqrt(n))
        void inc(int p, T delta) {  // increase value at position p
            int blockEnd = p - (p & blockMask) + blockMask + 1;
            for (; p < blockEnd; ++p) {
                t[p] += delta;
            }

            for (int p1 = (p >> logK); p1 < k2; ++p1) {
                b[p1] += delta;
            }
        }

        // works in O(1)
        T queryR(int r) {  // sum on interval [0, r)
            if (r == 0) {
                return 0;
            } else {
                r--;
            }
            return t[r] + b[r >> logK];
        }

};

template <class T>
class IndirectCmp {
    public:
        IndirectCmp(T const & x) : x(x) { }

    bool operator() (int i, int j) {
        return x[i] < x[j];
    }
    private:
        T const & x;
};

template<class T>
vector<int> order(T const & x) {
    vector<int> res(x.size());
    for (int i = 0; i < (int)x.size(); ++i) {
        res[i] = i;
    }

    IndirectCmp<T> cmp(x);
    sort(res.begin(), res.end(), cmp);
    return res;
}

vector<int> ranksFromOrder(vector<int> const& order) {
    vector<int> res(order.size());

    for (int i = 0; i < (int)order.size(); ++i) {
        res[order[i]] = i;
    }

    return res;
}


NumericVector gseaStats1(
        NumericVector const& stats,
        IntegerVector const& selectedStats,
        vector<int> const& selectedOrder,
        double gseaParam,
        bool rev = false) {
    int n = stats.size();
    double NR = 0;
    int k = selectedStats.size();

    NumericVector res(k);

    SegmentTree<int> xs(k + 1);
    SegmentTree<double> ys(k + 1);

    vector<int> selectedRanks = ranksFromOrder(selectedOrder);

    if (!rev) {
        int prev = -1;
        for (int i = 0; i < k; ++i) {
            int j = selectedOrder[i];
            int t = selectedStats[j] - 1;
            assert(t - prev >= 1);
            xs.inc(i, t - prev);
            prev = t;
        }
        xs.inc(k, n - 1 - prev);
    } else {
        int prev = n;
        for (int i = 0; i < k; ++i) {
            selectedRanks[i] = k - 1 - selectedRanks[i];
            int j = selectedOrder[k - 1 - i];
            int t = selectedStats[j] - 1;
            assert(prev - t>= 1);
            xs.inc(i, prev - t);
            prev = t;
        }
        xs.inc(k, prev - 0);

    }

    // previous node in stack
    vector<int> stPrev(k + 2);
    vector<int> stNext(k + 2);

    int k1 = (int)(sqrt(k + 1));
    int k2 = (k + 1) / k1 + 1;
    vector<int> blockSummit(k2);
    vector<int> blockStart(k2);
    vector<int> blockEnd(k2);


    for (int i = 0; i <= k + 1; i += k1) {
        int block = i / k1;
        blockStart[block] = i;
        blockEnd[block] = min(i + k1 - 1, k + 1);
        for (int j = 1; i + j <= blockEnd[block]; ++j) {
            stPrev[i + j] = i + j - 1;
            stNext[i + j - 1] = i + j;
        }
        stPrev[i] = i;
        stNext[blockEnd[block]] = blockEnd[block];
        blockSummit[block] = blockEnd[block];
    }

    double statEps = 1e-5;

    for (int i = 0; i < k; ++i) {
        int t = selectedStats[i] - 1;
        double xx = abs(stats[t]);
        if (xx > 0) {
            statEps = min(xx, statEps);
        }
    }
    statEps /= 1024;

    for (int i = 0; i < k; ++i) {
        int t = selectedStats[i] - 1;
        int tRank = selectedRanks[i];
        // cout << tRank << ":\n";
        // 0 values make problems, replacing with epsilon
        double adjStat = pow(max(abs(stats[t]), statEps), gseaParam);

        xs.inc(tRank, -1);
        ys.inc(tRank, adjStat);
        NR += adjStat;

        int m = i + 1;


        int curBlock = (tRank + 1) / k1;
        int bS = blockStart[curBlock];
        int bE = blockEnd[curBlock];

        // Redoing upper convext hull for invalidated block `curBlock`

        int curTop = max(tRank, bS);

        for (int j = tRank + 1; j <= bE; ++j) {
            int c = j;
            double xc = xs.queryR(c);
            double yc = ys.queryR(c);

            int b = curTop;

            double xb = xs.queryR(b);
            double yb = ys.queryR(b);

            while (stPrev[curTop] != curTop) {
                int a = stPrev[curTop];

                double xa = xs.queryR(a);
                double ya = ys.queryR(a);


                double pr = (xb - xa) * (yc - yb) - (yb - ya) * (xc - xb);
                if (yc - ya < eps) {
                    pr = 0;
                }
                if (pr <= 0) {
                    // right turn
                    break;
                }
                // left turn
                curTop = a;
                stNext[b] = -1;
                b = a;
                xb = xa;
                yb = ya;
            }

            stPrev[c] = curTop;
            stNext[curTop] = c;
            curTop = c;
            if (stNext[c] != -1) {
                break;
            }
        }

        for (int j = bS; j < bE; ++j) {
            // cout << "assert: " << j << " " << stNext[j] << " " << stPrev[stNext[j]] << "\n";
            // assert(stNext[j] == -1 || stPrev[stNext[j]] == j);
        }

        double coef = (double)(n - m) / NR;

        // Finding farthest points for `curBlock` from scratch

        // int topSummit = 0;
        double maxP = 0;

        // So that blockSummit[curBlock] is to the right of
        // actual summit as for every other block
        blockSummit[curBlock] = max(curTop, blockSummit[curBlock]);

        // Updating farthest points in valid blocks

        for (int block = 0; block < k2; ++block) {

            int curSummit = blockSummit[block];

            double curDist = ys.queryR(curSummit) * coef - xs.queryR(curSummit);

            while (1) {
                int nextSummit = stPrev[curSummit];
                double nextDist =
                    ys.queryR(nextSummit) * coef -
                    xs.queryR(nextSummit);

                if (nextDist <= curDist) {
                    break;
                }
                curDist = nextDist;
                curSummit = nextSummit;
            }

            blockSummit[block] = curSummit;
            /*
            if (curDist > maxP) {
                topSummit = curSummit;
            }
            */
            maxP = max(maxP, curDist);
        }


        maxP /= (double)(n - m);

        res[i] = maxP;

        // Checking correctness

        /*
        double maxP1 = 0;
        int topSummit1 = 0;

        for (int j = 1; j <= k + 1; ++j) {
            int c = j;
            double x = xs.queryR(c) / (double)(n - m);
            double y = ys.queryR(c) / NR;


            if (y - x > maxP1) {
                topSummit1 = j;
            }
            maxP1 = max(maxP1, y - x);
        }
        cout << "diff: " << maxP1  -  maxP << " " << maxP1 << " " << maxP << "\n";
        cout << "topSummit: " << topSummit << "\n";
        cout << "topSummit1: " << topSummit1 << "\n";
        assert(topSummit == topSummit1);
        */
    }
    return res;
}

NumericVector subvector(NumericVector const &from, IntegerVector const &indices) {
    NumericVector result(indices.size());
    for (int i = 0; i < indices.size(); ++i) {
        result[i] = from[indices[i] - 1];
    }
    return result;
}

NumericVector calcGseaStatCumulative(
        NumericVector const& stats,
        IntegerVector const& selectedStats, // Indexes start from one!
        double gseaParam,
        std::string scoreType
) {

    vector<int> selectedOrder = order(selectedStats);

    if (scoreType != "std" && scoreType != "pos" && scoreType != "neg"){
        throw std::invalid_argument("scoreType must take values from (\"std\", \"pos\", \"neg\")");
    }

    if (scoreType == "std"){
        NumericVector res = gseaStats1(stats, selectedStats, selectedOrder, gseaParam);
        NumericVector resDown = gseaStats1(stats, selectedStats, selectedOrder, gseaParam, true);
        for (int i = 0; i < (int)selectedStats.size(); ++i) {
            if (res[i] == resDown[i]) {
                res[i] = 0;
            } else if (res[i] < resDown[i]) {
                res[i] = -resDown[i];
            }
        }
        return res;
    }
    else if(scoreType == "pos"){
        NumericVector res = gseaStats1(stats, selectedStats, selectedOrder, gseaParam);
        return res;
    }
    else{
        NumericVector res = -gseaStats1(stats, selectedStats, selectedOrder, gseaParam, true);
        return res;
    }

}

NumericVector calcRandomGseaStatCumulative(
        NumericVector const& stats,
        int n,
        int k,
        double gseaParam,
        random_engine_t& rng,
        std::string scoreType
) {

    IntegerVector selectedStats = wrap(combination(1, n, k, rng));

    return calcGseaStatCumulative(stats, selectedStats, gseaParam, scoreType);
}

List calcGseaStatCumulativeBatch(
        NumericVector const& stats,
        double gseaParam,
        NumericVector const& pathwayScores,
        IntegerVector const& pathwaysSizes,
        int iterations,
        int seed,
        std::string scoreType) {

    int n = stats.size();
    int k = max(pathwaysSizes);
    int m = pathwaysSizes.size();

    NumericVector leEs(m);
    NumericVector geEs(m);
    NumericVector leZero(m);
    NumericVector geZero(m);
    NumericVector leZeroSum(m);
    NumericVector geZeroSum(m);

    NumericVector zeros(m);
    LogicalVector aux;
    NumericVector diff;

    random_engine_t rng(static_cast<uint32_t>(seed));

    for (int i = 0; i < iterations; ++i) {
        NumericVector randEs = calcRandomGseaStatCumulative(stats, n, k, gseaParam, rng, scoreType);
        NumericVector randEsP = subvector(randEs, pathwaysSizes);

        aux = randEsP <= pathwayScores;
        diff = wrap(aux);
        leEs = leEs + diff;

        aux = randEsP >= pathwayScores;
        diff = wrap(aux);
        geEs = geEs + diff;

        aux = randEsP <= zeros;
        diff = wrap(aux);
        leZero = leZero + diff;

        aux = randEsP >= zeros;
        diff = wrap(aux);
        geZero = geZero + diff;

        leZeroSum = leZeroSum + pmin(randEsP, zeros);
        geZeroSum = geZeroSum + pmax(randEsP, zeros);
    }
    return Rcpp::List::create(
            Rcpp::Named("leEs") = leEs,
            Rcpp::Named("geEs") = geEs,
            Rcpp::Named("leZero") = leZero,
            Rcpp::Named("geZero") = geZero,
            Rcpp::Named("leZeroSum") = leZeroSum,
            Rcpp::Named("geZeroSum") = geZeroSum
    );
}

NumericVector calcGseaStatBatchCpp(
        NumericVector const & stats,
        List const & selectedGenes,
        IntegerVector const & geneRanks) {
    int n = stats.size();
    int m = selectedGenes.size();
    NumericVector gseaStats(m);

    for (int i = 0; i < m; ++i)  {
        vector<int> S(as<vector<int> >(selectedGenes[i]));
        for (int j = 0; j < S.size(); ++j) {
            S[j] = geneRanks[S[j]-1];
        }
        sort(S.begin(), S.end());

        int k = S.size();
        gseaStats[i] = k;


        double NR = 0;
        for (int j = 0; j < k; ++j) {
            int tRank = S[j]; // 1-based indexing
            double cur = abs(stats[tRank - 1]);
            NR += cur;
        }

        double x = 0;
        double y = 0;

        double coef = (double)(n - k) / NR;

        double maxP = 0;
        double minP = 0;


        for (int j = 0; j < k; ++j) {
            int tRank = S[j]; // 1-based indexing
            double cur = abs(stats[tRank - 1]);

            x = tRank - j - 1;

            double bottom = y * coef - x;

            y += cur;

            double top = y * coef - x;
            maxP = max(maxP, top);
            minP = min(minP, bottom);
        }

        if(maxP > -minP) {
            gseaStats[i] = maxP;
        } else if (maxP < -minP) {
            gseaStats[i] = minP;
        } else {
            gseaStats[i] = 0;
        }
        gseaStats[i] /= (n - k);
    }

    return gseaStats;
}


