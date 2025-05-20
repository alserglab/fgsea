#include "esCalculation.h"
#include <algorithm>

score_t calcES(const vector<int64_t> &ranks, const vector<int> &p, int64_t NS) {
    // p must be sorted
    int n = (int) ranks.size();
    int k = (int) p.size();
    score_t res{NS, 0, n - k, 0};
    score_t cur{NS, 0, n - k, 0};
    int last = -1;
    for (int pos : p) {
        cur.coef_const += pos - last - 1;
        if (res.abs() < cur.abs()) {
            res = cur;
        }
        cur.coef_NS += ranks[pos];
        if (res.abs() < cur.abs()) {
            res = cur;
        }
        last = pos;
    }
    return res;
}

score_t calcES(const vector<int64_t> &ranks, const vector<int> &p) {
    // p must be sorted
    int64_t NS = 0;
    for (int pos : p) {
        NS += ranks[pos];
    }
    return calcES(ranks, p, NS);
}

score_t calcPositiveES(const vector<int64_t> &ranks, const vector<int> &p, int64_t NS) {
    // p must be sorted
    int n = (int) ranks.size();
    int k = (int) p.size();
    score_t res{NS, 0, n - k, 0};
    score_t cur{NS, 0, n - k, 0};
    int last = -1;
    for (int pos : p) {
        cur.coef_NS += ranks[pos];
        cur.coef_const += pos - last - 1;
        res = max(res, cur);
        last = pos;
    }
    return res;
}

score_t calcPositiveES(const vector<int64_t> &ranks, const vector<int> &p) {
    // p must be sorted
    int64_t NS = 0;
    for (int pos : p) {
        NS += ranks[pos];
    }
    return calcPositiveES(ranks, p, NS);
}

bool compareStat(const vector<double> &ranks, const vector<int> &p, double NS, double bound){
    // p must be sorted
    int n = (int) ranks.size();
    int k = (int) p.size();
    double cur = 0.0;
    double q1 = 1.0 / (n - k);
    double q2 = 1.0 / NS;
    int last = -1;
    for (int pos : p) {
        cur += q2 * ranks[pos] - q1 * (pos - last - 1);
        if (cur > bound) {
            return true;
        }
        last = pos;
    }
    return false;
}
