#include "util.h"
#include <vector>
#include <random>
#include <algorithm>

// generate k numbers from [a, b] - closed interval
// a should be non-negative, usually 0 or 1
std::vector<int> combination(const int &a, const int &b, const int &k, std::mt19937& rng) {
    // std::uniform_int_distribution<int> uni(a, b);
    uid_wrapper uni(a, b, rng);
    std::vector<int> v;
    v.reserve(k);

    int n = b - a + 1;
    std::vector<char> used(n);

    if (k < n * 1.0 / 2){
        for (int i = 0; i < k; i++) {
            for (int j = 0; ; j++) { // average < 2
                int x = uni();
                if (!used[x - a]) {
                    v.push_back(x);
                    used[x - a] = true;
                    break;
                }
            }
        }
    } else {
        for (int r = n - k; r < n; ++r){
            // int x = std::uniform_int_distribution<>(0, r)(rng);
            int x = uid_wrapper(0, r, rng)();
            if (!used[x]){
                v.push_back(a + x);
                used[x] = true;
            } else{
                v.push_back(a + r);
                used[r] = true;
            }
        }
        std::shuffle(v.begin(), v.end(), rng);
    }

    return v;
}

#ifdef USE_STD_UID
uid_wrapper::uid_wrapper(int _from, int _to, std::mt19937& _rng) : rng(_rng), uid(_from, _to) {}

int uid_wrapper::operator()() {
    return uid(rng);
}
#else
uid_wrapper::uid_wrapper(int _from, int _to, std::mt19937& _rng) : from(_from), len(_to - _from + 1), rng(_rng) {
    unsigned maxVal = rng.max();
    completePart = maxVal - maxVal % len;
}

int uid_wrapper::operator()() {
    unsigned x;
    do {
        x = rng();
    } while (x >= completePart);

    return from + x % len;
}
#endif

double betaMeanLog(unsigned long a, unsigned long b) {
    return boost::math::digamma(a) - boost::math::digamma(b + 1);
}
