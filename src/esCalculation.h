#ifndef FGSEAMULTILEVELCPP_ESCALCULATION_H
#define FGSEAMULTILEVELCPP_ESCALCULATION_H

#include <vector>
#include <cmath>

using namespace std;

struct score_t {
  //  score = coef_ns / NS - coef_const / diff
  int64_t NS;
  int64_t coef_NS;
  int64_t diff;   //  n - k
  int64_t coef_const;
  
  double getDouble() const {
      return 1.0 * coef_NS / NS - 1.0 * coef_const / diff;
  }

  int64_t getNumerator() const {
    return coef_NS * diff - coef_const * NS;
  }

  __int128_t Compare(score_t const& other) const {
    int64_t P1 = coef_NS * diff + NS * (other.coef_const - coef_const);
    int64_t Q1 = NS * diff;
    int64_t P2 = other.coef_NS;
    int64_t Q2 = other.NS;
    return __int128_t(P1) * Q2 - __int128_t(P2) * Q1;
  }

  bool operator<(score_t const& other) const {
    return Compare(other) < 0;
  }

  bool operator<=(score_t const& other) const {
    return Compare(other) <= 0;
  }

  bool operator>(score_t const& other) const {
    return Compare(other) > 0;
  }

  bool operator>=(score_t const& other) const {
    return Compare(other) >= 0;
  }

  bool operator==(score_t const& other) const {
    return Compare(other) == 0;
  }

  bool operator!=(score_t const& other) const {
    return Compare(other) != 0;
  }

  static int64_t getMaxNS() {
      return int64_t(1e9);
  }

  score_t operator-() const {
    return score_t{NS, -coef_NS, diff, -coef_const};
  }

  score_t abs() const {
    return std::max(*this, -(*this));
  }
};

score_t calcES(const vector<int64_t> &ranks, const vector<int> &p, int64_t NS);

score_t calcES(const vector<int64_t> &ranks, const vector<int> &p);

score_t calcPositiveES(const vector<int64_t> &ranks, const vector<int> &p, int64_t NS);

score_t calcPositiveES(const vector<int64_t> &ranks, const vector<int> &p);

bool compareStat(const vector<double> &ranks, const vector<int> &p, double NS, double bound);

#endif //FGSEAMULTILEVELCPP_ESCALCULATION_H
