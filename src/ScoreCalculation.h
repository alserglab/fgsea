#ifndef SCORE_CALC_H
#define SCORE_CALC_H
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>

double getVarianceEstimator(const std::vector<float> & elems);

double getScore(const std::vector<float> & profile);

std::vector<float> getProfile(const std::vector<float> & E,
                               const std::vector<unsigned> & indexes,
                               unsigned m);

void adjustProfile(const std::vector<float> & E,
                  const std::vector<float> & profile,
                  std::vector<float> & newProfile,
                  unsigned idNew, unsigned idOld,
                  unsigned m);

#endif //SCORE_CALC_H
