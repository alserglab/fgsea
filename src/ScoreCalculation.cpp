#include "ScoreCalculation.h"

double getVarianceEstimator(const std::vector<float> & elems){
    double varHat = 0.0;
    for (double elem : elems){
        varHat += elem * elem;
    }
    return varHat;
}


std::vector<float> getProfile(const std::vector<float> & E,
                               const std::vector<unsigned> & indexes,
                               unsigned m){
    // unsigned m = E[0].size(); // number or rows must be nonzero
    std::vector<float> profile(m);

    for (unsigned j = 0; j < m; j++){
        double value = 0.0;
        for (unsigned indx : indexes){
            value += E[indx*m + j];
        }
        profile[j] = value;
    }
    return profile;
}


double getScore(const std::vector<float> & profile){
    return(getVarianceEstimator(profile));
}



void adjustProfile(const std::vector<float> & E,
                                  const std::vector<float> & profile,
                                  std::vector<float> & newProfile,
                                  unsigned idNew, unsigned idOld,
                                  unsigned m){
    for (unsigned i = 0; i < newProfile.size(); i++){
        newProfile[i] = profile[i] - E[idOld*m + i] + E[idNew*m + i];
    }
}
