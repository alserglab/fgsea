#pragma once

#include <vector>
#include <random>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/random/mersenne_twister.hpp>

using random_engine_t = boost::mt19937;

// generate k numbers from [a, b] - closed interval
// a should be non-negative, usually 0 or 1
std::vector<int> combination(const int &a, const int &b, const int &k, random_engine_t& rng);

struct uid_wrapper {
	#ifdef USE_STD_UID
	random_engine_t& rng;
	std::uniform_int_distribution<int> uid;
	#else
	int from, len;
	random_engine_t& rng;
	unsigned completePart;
	#endif
	uid_wrapper(int, int, random_engine_t&);

	int operator()();
};

double betaMeanLog(unsigned long a, unsigned long b);

double multilevelError(int level, int sampleSize);
