#pragma once
#include <vector>
#include <random>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>

// generate k numbers from [a, b] - closed interval
// a should be non-negative, usually 0 or 1
std::vector<int> combination(const int &a, const int &b, const int &k, std::mt19937& rng);

struct uid_wrapper {
	#ifdef USE_STD_UID
	std::mt19937& rng;
	std::uniform_int_distribution<int> uid;
	#else
	int from, len;
	std::mt19937& rng;
	unsigned completePart;
	#endif
	uid_wrapper(int, int, std::mt19937&);

	int operator()();
};

double betaMeanLog(unsigned long a, unsigned long b);
