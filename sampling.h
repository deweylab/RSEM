#ifndef SAMPLING
#define SAMPLING

#include<ctime>
#include<cstdio>
#include<cassert>
#include<vector>
#include<set>

#include "boost/random.hpp"

typedef unsigned int seedType;
typedef boost::mt19937 engine_type;
typedef boost::gamma_distribution<> gamma_dist;
typedef boost::uniform_01<engine_type> uniform01;
typedef boost::variate_generator<engine_type&, gamma_dist> gamma_generator;

class engineFactory {
public:
	static engine_type *new_engine() {
		seedType seed;
		static engine_type seedEngine(time(NULL));
		static std::set<seedType> seedSet;			// empty set of seeds
		std::set<seedType>::iterator iter;

		do {
			seed = seedEngine();
			iter = seedSet.find(seed);
		} while (iter != seedSet.end());
		seedSet.insert(seed);

		return new engine_type(seed);
	}
};

// arr should be cumulative!
// interval : [,)
// random number should be in [0, arr[len - 1])
// If by chance arr[len - 1] == 0.0, one possibility is to sample uniformly from 0...len-1
int sample(uniform01& rg, std::vector<double>& arr, int len) {
  int l, r, mid;
  double prb = rg() * arr[len - 1];

  l = 0; r = len - 1;
  while (l <= r) {
    mid = (l + r) / 2;
    if (arr[mid] <= prb) l = mid + 1;
    else r = mid - 1;
  }

  if (l >= len) { printf("%d %lf %lf\n", len, arr[len - 1], prb); }
  assert(l < len);

  return l;
}

#endif
