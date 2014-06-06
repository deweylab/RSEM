#ifndef SAMPLING
#define SAMPLING

#include<ctime>
#include<cstdio>
#include<cassert>
#include<vector>
#include<set>

#include "boost/random.hpp"

typedef unsigned int seedType;
typedef boost::random::mt19937 engine_type;
typedef boost::random::uniform_01<> uniform_01_dist;
typedef boost::random::gamma_distribution<> gamma_dist;
typedef boost::random::variate_generator<engine_type&, uniform_01_dist> uniform_01_generator;
typedef boost::random::variate_generator<engine_type&, gamma_dist> gamma_generator;

class engineFactory {
public:
  static void init() { seedEngine = new engine_type(time(NULL)); }
  static void init(seedType seed) { seedEngine = new engine_type(seed); }

  static void finish() { if (seedEngine != NULL) delete seedEngine; }

	static engine_type *new_engine() {
		seedType seed;
		static std::set<seedType> seedSet;			// empty set of seeds
		std::set<seedType>::iterator iter;

		do {
			seed = (*seedEngine)();
			iter = seedSet.find(seed);
		} while (iter != seedSet.end());
		seedSet.insert(seed);

		return new engine_type(seed);
	}

 private:
	static engine_type *seedEngine;
};

engine_type* engineFactory::seedEngine = NULL;

// arr should be cumulative!
// interval : [,)
// random number should be in [0, arr[len - 1])
// If by chance arr[len - 1] == 0.0, one possibility is to sample uniformly from 0...len-1
int sample(uniform_01_generator& rg, std::vector<double>& arr, int len) {
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
