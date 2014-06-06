#ifndef SIMUL_H_
#define SIMUL_H_

#include<cassert>

#include "boost/random.hpp"

class simul {
public:

 simul(unsigned int seed) : engine(seed), rg(engine, boost::random::uniform_01<>()) {
  }

	// interval : [,)
	// random number should be in [0, arr[len - 1])
	// If by chance arr[len - 1] == 0.0, one possibility is to sample uniformly from 0 ... len - 1
	int sample(double* arr, int len) {
	  int l, r, mid;
	  double prb = random() * arr[len - 1];


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

	double random() { return rg(); };

private:
	boost::random::mt19937 engine;
	boost::random::variate_generator<boost::random::mt19937&, boost::random::uniform_01<> > rg;
};

#endif /* SIMUL_H_ */

