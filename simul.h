#ifndef SIMUL_H_
#define SIMUL_H_

#include <cassert>
#include <random>

class simul {
public:

 simul(unsigned int seed) : engine(seed), dist(0.0, 1.0) {
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

	  if (l >= len) {
	  	assert(arr[len - 1] == 0.0); 
	  	l = int(random() * len);
	  }

	  return l;
	}

	double random() { return dist(engine); }

private:
	std::mt19937 engine;
	std::uniform_real_distribution<double> dist;
};

#endif /* SIMUL_H_ */

