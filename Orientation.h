#ifndef ORIENTATION_H_
#define ORIENTATION_H_

#include<cstdio>
#include<cstring>
#include<cassert>

#include "simul.h"

class Orientation {
public:
	Orientation(double probF = 0.5) {
		prob[0] = probF;
		prob[1] = 1.0 - probF;
	}

	Orientation& operator= (const Orientation& rv) {
		if (this == &rv) return *this;
		memcpy(prob, rv.prob, sizeof(rv.prob));
		return *this;
	}

	//dir : 0 + 1 -
	double getProb(int dir) { return prob[dir]; }

	void read(FILE* fi) {
		assert(fscanf(fi, "%lf", &prob[0]) == 1);
		prob[1] = 1.0 - prob[0];
	}

	void write(FILE* fo) {
		fprintf(fo, "%.10g\n", prob[0]);
	}


	int simulate(simul* sampler) { return (sampler->random() < prob[0] ? 0 : 1); }

private:
	double prob[2]; //0 + 1 -
};

#endif /* ORIENTATION_H_ */
