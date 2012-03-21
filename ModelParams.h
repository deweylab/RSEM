#ifndef MODELPARAMS_H_
#define MODELPARAMS_H_

#include<cstdio>
#include<cstring>

#include "utils.h"
#include "Refs.h"

struct ModelParams {
	int M;
	READ_INT_TYPE N[3];
	int minL, maxL;
	bool estRSPD; // true if user wants to estimate RSPD; false if use uniform distribution
	int B; // number of bins in RSPD
	int mate_minL, mate_maxL;
	double probF; //probability of forward strand
	double mean, sd;
	Refs *refs;

	int seedLen;

	//default parameters
	ModelParams() {
		minL = 1; maxL = 1000;
		estRSPD = false;
		B = 20; // default bin size if estRSPD is true
		mate_minL = 1; mate_maxL = 1000;
		probF = 0.5;
		mean = -1; sd = 0;

		M = 0;
		memset(N, 0, sizeof(N));
		refs = NULL;

		seedLen = 0;
	}
};
#endif /* MODELPARAMS_H_ */
