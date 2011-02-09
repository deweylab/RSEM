#ifndef SIMUL_MERSENNE_H_
#define SIMUL_MERSENNE_H_

#include<ctime>

#include "randomc.h"
#include "simul.h"

class simul_mersenne : public simul {
public:
	simul_mersenne() : rg(time(NULL)) { }

	virtual ~simul_mersenne() {}
	virtual double random() { return rg.Random(); }

private:
	CRandomMersenne rg;
};

#endif /* SIMUL_MERSENNE_H_ */
