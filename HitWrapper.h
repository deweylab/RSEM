#ifndef HITWRAPPER_H_
#define HITWRAPPER_H_

#include "utils.h"
#include "HitContainer.h"

// assume each hit vector contains at least one hit

template<class HitType>
class HitWrapper {
public:
	HitWrapper(int nThreads, HitContainer<HitType> **hitvs) {
		this->nThreads = nThreads;
		this->hitvs = hitvs;
		i = 0; j = 0;
	}

	HitType* getNextHit() {
		HitType *res;

		if (i >= nThreads) return NULL;
		res = &(hitvs[i]->getHitAt(j));
		++j;
		if (j >= hitvs[i]->getNHits()) { ++i; j = 0; }

		return res;
	}

private:
	int i, nThreads;
	HIT_INT_TYPE j;
	HitContainer<HitType> **hitvs;
};

#endif /* HITWRAPPER_H_ */
