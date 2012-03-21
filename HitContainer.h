#ifndef HITCONTAINER_H_
#define HITCONTAINER_H_

#include<cassert>
#include<iostream>
#include<vector>
#include<algorithm>

#include "utils.h"
#include "GroupInfo.h"

template<class HitType>
class HitContainer {
public:
	HitContainer() {
		clear();
	}

	void clear() {
		n = nhits = 0;
		s.clear();
		hits.clear();

		s.push_back(0);
	}

	bool read(std::istream&); // each time a read
	void write(std::ostream&); // write all reads' hit out

	void push_back(const HitType& hit)  {
		hits.push_back(hit);
		++nhits;
	}

	//update read information vector etc
	void updateRI() {
		if (nhits > s.back()) {  //Do not change if last read does not have hits
			s.push_back(nhits);
			++n;
		}
	}

	READ_INT_TYPE getN() { return n; }

	HIT_INT_TYPE getNHits() { return nhits; }

	READ_INT_TYPE calcNumGeneMultiReads(const GroupInfo&);
	READ_INT_TYPE calcNumIsoformMultiReads();

	HIT_INT_TYPE getSAt(READ_INT_TYPE pos) { assert(pos >= 0 && pos <= n); return s[pos]; }

	HitType& getHitAt(HIT_INT_TYPE pos) { assert(pos >= 0 && pos < nhits); return hits[pos]; }

private:
	READ_INT_TYPE n; // n reads in total
	HIT_INT_TYPE nhits; // # of hits
	std::vector<HIT_INT_TYPE> s;
	std::vector<HitType> hits;
};

//Each time only read one read's hits. If you want to start over, must call clear() first!
template<class HitType>
bool HitContainer<HitType>::read(std::istream& in) {
	HIT_INT_TYPE tot;

	if (!(in>>tot)) return false;
	assert(tot > 0);
	for (HIT_INT_TYPE i = 0; i < tot; i++) {
		HitType hit;
		if (!hit.read(in)) return false;
		hits.push_back(hit);
	}

	nhits = nhits + tot;
	++n;
	s.push_back(nhits);

	return true;
}

template<class HitType>
void HitContainer<HitType>::write(std::ostream& out) {
	if (n <= 0) return;
	for (READ_INT_TYPE i = 0; i < n; i++) {
		out<<s[i + 1] - s[i];
		for (HIT_INT_TYPE j = s[i]; j < s[i + 1]; j++) {
			hits[j].write(out);
		}
		out<<std::endl;
	}
}

template<class HitType>
READ_INT_TYPE HitContainer<HitType>::calcNumGeneMultiReads(const GroupInfo& gi) {
	READ_INT_TYPE res = 0;
	int *sortgids = NULL;

	for (READ_INT_TYPE i = 0; i < n; i++) {
		HIT_INT_TYPE num = s[i + 1] - s[i];
		sortgids = new int[num];
		for (HIT_INT_TYPE j = s[i]; j < s[i + 1]; j++) sortgids[j] = gi.gidAt(hits[j].getSid());
		std::sort(sortgids, sortgids + num);
		if (std::unique(sortgids, sortgids + num) - sortgids > 1) ++res;
		delete[] sortgids;
	}

	return res;
}

template<class HitType>
READ_INT_TYPE HitContainer<HitType>::calcNumIsoformMultiReads() {
	READ_INT_TYPE res = 0;
	for (READ_INT_TYPE i = 0; i < n; i++)
		if (s[i + 1] - s[i] > 1) ++res;
	return res;
}

#endif /* HITCONTAINER_H_ */
