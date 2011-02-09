#ifndef HITCONTAINER_H_
#define HITCONTAINER_H_

#include<cassert>
#include<iostream>
#include<vector>

#include<algorithm>
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

	int getN() { return n; }

	int getNHits() { return nhits; }

	int calcNumGeneMultiReads(const GroupInfo&);
	int calcNumIsoformMultiReads();

	int getSAt(int pos) { assert(pos >= 0 && pos <= n); return s[pos]; }

	HitType& getHitAt(int pos) { assert(pos >= 0 && pos < nhits); return hits[pos]; }

private:
	int n; // n reads in total
	int nhits; // # of hits
	std::vector<int> s;
	std::vector<HitType> hits;
};

//Each time only read one read's hits. If you want to start over, must call clear() first!
template<class HitType>
bool HitContainer<HitType>::read(std::istream& in) {
	int tot;

	if (!(in>>tot)) return false;
	assert(tot > 0);
	for (int i = 0; i < tot; i++) {
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
	for (int i = 0; i < n; i++) {
		out<<s[i + 1] - s[i];
		for (int j = s[i]; j < s[i + 1]; j++) {
			hits[j].write(out);
		}
		out<<std::endl;
	}
}

template<class HitType>
int HitContainer<HitType>::calcNumGeneMultiReads(const GroupInfo& gi) {
	int res = 0;
	int *sortgids = NULL;

	for (int i = 0; i < n; i++) {
		int num = s[i + 1] - s[i];
		sortgids = new int[num];
		for (int j = s[i]; j < s[i + 1]; j++) sortgids[j] = gi.gidAt(hits[j].getSid());
		std::sort(sortgids, sortgids + num);
		if (std::unique(sortgids, sortgids + num) - sortgids > 1) ++res;
		delete[] sortgids;
	}

	return res;
}

template<class HitType>
int HitContainer<HitType>::calcNumIsoformMultiReads() {
	int res = 0;
	for (int i = 0; i < n; i++)
		if (s[i + 1] - s[i] > 1) ++res;
	return res;
}

#endif /* HITCONTAINER_H_ */
