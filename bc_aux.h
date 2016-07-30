#ifndef BC_AUX_H_
#define BC_AUX_H_

#include<map>

#include <stdint.h>
#include "htslib/sam.h"

struct SingleEndT {
	bam1_t *b;

	SingleEndT(bam1_t *b) {
		this->b = b;
	}

	int getSign(bool value) const { return value ? -1 : 1; }

	int compare(const SingleEndT& o) const {
		int strand1, strand2;
		uint32_t *p1, *p2;

		if (b->core.tid != o.b->core.tid) return getSign(b->core.tid < o.b->core.tid);
		if (b->core.pos != o.b->core.pos) return getSign(b->core.pos < o.b->core.pos);
		strand1 = bam_is_rev(b); strand2 = bam_is_rev(o.b);
		if (strand1 != strand2) return getSign(strand1 < strand2);
		if (b->core.n_cigar != o.b->core.n_cigar) return getSign(b->core.n_cigar < o.b->core.n_cigar);
		p1 = bam_get_cigar(b); p2 = bam_get_cigar(o.b);
		for (int i = 0; i < (int)b->core.n_cigar; ++i) {
			if (*p1 != *p2) return getSign(*p1 < *p2);
			++p1; ++p2;
		}

		return 0;
	}

	bool operator< (const SingleEndT& o) const {
		return compare(o) < 0;
	}
};

struct PairedEndT {
	SingleEndT mate1, mate2;

	PairedEndT(const SingleEndT& mate1, const SingleEndT& mate2) : mate1(mate1), mate2(mate2) {
	}

	bool operator< (const PairedEndT& o) const {
		int value = mate1.compare(o.mate1);
		return value < 0 || (value == 0 && mate2 < o.mate2);
	}
};

class CollapseMap {
public:
	CollapseMap() { isPaired = false; smap.clear(); pmap.clear(); }

	void init(bool isPaired) {
		this->isPaired = isPaired;
		isPaired ? pmap.clear() : smap.clear();
	}

	void insert(bam1_t *b, bam1_t *b2, float prb) {
		if (!isPaired) {
			smapIter = smap.find(SingleEndT(b));
			if (smapIter == smap.end()) { smap[SingleEndT(bam_dup1(b))] = prb; }
			else smapIter->second += prb;
		}
		else {
			pmapIter = pmap.find(PairedEndT(SingleEndT(b), SingleEndT(b2)));
			if (pmapIter == pmap.end()) { pmap[PairedEndT(SingleEndT(bam_dup1(b)), SingleEndT(bam_dup1(b2)))] = prb; }
			else pmapIter->second += prb;
		}
	}

	//once this function is called, "insert" cannot be called anymore
	bool empty(bool& par) {
		bool value;

		par = isPaired;
		if (!isPaired) { value = smap.empty(); smapIter = smap.begin(); }
		else { value = pmap.empty(); pmapIter = pmap.begin(); }

		return value;
	}

	bool next(bam1_t*& b, bam1_t*& b2, float& prb) {
		bool value;

		if (!isPaired) {
			value = smapIter != smap.end();
			if (value) {
				b = smapIter->first.b;
				prb = smapIter->second;
				smapIter++;
			}
		}
		else {
			value = pmapIter != pmap.end();
			if (value) {
				b = pmapIter->first.mate1.b;
				b2 = pmapIter->first.mate2.b;
				prb = pmapIter->second;
				pmapIter++;
			}
		}

		return value;
	}

private:
	bool isPaired;

	std::map<SingleEndT, float> smap;
	std::map<SingleEndT, float>::iterator smapIter;

	std::map<PairedEndT, float> pmap;
	std::map<PairedEndT, float>::iterator pmapIter;
};

#endif /* BC_AUX_H_ */
