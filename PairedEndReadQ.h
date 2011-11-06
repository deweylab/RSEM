#ifndef PAIREDENDREADQ_H_
#define PAIREDENDREADQ_H_

#include<cassert>
#include<iostream>
#include<string>

#include "Read.h"
#include "SingleReadQ.h"

class PairedEndReadQ : public Read {
public:
	PairedEndReadQ() : mate1(), mate2() {}
	PairedEndReadQ(const SingleReadQ& mate1, const SingleReadQ& mate2) {
		this->mate1 = mate1;
		this->mate2 = mate2;
		this->name = mate1.getName();

		calc_lq();
	}

	bool read(int argc, std::istream* argv[], int flags = 7);
	void write(int argc, std::ostream* argv[]);

	const SingleReadQ& getMate1() const { return mate1; }
	const SingleReadQ& getMate2() const { return mate2; }
	const SingleReadQ& getMate(int i) const {
		if (i == 1) return mate1;
		else return mate2;
	}

private:
	SingleReadQ mate1, mate2;

	void calc_lq();
};

bool PairedEndReadQ::read(int argc, std::istream* argv[], int flags) {
	bool success;
    std::istream *inpMate1[1], *inpMate2[1];

	assert(argc == 2);
	inpMate1[0] = argv[0]; inpMate2[0] = argv[1];
	success = mate1.read(1, inpMate1, flags) && mate2.read(1, inpMate2, flags);
	name = "";
	if (flags & 4) { name = mate1.getName(); } //May chop 1 char later if we want

	if (flags & 1) calc_lq();

	return success;
}

void PairedEndReadQ::write(int argc, std::ostream* argv[]) {
	std::ostream *outMate1[1], *outMate2[1];

	assert(argc == 2);
	outMate1[0] = argv[0]; outMate2[0] = argv[1];
	mate1.write(1, outMate1);
	mate2.write(1, outMate2);
}

void PairedEndReadQ::calc_lq() {
	low_quality = mate1.isLowQuality() && mate2.isLowQuality();
	if (mate1.getReadLength() < OLEN || mate2.getReadLength() < OLEN) low_quality = true;
}

#endif /* PAIREDENDREADQ_H_ */
