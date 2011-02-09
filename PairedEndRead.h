#ifndef PAIREDENDREAD
#define PAIREDENDREAD

#include<cassert>
#include<iostream>
#include<string>

#include "Read.h"
#include "SingleRead.h"

class PairedEndRead : public Read {
 public:
  PairedEndRead() : mate1(), mate2() {}
  PairedEndRead(const SingleRead& mate1, const SingleRead& mate2) {
	  this->mate1 = mate1;
	  this->mate2 = mate2;
	  this->name = mate1.getName();

	  calc_lq();
  }

  bool read(int argc, std::istream* argv[], int flags = 7);
  void write(int argc, std::ostream* argv[]);

  const SingleRead& getMate1() const { return mate1; }
  const SingleRead& getMate2() const { return mate2; }
  const SingleRead& getMate(int i) const {
	  if (i == 1) return mate1;
	  else return mate2;
  }

 private:
  SingleRead mate1, mate2;

  void calc_lq();
};

bool PairedEndRead::read(int argc, std::istream* argv[], int flags) {
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

void PairedEndRead::write(int argc, std::ostream *argv[]) {
	std::ostream *outMate1[1], *outMate2[1];

	assert(argc == 2);
	outMate1[0] = argv[0]; outMate2[0] = argv[1];
	mate1.write(1, outMate1);
	mate2.write(1, outMate2);
}

void PairedEndRead::calc_lq() {
	low_quality = mate1.isLowQuality() && mate2.isLowQuality();
}

#endif
