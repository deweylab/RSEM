#ifndef SINGLEREAD
#define SINGLEREAD

#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cassert>
#include<iostream>
#include<string>

#include "utils.h"
#include "Read.h"

class SingleRead : public Read {
public:
	SingleRead() { readseq = ""; len = 0; }
	SingleRead(const std::string& name, const std::string& readseq) {
		this->name = name;
		this->readseq = readseq;
		this->len = readseq.length();
	}

	bool read(int argc, std::istream* argv[], int flags = 7);
	void write(int argc, std::ostream* argv[]);

	const int getReadLength() const { return len; /*readseq.length();*/ } // If need memory and .length() are guaranteed O(1), use statement in /* */
	const std::string& getReadSeq() const { return readseq; }

	void calc_lq(bool, int); // calculate if this read is low quality. Without calling this function, isLowQuality() will always be false

private:
	int len; // read length
	std::string readseq; // read sequence
};

//If return false, you should not trust the value of any member
bool SingleRead::read(int argc, std::istream* argv[], int flags) {
	std::string line;

	assert(argc == 1);
	if (!getline((*argv[0]), line)) return false;
	if (line[0] != '>') { fprintf(stderr, "Read file does not look like a FASTA file!"); exit(-1); }
	name = "";
	if (flags & 4) { name = line.substr(1); }
	if (!getline((*argv[0]), readseq)) return false;
	len = readseq.length(); // set read length
	if (!(flags & 1)) { readseq = ""; }

	return true;
}

void SingleRead::write(int argc, std::ostream* argv[]) {
	assert(argc == 1);
	(*argv[0])<<">"<<name<<std::endl<<readseq<<std::endl;
}

//calculate if this read is low quality
void SingleRead::calc_lq(bool hasPolyA, int seedLen) {
	low_quality = false;
	if (len < seedLen) { low_quality = true; return; }

	// if no polyA, no need to do the following calculation
	if (!hasPolyA) return;

	assert(readseq != "");

	int numA = 0, numT = 0, numAO = 0, numTO = 0; // numAO : number of A in overlap seed region
	int threshold_1, threshold_2;

	threshold_1 = int(0.9 * len - 1.5 * sqrt(len * 1.0) + 0.5);
	threshold_2 = (OLEN - 1) / 2 + 1;
	for (int i = 0; i < len; i++) {
		if (readseq[i] == 'A') {
			++numA;
			if (i < OLEN) ++numAO;
		}
		if (readseq[i] == 'T') {
			++numT;
			if (i >= len - OLEN) ++numTO;
		}
	}
  
	if (numA >= threshold_1) {
		low_quality = (numAO >= threshold_2);
	}
	else if (numT >= threshold_1) {
		low_quality = (numTO >= threshold_2);
	}
	else low_quality = false;
}

#endif
