#ifndef SINGLEREADQ
#define SINGLEREADQ

#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cassert>
#include<string>
#include<iostream>

#include "utils.h"
#include "Read.h"

class SingleReadQ : public Read {
public:
	SingleReadQ() { readseq = qscore = ""; len = 0; }
	SingleReadQ(const std::string& name, const std::string& readseq, const std::string& qscore) {
		this->name = name;
		this->readseq = readseq;
		this->qscore = qscore;
		this->len = readseq.length();
	}

	bool read(int argc, std::istream* argv[], int flags = 7);
	void write(int argc, std::ostream* argv[]);

	int getReadLength() const { return len; }
	const std::string& getReadSeq() const { return readseq; }
	const std::string& getQScore() const { return qscore; }

	void calc_lq(bool, int); // calculate if this read is low quality. Without calling this function, isLowQuality() will always be false

private:
	int len; // read length
	std::string readseq, qscore; // qscore : quality scores
};

bool SingleReadQ::read(int argc, std::istream* argv[], int flags) {
	std::string line;

	assert(argc == 1);
	if (!getline((*argv[0]), line)) return false;
	if (line[0] != '@') { fprintf(stderr, "Read file does not look like a FASTQ file!\n"); exit(-1); }
	name = "";
	if (flags & 4) { name = line.substr(1); }
	if (!getline((*argv[0]), readseq)) return false;
	len = readseq.length();
	if (!(flags & 1)) { readseq = ""; }
	if (!getline((*argv[0]), line)) return false;
	if (line[0] != '+') { fprintf(stderr, "Read file does not look like a FASTQ file!\n"); exit(-1); }
	if (!getline((*argv[0]), qscore)) return false;
	if (!(flags & 2)) { qscore = ""; }

	return true;
}

void SingleReadQ::write(int argc, std::ostream* argv[]) {
	assert(argc == 1);
	(*argv[0])<<"@"<<name<<std::endl<<readseq<<std::endl<<"+\n"<<qscore<<std::endl;
}

//calculate if this read is low quality
void SingleReadQ::calc_lq(bool hasPolyA, int seedLen) {
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
