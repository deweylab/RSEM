/*
 * transcripts are numbered from 1. 0 is reserved for noise isoform
 */
#ifndef TRANSCRIPTS_H_
#define TRANSCRIPTS_H_

#include<cstdio>
#include<cstdlib>
#include<cassert>
#include<fstream>
#include<vector>
#include<algorithm>

#include "Transcript.h"


class Transcripts {
public:
	Transcripts(int type = 0) {
		M = 0; this->type = type;
		transcripts.clear();
		transcripts.push_back(Transcript());
	}

	int getM() { return M; }
	int getType() { return type; }

	const Transcript& getTranscriptAt(int pos) {
		assert(pos > 0 && pos <= M);
		return transcripts[pos];
	}

	void add(const Transcript& transcript) {
		transcripts.push_back(transcript);
		++M;
	}

	void sort() {
		std::sort(transcripts.begin(), transcripts.end());
	}

	void readFrom(const char*);
	void writeTo(const char*);

private:
	int M, type; // type 0 from genome , 1 standalone transcriptome
	std::vector<Transcript> transcripts;
};

void Transcripts::readFrom(const char* inpF) {
	std::string line;
	std::ifstream fin(inpF);

	if (!fin.is_open()) { fprintf(stderr, "Cannot open %s! It may not exist.\n", inpF); exit(-1); }

	fin>>M>>type;
	getline(fin, line);
	transcripts.clear();
	transcripts.resize(M + 1);
	for (int i = 1; i <= M; i++) {
		transcripts[i].read(fin);
	}
	fin.close();
}

void Transcripts::writeTo(const char* outF) {
	std::ofstream fout(outF);
	fout<<M<<" "<<type<<std::endl;
	for (int i = 1; i <= M; i++) {
		transcripts[i].write(fout);
	}
	fout.close();
}

#endif /* TRANSCRIPTS_H_ */
