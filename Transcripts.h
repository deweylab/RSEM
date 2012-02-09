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
#include<map>
#include<string>

#include "my_assert.h"
#include "Transcript.h"


class Transcripts {
public:
	Transcripts(int type = 0) {
		M = 0; this->type = type;
		transcripts.clear();
		transcripts.push_back(Transcript());

		e2i.clear(); i2e.clear();
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

	//Eid: external sid
	int getInternalSid(int eid) {
		assert(eid > 0 && eid <= M);
		return e2i[eid];
	}

	const Transcript& getTranscriptViaEid(int eid) {
		return transcripts[getInternalSid(eid)];
	}

	void buildMappings(int, char**);

private:
	int M, type; // type 0 from genome , 1 standalone transcriptome
	std::vector<Transcript> transcripts;

	std::vector<int> e2i, i2e; // external sid to internal sid, internal sid to external sid
};

void Transcripts::readFrom(const char* inpF) {
	std::string line;
	std::ifstream fin(inpF);

	if (!fin.is_open()) { fprintf(stderr, "Cannot open %s! It may not exist.\n", inpF); exit(-1); }

	fin>>M>>type;
	getline(fin, line);
	transcripts.assign(M + 1, Transcript());
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

void Transcripts::buildMappings(int n_targets, char** target_name) {
	std::map<std::string, int> dict;
	std::map<std::string, int>::iterator iter;

	general_assert(n_targets == M, "Number of transcripts does not match! Please check if the reads are aligned to a transcript set (instead of a genome)!");

	dict.clear();
	for (int i = 1; i <= M; i++) {
		const std::string& tid = transcripts[i].getTranscriptID();
		iter = dict.find(tid);
		assert(iter == dict.end());
		dict[tid] = i;
	}

	e2i.assign(M + 1, 0);
	i2e.assign(M + 1, 0);
	for (int i = 0; i < n_targets; i++) {
		iter = dict.find(std::string(target_name[i]));
		general_assert(iter != dict.end(), "RSEM can not recognize transcript " + cstrtos(target_name[i]) + "!");
		general_assert(iter->second > 0, "Reference sequence name " + cstrtos(target_name[i]) + " is duplicated!");
		e2i[i + 1] = iter->second;
		i2e[iter->second] = i + 1;
		iter->second = -1;
	}
}

#endif /* TRANSCRIPTS_H_ */
