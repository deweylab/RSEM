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

#include "utils.h"
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

	// used in shrinking the transcripts
	void setM(int M) { this->M = M; transcripts.resize(M + 1); } 
	
	void move(int from, int to) {
	  assert(from >= to);
	  if (from > to) transcripts[to] = transcripts[from];
	}
	
	int getType() { return type; }
	void setType(int type) { this->type = type; }

	bool isAlleleSpecific() { return type == 2; }

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

	void buildMappings(int, char**, const char* = NULL);

private:
	int M, type; // type 0 from genome, 1 standalone transcriptome, 2 allele-specific 
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

void Transcripts::buildMappings(int n_targets, char** target_name, const char* imdName) {
	std::map<std::string, int> dict;
	std::map<std::string, int>::iterator iter;
	std::vector<bool> appeared;

	general_assert(n_targets > 0, "The SAM/BAM file declares less than one reference sequence!");
	general_assert(n_targets <=  M, "The SAM/BAM file declares more reference sequences (" + itos(n_targets) + ") than RSEM knows (" + itos(M) + ")!");
	if (n_targets < M) fprintf(stderr, "Warning: The SAM/BAM file declares less reference sequences (%d) than RSEM knows (%d)! Please make sure that you aligned your reads against transcript sequences instead of genome.\n", n_targets, M);

	dict.clear();
	for (int i = 1; i <= M; i++) {
		const std::string& tid = isAlleleSpecific() ? transcripts[i].getSeqName() : transcripts[i].getTranscriptID();
		iter = dict.find(tid);
		general_assert(iter == dict.end(), "RSEM's indices might be corrupted, " + tid + " appears more than once!");
		dict[tid] = i;
	}

	e2i.assign(M + 1, 0);
	i2e.assign(M + 1, 0);
	appeared.assign(M + 1, false);
	for (int i = 0; i < n_targets; i++) {
		iter = dict.find(std::string(target_name[i]));
		general_assert(iter != dict.end(), "RSEM can not recognize reference sequence name " + cstrtos(target_name[i]) + "!");
		general_assert(iter->second > 0, "Reference sequence name " + cstrtos(target_name[i]) + " appears more than once in the SAM/BAM file!");
		e2i[i + 1] = iter->second;
		i2e[iter->second] = i + 1;
		iter->second = -1;
		appeared[e2i[i + 1]] = true;
	}

	if (imdName != NULL) {
	  char omitF[STRLEN];
	  sprintf(omitF, "%s.omit", imdName);
	  FILE *fo = fopen(omitF, "w");
	  for (int i = 1; i <= M; i++) 
	    if (!appeared[i]) fprintf(fo, "%d\n", i);
	  fclose(fo);
	}
}

#endif /* TRANSCRIPTS_H_ */
