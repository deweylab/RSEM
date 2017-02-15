#include <cstdio>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>

#include "utils.h"
#include "Transcript.h"
#include "Transcripts.h"

using namespace std;


struct GeneRecord {
	string gid, breakpoint;
	bool isLeft;
	char dir;
	int coord;

	GeneRecord(const string& gid, const string& breakpoint, bool isLeft, char dir, int coord) : gid(gid), breakpoint(breakpoint), isLeft(isLeft), dir(dir), coord(coord) {
	}
};

struct SeqType {
	string name, seq;

	SeqType(const string& name, const string& seq) : name(name), seq(seq) {}

	bool operator< (const SeqType& o) const {
		return seq < o.seq;
	}

	bool operator== (const SeqType& o) const {
		return seq == o.seq;
	}
};

char tiF[STRLEN], seqF[STRLEN];

vector<GeneRecord> records;
map<string, vector<int> > gidMap;

int M;
Transcripts transcripts;
vector<SeqType> seqs;
vector<set<SeqType> > fuseCands;

void parseInput(char* inpF) {
	ifstream fin(inpF);
	string line;
	int no;

	no = 0;
	records.clear();
	gidMap.clear();
	getline(fin, line);
	while (getline(fin, line)) {
		istringstream strin(line);

		string token, gid;
		char dir;
		int coord;
		size_t fr, to;

		strin>> token>> token>> token>> token;
		strin>> token;

		fr = token.find_first_of('^') + 1;
		to = token.find_last_of('.');
		gid = token.substr(fr, to - fr);

		strin>> token;
		fr = token.find_first_of(':') + 1;
		to = token.find_last_of(':');
		dir = token[to + 1];
		coord = stoi(token.substr(fr, to - fr));

		records.push_back(GeneRecord(gid, token, true, dir, coord));
		auto iter = gidMap.find(gid);
		if (iter == gidMap.end()) {
			gidMap.emplace(gid, vector<int>(1, no));
		}
		else iter->second.push_back(no);
		++no;


		strin>> token;
		fr = token.find_first_of('^') + 1;
		to = token.find_last_of('.');
		gid = token.substr(fr, to - fr);

		strin>> token;
		fr = token.find_first_of(':') + 1;
		to = token.find_last_of(':');
		dir = token[to + 1];
		coord = stoi(token.substr(fr, to - fr));

		records.push_back(GeneRecord(gid, token, false, dir, coord));
		iter = gidMap.find(gid);
		if (iter == gidMap.end()) {
			gidMap.emplace(gid, vector<int>(1, no));
		}
		else iter->second.push_back(no);
		++no;

	}
	fin.close();
	cout<< "parseInput is finished."<< endl;
}

void loadReference(char* reference_name) {
	sprintf(tiF, "%s.ti", reference_name);
	transcripts.readFrom(tiF);
	M = transcripts.getM();

	sprintf(seqF, "%s.transcripts.fa", reference_name);
	ifstream fin(seqF);
	string tag, seq;
	seqs.clear();
	seqs.emplace_back("", "");
	for (int i = 1; i <= M; ++i) {
		getline(fin, tag);
		getline(fin, seq);
		seqs.emplace_back(tag.substr(1), seq);
	}
	fin.close();

	cout<< "loadReference is finished."<< endl;
}

void enumerateFusedTranscripts(char* outF) {
	fuseCands.clear();
	fuseCands.resize(records.size());
	for (int i = 1; i <= M; ++i) {
		const Transcript& transcript = transcripts.getTranscriptAt(i);
		auto iter = gidMap.find(transcript.getGeneID());
		if (iter == gidMap.end()) continue;
		for (auto&& no : iter->second) {
			if (transcript.getStrand() != records[no].dir) continue;

			int len = 0, j;
			bool bingo = false;
			const std::vector<Interval>& exons = transcript.getStructure();
			if ((records[no].isLeft && records[no].dir == '+') || (!records[no].isLeft && records[no].dir == '-')) {
				for (j = 0; j < exons.size(); ++j) {
					len += exons[j].end - exons[j].start + 1;
					if (exons[j].end == records[no].coord) break;
				}
				bingo = j < exons.size();
			}
			else {
				for (j = exons.size() - 1; j >= 0; --j) {
					len += exons[j].end - exons[j].start + 1;
					if (exons[j].start == records[no].coord) break;
				}
				bingo = j >= 0;
			}

			if (bingo) {
				string name = seqs[i].name + ":" + records[no].breakpoint;
				string sequence = (records[no].isLeft ? seqs[i].seq.substr(0, len) : seqs[i].seq.substr(transcript.getLength() - len, len));

				fuseCands[no].emplace(name, sequence);
			}
		}
	}

	ofstream fout(outF);
	int nFuses = 0;
	set<string> seen;
	for (size_t i = 0; i < records.size(); i += 2) {
		for (auto&& left : fuseCands[i]) 
			for (auto && right : fuseCands[i + 1]) {
				string sequence = left.seq + right.seq;
				if (seen.find(sequence) == seen.end()) {
					seen.insert(sequence);
					fout<< ">Fuse^"<< left.name<< "^"<< right.name<< endl << sequence<< endl;
					++nFuses;
				}
			}
	}
	for (int i = 1; i <= M; ++i) 
		fout<< ">"<< seqs[i].name<< endl<< seqs[i].seq<< endl;
	fout.close();
	
	cout<< nFuses<< " candidates are found."<< endl;	
}


int main(int argc, char* argv[]) {
	if (argc != 4) {
		printf("Usage: genFusionTranscripts fusion_input.txt reference_name output.fa\n");
		exit(-1);
	}

	parseInput(argv[1]);
	loadReference(argv[2]);
	enumerateFusedTranscripts(argv[3]);

	return 0;
}
