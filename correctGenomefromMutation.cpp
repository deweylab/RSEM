#include <cstdio>
#include <cctype>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>

#include "utils.h"

using namespace std;

bool verbose = true;

struct MAFLine {

	MAFLine(const string& chr, const string& event, int pos, const string& ref, const string& allele1, const string& allele2) : chr(chr), event(event), pos(pos), ref(ref), allele1(allele1), allele2(allele2) {}

	string chr, event;
	int pos; // 0-based
	string ref, allele1, allele2;
};

typedef unordered_map<string, vector<MAFLine> > ChrMap;
typedef ChrMap::iterator ChrMapIter;

ChrMap subs, ins_del;

vector<char> chrseq, newseq;

void loadMAF(const char* inpF, ChrMap& chrMap) {
	ifstream fin(inpF);
	string line;
	vector<string> fields;
	ChrMapIter it;

	chrMap.clear();

	getline(fin, line);
	while (getline(fin, line)) {
		split(line, '\t', fields);
		it = chrMap.find(fields[4]);
		if (it == chrMap.end()) it = chrMap.emplace(fields[4], vector<MAFLine>()).first;
		it->second.emplace_back(fields[4], fields[9], stoi(fields[5]) - 1, fields[10], fields[11], fields[12]);
	}

	fin.close();

	if (verbose) printf("%s is loaded.\n", inpF);
}

void correctGenome(const char* inpF, const char* outF) {
	ifstream fin(inpF);
	ofstream fout(outF);

	string header, line, chrn;
	vector<string> fields;
	ChrMapIter it1, it2;	

	getline(fin, line);

	do {
		assert(line[0] == '>');
		header = line;
		split(header, ' ', fields);
		chrn = fields[1];

		chrseq.clear();
		while (getline(fin, line) && line[0] != '>') {
			for (auto&& c : line) if (isalpha(c)) chrseq.push_back(toupper(c));				
		}

		// check consistency
		it1 = subs.find(chrn);
		if (it1 != subs.end()) {
			for (auto&& mafl : it1->second) {
				for (int i = 0; i < mafl.ref.length(); ++i) assert(mafl.ref[i] == chrseq[mafl.pos + i]);
			}
		}

		// check consistency
		it2 = ins_del.find(chrn);
		if (it2 != ins_del.end()) {
			for (auto&& mafl : it2->second) if (mafl.event == "DEL") {
				for (int i = 0; i < mafl.ref.length(); ++i) assert(mafl.ref[i] == chrseq[mafl.pos + i]);
			}
		}

		if (it1 != subs.end()) {
			for (auto&& mafl : it1->second) 
				for (int i = 0; i < mafl.allele1.length(); ++i) chrseq[mafl.pos + i] = toupper(mafl.allele1[i]);
		}


		if (it2 != subs.end()) {
			newseq.clear();
			int pos = 0;
			for (auto&& mafl : it2->second) {
				while (pos < mafl.pos) newseq.push_back(chrseq[pos]), ++pos;
				if (mafl.event == "INS") {
					newseq.push_back(chrseq[pos]);
					for (auto&& c : mafl.allele1) newseq.push_back(toupper(c));
					++pos;
				}
				else pos += mafl.ref.length();
			}
			while (pos < chrseq.size()) newseq.push_back(chrseq[pos]), ++pos;
		}
		else newseq = chrseq; 

		fout<< header<< endl;
		for (int i = 0; i < newseq.size(); ++i) {
			fout<< newseq[i]; 
			if ((i + 1) % 80 == 0) fout<< endl;
		}
		if (newseq.size() % 80 != 0) fout<< endl;

	} while (fin.good());


	fin.close();
	fout.close();

	if (verbose) printf("correctGenome finished.\n");
}

int main(int argc, char* argv[]) {
	if (argc != 5) { 
		printf("Usage: correctGenomefromMutation subs.maf ins_del.maf input_genome.fa output_genome.fa\n");
		exit(-1);
	}

	loadMAF(argv[1], subs);
	loadMAF(argv[2], ins_del);

	correctGenome(argv[3], argv[4]);

	return 0;
}