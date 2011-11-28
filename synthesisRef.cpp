#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<fstream>
#include<sstream>
#include<map>

#include "utils.h"
#include "Transcript.h"
#include "Transcripts.h"

using namespace std;

int M;

map<string, string> name2seq;
map<string, string>::iterator iter;

Transcripts transcripts(1); // no genome, just transcript set
char groupF[STRLEN], tiF[STRLEN], refFastaF[STRLEN], chromListF[STRLEN];

bool hasMappingFile;
char mappingFile[STRLEN];

map<string, string> mi_table; // mapping info table
map<string, string>::iterator mi_iter; //mapping info table's iterator

void loadMappingInfo(char* mappingF) {
  ifstream fin(mappingF);
  string line, key, value;

  if (!fin.is_open()) {
	  fprintf(stderr, "Cannot open %s! It may not exist.\n", mappingF);
	  exit(-1);
  }

  mi_table.clear();
  while (getline(fin, line)) {
    line = cleanStr(line);
    if (line[0] == '#') continue;
    istringstream strin(line);
    strin>>value>>key;
    mi_table[key] = value;
  }

  fin.close();
}

char check(char c) {
	if (!isalpha(c)) { fprintf(stderr, "Sequence contains unknown letter '%c'!\n", c); exit(-1); }
	//assert(isalpha(c));
	if (isupper(c) && c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'N';
	if (islower(c) && c != 'a' && c != 'c' && c != 'g' && c != 't') c = 'n';
	return c;
}

void writeResults(char* refName) {
	ofstream fout, fout2;
	string cur_gene_id, name;

	sprintf(groupF, "%s.grp", refName);
	sprintf(tiF, "%s.ti", refName);
	sprintf(refFastaF, "%s.transcripts.fa", refName);
	sprintf(chromListF, "%s.chrlist", refName);

	transcripts.writeTo(tiF);
	if (verbose) { printf("Transcript Information File is generated!\n"); }

	fout.open(groupF);
	cur_gene_id = "";
	for (int i = 1; i <= M; i++) {
		const Transcript& transcript = transcripts.getTranscriptAt(i);
		if (cur_gene_id != transcript.getGeneID()) {
			fout<<i<<endl;
			cur_gene_id = transcript.getGeneID();
		}
	}
	fout<<M + 1<<endl;
	fout.close();
	if (verbose) { printf("Group File is generated!\n"); }

	fout2.open(chromListF);
	fout.open(refFastaF);
	for (int i = 1; i <= M; i++) {
		name = transcripts.getTranscriptAt(i).getTranscriptID();
		iter = name2seq.find(name);
		if (iter == name2seq.end()) {
			fprintf(stderr, "Cannot recognize transcript ID %s!\n", name.c_str());
			exit(-1);
		}
		fout<<">"<<name<<endl;
		fout<<iter->second<<endl;

		fout2<<name<<'\t'<<iter->second.length()<<endl;
	}
	fout.close();
	fout2.close();
	
	if (verbose) { 
	  printf("Chromosome List File is generated!\n"); 
	  printf("Extracted Sequences File is generated!\n"); 
	}
}

int main(int argc, char* argv[]) {
  if (argc < 5 || ((hasMappingFile = atoi(argv[3])) && argc < 6)) {
		printf("Usage: synthesisRef refName quiet hasMappingFile [mappingFile] reference_file_1 [reference_file_2 ...]\n");
		exit(-1);
	}

	verbose = !atoi(argv[2]);

	if (hasMappingFile) { loadMappingInfo(argv[4]); }

	int start = hasMappingFile ? 5 : 4;

	ifstream fin;
	string line, gseq;
	string seqname, gene_id;
	void* pt;

	vector<Interval> vec;

	M = 0;
	name2seq.clear();
	for (int i = start; i < argc; i++) {
		fin.open(argv[i]);
		if (!fin.is_open()) { fprintf(stderr, "Cannot open %s! It may not exist.\n", argv[i]); exit(-1); }
		pt = getline(fin, line);
		while (pt != 0 && line[0] == '>') {
			istringstream strin(line.substr(1));
			strin>>seqname;

			gseq = "";
			while((pt = getline(fin, line)) && line[0] != '>') {
			    gseq += line;
			}

			int len = gseq.length();
			assert(len > 0);
			for (int j = 0; j < len; j++) gseq[j] = check(gseq[j]);

			name2seq[seqname] = gseq;

			if (hasMappingFile) {
			      mi_iter = mi_table.find(seqname);
			      if (mi_iter == mi_table.end()) {
			    	  fprintf(stderr, "Mapping Info is not correct, cannot find %s's gene_id!\n", seqname.c_str());
			    	  exit(-1);
			      }
			      //assert(iter != table.end());
			      gene_id = mi_iter->second;
			}
			else gene_id = seqname;

			vec.clear();
			vec.push_back(Interval(1, len));
			transcripts.add(Transcript(seqname, gene_id, seqname, '+', vec, ""));
			++M;

			if (verbose && M % 1000000 == 0) { printf("%d sequences are processed!\n", M); }
		}
		fin.close();
	}

	if (M < 1) {
		fprintf(stderr, "Number of transcripts in the reference is less than 1!\n");
		exit(-1);
	}

	assert(M == transcripts.getM());
	transcripts.sort();

	writeResults(argv[1]);

	return 0;
}
