#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<fstream>
#include<sstream>
#include<map>
#include<vector>

#include "utils.h"
#include "my_assert.h"
#include "Transcript.h"
#include "Transcripts.h"

using namespace std;

bool verbose = true;

int M;

map<string, string> name2seq;
map<string, string>::iterator iter;

Transcripts transcripts(1); // no genome, just transcript set
char groupF[STRLEN], tiF[STRLEN], refFastaF[STRLEN];
char gtF[STRLEN], taF[STRLEN]; // group info between gene and transcript, transcript and allele

int hasMappingFile;
char mappingFile[STRLEN];

map<string, string> mi_table; // mapping info table
map<string, string>::iterator mi_iter; //mapping info table's iterator

map<string, string> mi_table2; // allele_id to transcript_id
map<string, string>::iterator mi_iter2; // corresponding iterator

void loadMappingInfo(int file_type, char* mappingF) {
  ifstream fin(mappingF);
  string line, key, value, value2;

  general_assert(fin.is_open(), "Cannot open " + cstrtos(mappingF) + "! It may not exist.");

  switch(file_type) {
  case 1: 
    mi_table.clear();
    while (getline(fin, line)) {
      line = cleanStr(line);
      if (line[0] == '#') continue;
      istringstream strin(line);
      strin>>value>>key;
      mi_table[key] = value;
    }
    break;
  case 2: 
    mi_table.clear();
    mi_table2.clear();
    while (getline(fin, line)) {
      line = cleanStr(line);
      if (line[0] == '#') continue;
      istringstream strin(line);
      strin>> value>> value2>> key;
      mi_table[key] = value;
      mi_table2[key] = value2;
    }
    break;
  default: assert(false);
  }

  fin.close();
}

void writeResults(int option, char* refName) {
	ofstream fout;
	string cur_gene_id, cur_transcript_id, name;
	vector<int> gi, gt, ta;

	sprintf(tiF, "%s.ti", refName);
	transcripts.writeTo(tiF);
	if (verbose) { printf("Transcript Information File is generated!\n"); }

	cur_gene_id = ""; gi.clear(); 
	if (option == 2) { cur_transcript_id = ""; gt.clear(); ta.clear(); }
	for (int i = 1; i <= M; i++) {
		const Transcript& transcript = transcripts.getTranscriptAt(i);
		if (cur_gene_id != transcript.getGeneID()) {
		  gi.push_back(i);
		  if (option == 2) gt.push_back((int)ta.size());
		  cur_gene_id = transcript.getGeneID();
		}
		if ((option == 2) && (cur_transcript_id != transcript.getTranscriptID())) {
		    ta.push_back(i);
		    cur_transcript_id = transcript.getTranscriptID();
		}
	}
	gi.push_back(M + 1);
	if (option == 2) { gt.push_back((int)ta.size()); ta.push_back(M + 1); }

	sprintf(groupF, "%s.grp", refName);
	fout.open(groupF);
	for (int i = 0; i < (int)gi.size(); i++) fout<< gi[i]<< endl;
	fout.close();
	if (verbose) { printf("Group File is generated!\n"); }

	if (option == 2) {
	  sprintf(gtF, "%s.gt", refName);
	  fout.open(gtF);
	  for (int i = 0; i < (int)gt.size(); i++) fout<< gt[i]<< endl;
	  fout.close();
	  sprintf(taF, "%s.ta", refName);
	  fout.open(taF);
	  for (int i = 0; i < (int)ta.size(); i++) fout<< ta[i]<< endl;
	  fout.close();
	  if (verbose) { printf("Allele-specific group files are generated!\n"); }
	}

	sprintf(refFastaF, "%s.transcripts.fa", refName);
	fout.open(refFastaF);
	for (int i = 1; i <= M; i++) {
		name = transcripts.getTranscriptAt(i).getSeqName();
		iter = name2seq.find(name);
		general_assert(iter != name2seq.end(), "Cannot recognize sequence ID" + name + "!");
		fout<<">"<<name<<endl;
		fout<<iter->second<<endl;
	}
	fout.close();
	
	if (verbose) { 
	  printf("Extracted Sequences File is generated!\n"); 
	}
}

struct CursorPos {
  char *filename;
  int line_no, pos;
} cursor;

inline char check(char c) {
  general_assert(isalpha(c), "FASTA file " + cstrtos(cursor.filename) + " contains an unknown character, " + \
		 ctos(c) + " (ASCII code " + itos(c) + "), at line " + itos(cursor.line_no) + ", position " + itos(cursor.pos + 1) + "!");
  if (isupper(c) && c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'N';
  if (islower(c) && c != 'a' && c != 'c' && c != 'g' && c != 't') c = 'n';
  return c;
}

int main(int argc, char* argv[]) {
  if (argc < 5 || ((hasMappingFile = atoi(argv[3])) && argc < 6)) {
		printf("Usage: synthesisRef refName quiet hasMappingFile<0,no;1,yes;2,allele-specific> [mappingFile] reference_file_1 [reference_file_2 ...]\n");
		exit(-1);
	}

	verbose = !atoi(argv[2]);

	if (hasMappingFile) { loadMappingInfo(hasMappingFile, argv[4]); }

	// allele-specific
	if (hasMappingFile == 2) { transcripts.setType(2); }

	int start = hasMappingFile ? 5 : 4;

	ifstream fin;
	string line, gseq;
	string seqname, gene_id, transcript_id;
	int seqlen, len;
	
	vector<Interval> vec;

	M = 0;
	name2seq.clear();
	for (int i = start; i < argc; i++) {
		fin.open(argv[i]);
		general_assert(fin.is_open(), "Cannot open " + cstrtos(argv[i]) + "! It may not exist.");

		cursor.filename = argv[i]; cursor.line_no = cursor.pos = 0;
		
		getline(fin, line);
		while ((fin) && (line[0] == '>')) {
			istringstream strin(line.substr(1));
			strin>>seqname;
			++cursor.line_no;
			
			gseq = ""; seqlen = 0;
			while((getline(fin, line)) && (line[0] != '>')) {
			  ++cursor.line_no;
			  len = line.length();
			  for (cursor.pos = 0; cursor.pos < len; ++cursor.pos) line[cursor.pos] = check(line[cursor.pos]);
			  seqlen += len;
			  gseq += line;
			}
			assert(seqlen > 0);
			name2seq[seqname] = gseq;

			transcript_id = seqname;
			gene_id = seqname;

			if (hasMappingFile) {
			      mi_iter = mi_table.find(seqname);
			      general_assert(mi_iter != mi_table.end(), "Mapping Info is not correct, cannot find " + seqname + "'s gene_id!");
			      gene_id = mi_iter->second;
			      if (hasMappingFile == 2) {
				mi_iter2 = mi_table2.find(seqname);
				general_assert(mi_iter2 != mi_table2.end(), "Mapping Info is not correct, cannot find allele " + seqname + "'s transcript_id!");
				transcript_id = mi_iter2->second;
			      }
			}
			
			vec.clear();
			vec.push_back(Interval(1, seqlen));
			transcripts.add(Transcript(transcript_id, gene_id, seqname, '+', vec, ""));
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

	writeResults(hasMappingFile, argv[1]);

	return 0;
}
