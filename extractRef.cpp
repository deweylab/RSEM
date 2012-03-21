#include<cstdio>
#include<cstring>
#include<cctype>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<map>
#include<vector>
#include<algorithm>

#include "utils.h"
#include "GTFItem.h"
#include "Transcript.h"
#include "Transcripts.h"

using namespace std;

struct ChrInfo {
	string name;
	size_t len;

	ChrInfo(const string& name, size_t len) {
		this->name = name;
		this->len = len;
	}

	bool operator< (const ChrInfo& o) const {
		return name < o.name;
	}
};

int M;

vector<GTFItem> items;
vector<string> seqs;
vector<int> starts; // used to generate .grp
map<string, vector<int> > sn2tr; // map from seqname to transcripts
map<string, vector<int> >::iterator iter;
vector<ChrInfo> chrvec;

Transcripts transcripts;

char groupF[STRLEN], tiF[STRLEN], refFastaF[STRLEN];
char chromListF[STRLEN];

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

bool buildTranscript(int sp, int ep) {
	int cur_s, cur_e; // current_start, current_end
	vector<Interval> vec;

	string transcript_id = items[sp].getTranscriptID();
	string gene_id = items[sp].getGeneID();
	char strand = items[sp].getStrand();
	string seqname = items[sp].getSeqName();
	string left = items[sp].getLeft();

	vec.clear();
	cur_s = cur_e = -1;
	for (int i = sp; i <= ep; i++) {
		int start = items[i].getStart();
		int end = items[i].getEnd();

		if (strand != items[i].getStrand()) {
		  fprintf(stderr, "According to the GTF file given, a transcript has exons from different orientations!\n");
		  exit(-1);
		}
		if (seqname != items[i].getSeqName()) {
		  fprintf(stderr, "According to the GTF file given, a transcript has exons on multiple chromosomes!\n");
		  exit(-1);
		}

		if (cur_e + 1 < start) {
			if (cur_s > 0) vec.push_back(Interval(cur_s, cur_e));
			cur_s = start;
		}
		cur_e = (cur_e < end ? end : cur_e);
	}
	if (cur_s > 0) vec.push_back(Interval(cur_s, cur_e));

	transcripts.add(Transcript(transcript_id, gene_id, seqname, strand, vec, left));

	return true;
}

void parse_gtf_file(char* gtfF) {
	ifstream fin(gtfF);
	string line, curgid, tid, gid; //  curgid: current gene id;
	GTFItem item;

	if (!fin.is_open()) { fprintf(stderr, "Cannot open %s! It may not exist.\n", gtfF); exit(-1); }

	int cnt = 0;

	items.clear();
 	while (getline(fin, line)) {
 		if (line[0] == '#') continue; // if this line is comment, jump it
 		item.parse(line);
 		string feature = item.getFeature();
 		if (feature == "exon") {
 			if (item.getStart() > item.getEnd()) {
 				fprintf(stderr, "Warning: exon's start position is larger than its end position! This exon is discarded.\n");
 				fprintf(stderr, "\t%s\n\n", line.c_str());
 			}
 			else if (item.getStart() < 1) {
 				fprintf(stderr, "Warning: exon's start position is less than 1! This exon is discarded.\n");
 				fprintf(stderr, "\t%s\n\n", line.c_str());
 			}
 			else {
 		 		if (hasMappingFile) {
 		 			tid = item.getTranscriptID();
					mi_iter = mi_table.find(tid);
					if (mi_iter == mi_table.end()) {
					  fprintf(stderr, "Mapping Info is not correct, cannot find %s's gene_id!\n", tid.c_str());
					  exit(-1);
					}
					//assert(iter != table.end());
					gid = mi_iter->second;
					item.setGeneID(gid);
 		 		}
 				items.push_back(item);
 			}
 		}

 		++cnt;
 		if (verbose && cnt % 200000 == 0) { printf("Parsed %d lines\n", cnt); }
	}
	fin.close();

	sort(items.begin(), items.end());

	starts.clear();
	sn2tr.clear();
	curgid = "";

	int sp = 0, ep; // start pointer, end pointer
	int nItems = items.size();

	while (sp < nItems) {
		tid = items[sp].getTranscriptID();
		gid = items[sp].getGeneID();

		ep = sp + 1;
		while (ep < nItems && items[ep].getTranscriptID() == tid) ep++;
		ep--;

		buildTranscript(sp, ep);

		int sid = transcripts.getM();
		const Transcript& transcript = transcripts.getTranscriptAt(sid);

		if (curgid != gid) {
			starts.push_back(sid);
			curgid = gid;
		}
		iter = sn2tr.find(transcript.getSeqName());
		if (iter == sn2tr.end()) {
			vector<int> vec(1, sid);
			sn2tr[transcript.getSeqName()] = vec;
		}
		else {
			iter->second.push_back(sid);
		}

		sp = ep + 1;
	}

	M = transcripts.getM();
	starts.push_back(M + 1);
	items.clear();

	if (M < 1) {
		fprintf(stderr, "Number of transcripts in the reference is less than 1!\n");
		exit(-1);
	}

	if (verbose) { printf("Parsing gtf File is done!\n"); }
}

char check(char c) {
	if (!isalpha(c)) { fprintf(stderr, "Sequence contains unknown letter '%c'!\n", c); exit(-1); }
	//assert(isalpha(c));
	if (isupper(c) && c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'N';
	if (islower(c) && c != 'a' && c != 'c' && c != 'g' && c != 't') c = 'n';
	return c;
}

void writeResults(char* refName) {
	int s;
	ofstream fout;

	sprintf(groupF, "%s.grp", refName);
	sprintf(tiF, "%s.ti", refName);
	sprintf(refFastaF, "%s.transcripts.fa", refName);
	sprintf(chromListF, "%s.chrlist", refName);


	fout.open(groupF);
	s = starts.size();
	for (int i = 0; i < s; i++) fout<<starts[i]<<endl;
	fout.close();
	if (verbose) { printf("Group File is generated!\n"); }

	transcripts.writeTo(tiF);
	if (verbose) { printf("Transcript Information File is generated!\n"); }

	fout.open(chromListF);
	s = chrvec.size();
	for (int i = 0; i < s; i++) {
		fout<<chrvec[i].name<<'\t'<<chrvec[i].len<<endl;
	}
	fout.close();
	if (verbose) { printf("Chromosome List File is generated!\n"); }

	fout.open(refFastaF);
	for (int i = 1; i <= M; i++) {
		fout<<">"<<transcripts.getTranscriptAt(i).getTranscriptID()<<endl;
		fout<<seqs[i]<<endl;
	}
	fout.close();
	if (verbose) { printf("Extracted Sequences File is generated!\n"); }
}

int main(int argc, char* argv[]) {
  if (argc < 6 || ((hasMappingFile = atoi(argv[4])) && argc < 7)) {
		printf("Usage: rsem-extract-reference-transcripts refName quiet gtfF hasMappingFile [mappingFile] chromosome_file_1 [chromosome_file_2 ...]\n");
		exit(-1);
	}

	verbose = !atoi(argv[2]);
	if (hasMappingFile) {
		loadMappingInfo(argv[5]);
	}
	parse_gtf_file(argv[3]);

	ifstream fin;
	string line, gseq, seqname;
	void* pt;

	chrvec.clear();

	seqs.clear();
	seqs.resize(M + 1, "");
	int start = hasMappingFile ? 6 : 5;
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

		    size_t len = gseq.length();
		    assert(len > 0);
		    for (size_t j = 0; j < len; j++) gseq[j] = check(gseq[j]);

		    iter = sn2tr.find(seqname);
		    if (iter == sn2tr.end()) continue;

		    chrvec.push_back(ChrInfo(seqname, len));

		    vector<int>& vec = iter->second;
		    int s = vec.size();
		    for (int j = 0; j < s; j++) {
		    	assert(vec[j] > 0 && vec[j] <= M);
		    	transcripts.getTranscriptAt(vec[j]).extractSeq(gseq, seqs[vec[j]]);
		    }
		}
		fin.close();

		if (verbose) { printf("%s is processed!\n", argv[i]); }
	}

	for (int i = 1; i <= M; i++) {
		if (seqs[i] == "") {
			const Transcript& transcript = transcripts.getTranscriptAt(i);
			fprintf(stderr, "Cannot extract transcript %s's sequence from chromosome %s, whose information might not be provided! Please check if the chromosome directory is set correctly or the list of chromosome files is complete.\n", \
					transcript.getTranscriptID().c_str(), transcript.getGeneID().c_str());
			exit(-1);
		}
	}

	sort(chrvec.begin(), chrvec.end());

	if (verbose) { printf("Extracting sequences is done!\n"); }

	writeResults(argv[1]);

	return 0;
}
