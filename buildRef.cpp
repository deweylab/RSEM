/* Copyright (c) 2017
   Bo Li (The Broad Institute of MIT and Harvard)
   libo@broadinstitute.org

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.   

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA
*/

#include <cstdio>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>

#include "utils.h"
#include "my_assert.h"
#include "GTFItem.h"
#include "Transcript.hpp"
#include "Transcripts.hpp"
#include "Refs.hpp"

using namespace std;

bool verbose = true; // define verbose

// hasGTF: if has GTF file ; hasMapping: if has mapping file; isAllele: if allele-specific; rmdup : if remove duplicate sequences; appendPolyA : if append polyA tails; n2g_idx : if generate n2g index
bool hasGTF, hasMapping, isAllele, rmdup, appendPolyA, n2g_idx; 

int num_files; // Number of reference input files
int mappingPos, file_pos; // position in argv vector for mapping file and reference files
int polyALen; // length of poly(A)s appended

char gtfF[STRLEN];
char tiF[STRLEN], refFastaF[STRLEN];
char transListF[STRLEN], chromListF[STRLEN];
char groupF[STRLEN], gtF[STRLEN], taF[STRLEN]; 
char n2g_idxF[STRLEN];

// Mapping file related data structures
map<string, string> mi_table; // mapping info table
map<string, string>::iterator mi_iter; //mapping info table's iterator
map<string, string> mi_table2; // allele_id to transcript_id
map<string, string>::iterator mi_iter2; // corresponding iterator



int M;
Transcripts transcripts;
Refs refs;
vector<Interval> vec;

vector<GTFItem> items;
map<string, vector<int> > sn2tr; // map from chromosome name to transcripts
map<string, vector<int> >::iterator sn2tr_iter;






// for GTF extracting
vector<string> seqs;
// for transcripts input
map<string, string> name2seq;
map<string, string>::iterator n2s_iter;


// test if we should skip this line
inline bool skip(const string& line) {
	size_t pos = 0, len = line.length();
	while (pos < len && isspace(line[pos])) ++pos;
	return pos >= len || line[pos] == '#'; // skipping if empty line or commented line
}


// chromosome info

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

vector<ChrInfo> chrvec;



void loadMappingInfo(char* mappingF, bool isAllele) {
	ifstream fin(mappingF);
	string line, key, value, value2;

	general_assert(fin.is_open(), "Cannot open " + cstrtos(mappingF) + "! It may not exist.");

	if (!isAllele) {
		mi_table.clear();
		while (getline(fin, line)) {
			if (skip(line)) continue;
			istringstream strin(line);
			strin>> value>> key; // gene_id transcript_id
			mi_table[key] = value;
		}
	}
	else {
		mi_table.clear();
		mi_table2.clear();
		while (getline(fin, line)) {
			if (skip(line)) continue;
			istringstream strin(line);
			strin>> value>> value2>> key; // gene_id transcript_id allele_id
			mi_table[key] = value;
			mi_table2[key] = value2;
		}
	}

	fin.close();
}

// Extracting transcript sequences from the genome



void buildTranscript(int sp, int ep) {
	int cur_s, cur_e; // current_start, current_end
	
	string transcript_id = items[sp].getTranscriptID();
	string gene_id = items[sp].getGeneID();
	string gene_name = "", transcript_name = "";
	
	char strand = items[sp].getStrand();
	string seqname = items[sp].getSeqName();
	string left = items[sp].getLeft();
	
	vec.clear();
	cur_s = cur_e = -1;
	for (int i = sp; i <= ep; ++i) {
		int start = items[i].getStart();
		int end = items[i].getEnd();
		
		general_assert(strand == items[i].getStrand(), "According to the GTF file given, a transcript has exons from different orientations!");
		general_assert(seqname == items[i].getSeqName(), "According to the GTF file given, a transcript has exons on multiple chromosomes!");

		if (items[i].getGeneName() != "") {
			if (gene_name == "") gene_name = items[i].getGeneName();
			else general_assert(gene_name == items[i].getGeneName(), "A transcript is associated with multiple gene names!");
		}
		if (items[i].getTranscriptName() != "") {
			if (transcript_name == "") transcript_name = items[i].getTranscriptName();
			else general_assert(transcript_name == items[i].getTranscriptName(), "A transcript is associated with multiple transcript names!");
		}
		
		if (cur_e + 1 < start) {
			if (cur_s > 0) vec.push_back(Interval(cur_s, cur_e));
			cur_s = start;
		}
		cur_e = (cur_e < end ? end : cur_e);
	}
	if (cur_s > 0) vec.push_back(Interval(cur_s, cur_e));

	//  if (gene_name != "") gene_id += "_" + gene_name;
	//  if (transcript_name != "") transcript_id += "_" + transcript_name;
	
	transcripts.add(Transcript(transcript_id, gene_id, seqname, strand, vec, left));
}

void parse_gtf_file(char* gtfF) {
	ifstream fin(gtfF);
	string line, tid, gid;
	GTFItem item;

	general_assert(fin.is_open(), "Cannot open " + cstrtos(gtfF) + "! It may not exist.");

	int cnt = 0;
	
	items.clear();
	while (getline(fin, line)) {
		if (skip(line)) continue;
		item.parse(line);
		string feature = item.getFeature();
		if (feature == "exon") {
			if (item.getStart() > item.getEnd()) {
	printf("Warning: exon's start position is larger than its end position! This exon is discarded.\n");
	printf("\t%s\n\n", line.c_str());
			}
			else if (item.getStart() < 1) {
	printf("Warning: exon's start position is less than 1! This exon is discarded.\n");
	printf("\t%s\n\n", line.c_str());
			}
			else {
	item.parseAttributes(line);
	if (mappingType > 0) {
		tid = item.getTranscriptID();
		mi_iter = mi_table.find(tid);
		general_assert(mi_iter != mi_table.end(), "Mapping Info is not correct, cannot find " + tid + "'s gene_id!");
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
	
	int sp = 0, ep; // start pointer, end pointer
	int nItems = items.size();
	
	sn2tr.clear();
	while (sp < nItems) {
		tid = items[sp].getTranscriptID();
		
		ep = sp + 1;
		while (ep < nItems && items[ep].getTranscriptID() == tid) ++ep;
		--ep;
		
		buildTranscript(sp, ep);
		
		int sid = transcripts.getM();
		const Transcript& transcript = transcripts.getTranscriptAt(sid);
		
		sn2tr_iter = sn2tr.find(transcript.getSeqName());
		if (sn2tr_iter == sn2tr.end()) {
			vector<int> vec(1, sid);
			sn2tr[transcript.getSeqName()] = vec;
		}
		else {
			sn2tr_iter->second.push_back(sid);
		}
		
		sp = ep + 1;
	}
	
	items.clear();
	
	if (verbose) { printf("Parsing GTF File is done!\n"); }
}



inline string n2g(const string& seq) {
	string newseq = seq;
	int len = newseq.length();
	
	for (int i = 0; i < len; ++i)
		if (newseq[i] == 'N') newseq[i] = 'G';
	
	return newseq;
}

void writeToDisk(char* refName) {
	ofstream fout;

	sprintf(tiF, "%s.ti", refName);
	transcripts.writeTo(tiF);
	if (verbose) { printf("Transcript Information File is generated!\n"); }
	
	sprintf(refFastaF, "%s.transcripts.fa", refName);
	refs.writeTo(refFastaF);

	sprintf(transListF, "%s.translist", refName);
	refs.writeTransListTo(transListF);

	sprintf(chromListF, "%s.chrlist", refName);
	fout.open(chromListF);
	for (int i = 0; i < (int)chrvec.size(); ++i)
		fout<< chrvec[i].name<< '\t'<< chrvec[i].len<< endl;
	fout.close();
	if (verbose) { printf("Chromosome List File is generated!\n"); }
	
	string cur_gene_id, cur_transcript_id, name;
	vector<int> gi, gt, ta;

	cur_gene_id = ""; gi.clear(); 
	if (mappingType == 2) { cur_transcript_id = ""; gt.clear(); ta.clear(); }
	for (int i = 1; i <= M; ++i) {
		const Transcript& transcript = transcripts.getTranscriptAt(i);
		if (cur_gene_id != transcript.getGeneID()) {
			gi.push_back(i);
			if (mappingType == 2) gt.push_back((int)ta.size());
			cur_gene_id = transcript.getGeneID();
		}
		if ((mappingType == 2) && (cur_transcript_id != transcript.getTranscriptID())) {
			ta.push_back(i);
			cur_transcript_id = transcript.getTranscriptID();
		}
	}
	
	gi.push_back(M + 1);
	if (mappingType == 2) { gt.push_back((int)ta.size()); ta.push_back(M + 1); }

	sprintf(groupF, "%s.grp", refName);
	fout.open(groupF);
	for (int i = 0; i < (int)gi.size(); ++i) fout<< gi[i]<< endl;
	fout.close();
	if (verbose) { printf("Group File is generated!\n"); }

	if (mappingType == 2) {
		sprintf(gtF, "%s.gt", refName);
		fout.open(gtF);
		for (int i = 0; i < (int)gt.size(); ++i) fout<< gt[i]<< endl;
		fout.close();
		sprintf(taF, "%s.ta", refName);
		fout.open(taF);
		for (int i = 0; i < (int)ta.size(); ++i) fout<< ta[i]<< endl;
		fout.close();
		if (verbose) { printf("Allele-specific group files are generated!\n"); }
	}

	if (n2g_idx) {
		sprintf(n2g_idxF, "%s.n2g.idx.fa", refName);
		fout.open(n2g_idxF);
		for (int i = 1; i <= M; ++i) 
			fout<< '>'<< refs.getRef(i)->getName()<< endl<< n2g(refs.getRef(i)->getSeq())<< endl;
		fout.close();
		if (verbose) printf("%s is generated!\n", n2g_idxF);
	}
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		printf("Usage: PROBer-build-reference refName [--gtf gtfF] [--mapping mappingF] [--allele-specific] [--remove-duplicates] [--polyA-length length] [--n2g-index] [-q] [--files num_of_files file_1 file_2 ...]\n");
		exit(-1);
	}

	hasGTF = false;
	hasMapping = false;
	isAllele = false;
	rmdup = false;
	appendPolyA = false;
	n2g_idx = false;
	
	int argpos = 2;
	while (argpos < argc) {
		if (!strcmp(argv[argpos], "--gtf")) {
			hasGTF = true;
			strcpy(gtfF, argv[++argpos]);
		}
		if (!strcmp(argv[argpos], "--mapping")) {
			hasMapping = true;
			mappingPos = ++argpos;
		}
		if (!strcmp(argv[argpos], "--allele-specific")) isAllele = true;
		if (!strcmp(argv[argpos], "--remove-duplicates")) rmdup = true;
		if (!strcmp(argv[argpos], "--polyA-length")) {
			appendPolyA = true;
			polyALen = atoi(argv[++argpos]);
		}
		if (!strcmp(argv[argpos], "--n2g-index")) n2g_idx = true;
		if (!strcmp(argv[argpos], "-q")) verbose = false;
		if (!strcmp(argv[argpos], "--files")) {
			num_files = atoi(argv[++argpos]);
			file_pos = argpos + 1; // the position in argv for the first file
			argpos += num_files;
		}
		++argpos;
	}

	if (hasMapping) loadMappingInfo(argv[mappingPos], isAllele);

	ifstream fin;
	string line, gseq, tseq; // gseq, genomic sequence; tseq, transcript sequence
	string seqname, gene_id, transcript_id;
	
	if (hasGTF) {
		transcripts.setType(0);
		general_assert(!isAllele, "RSEM could not extract allele-specific transcript sequences from a genome!");
		parse_gtf_file(gtfF);

		M = transcripts.getM();
		general_assert(M > 0, "The reference contains no transcripts!");
		seqs.assign(M + 1, "");
		
		chrvec.clear();
		
		for (int i = 0; i < num_files; ++i, ++file_pos) {
			fin.open(argv[file_pos]);
			general_assert(fin.is_open(), "Cannot open " + cstrtos(argv[file_pos]) + "! It may not exist.");
			getline(fin, line);
			while ((fin) && (line[0] == '>')) {
	istringstream strin(line.substr(1));
	strin>>seqname;
	
	gseq = "";
	while((getline(fin, line)) && (line[0] != '>')) {
		gseq += line;
	}
	assert(gseq.length() > 0);
			
	sn2tr_iter = sn2tr.find(seqname);
	if (sn2tr_iter == sn2tr.end()) continue;
	
	chrvec.push_back(ChrInfo(seqname, gseq.length()));
	
	vector<int>& vec = sn2tr_iter->second;
	int s = vec.size();
	for (int j = 0; j < s; ++j) {
		assert(vec[j] > 0 && vec[j] <= M);
		transcripts.getTranscriptAt(vec[j]).extractSeq(gseq, seqs[vec[j]]);
	}
			}
			fin.close();

			if (verbose) { printf("%s is processed!\n", argv[file_pos]); } 
		}
		
		sort(chrvec.begin(), chrvec.end());

		// Shrink and build up Refs
		int curp = 0;
		for (int i = 1; i <= M; ++i) {
			const Transcript& transcript = transcripts.getTranscriptAt(i);
			if (seqs[i] == "") 
	printf("Warning: Cannot extract transcript %s because the chromosome it locates -- %s -- is absent!\n", transcript.getTranscriptID().c_str(), transcript.getSeqName().c_str());
			else {
	refs.addRef(transcript.getTranscriptID(), seqs[i]); // insert RefSeqs
	++curp;
	transcripts.move(i, curp);
			}
		}
		printf("%d transcripts are extracted and %d transcripts are omitted.\n", curp, M - curp);
		
		transcripts.setM(curp);
		M = transcripts.getM();
		general_assert(M > 0, "The reference contains no transcripts!");
		assert(refs.getM() == M);
	}
	else {
		transcripts.setType(mappingType != 2 ? 1 : 2);
		
		M = 0;
		name2seq.clear();
		for (int i = 0; i < num_files; ++i, ++file_pos) {
			fin.open(argv[file_pos]);
			general_assert(fin.is_open(), "Cannot open " + cstrtos(argv[file_pos]) + "! It may not exist."); 
			getline(fin, line);
			while ((fin) && (line[0] == '>')) {
	istringstream strin(line.substr(1));
	strin>>seqname;
	
	tseq = "";
	while((getline(fin, line)) && (line[0] != '>')) {
		tseq += line;
	}
	assert(tseq.length() > 0);
	name2seq[seqname] = tseq;
	transcript_id = seqname;
	gene_id = seqname;
	
	if (mappingType > 0) {
		mi_iter = mi_table.find(seqname);
		general_assert(mi_iter != mi_table.end(), "Mapping Info is not correct, cannot find " + seqname + "'s gene_id!");
		gene_id = mi_iter->second;
		if (mappingType == 2) {
			mi_iter2 = mi_table2.find(seqname);
			general_assert(mi_iter2 != mi_table2.end(), "Mapping Info is not correct, cannot find allele " + seqname + "'s transcript_id!");
			transcript_id = mi_iter2->second;
		}
	}
	
	vec.assign(1, Interval(1, (int)tseq.length()));
	transcripts.add(Transcript(transcript_id, gene_id, seqname, '+', vec, ""));
	++M;

	if (verbose && M % 1000000 == 0) { printf("%d sequences are processed!\n", M); }
			}
			fin.close();
		}

		assert(M == transcripts.getM());
		general_assert(M > 0, "The reference contains no transcripts!");
		transcripts.sort();
	
		// build refs
		for (int i = 1; i <= M; ++i) {
			seqname = transcripts.getTranscriptAt(i).getSeqName();
			n2s_iter = name2seq.find(seqname);
			general_assert(n2s_iter != name2seq.end(), "Cannot recognize sequence ID" + seqname + "!");
			refs.addRef(seqname, n2s_iter->second);
		}
	}

	writeToDisk(argv[1]); // write out generated indices
	
	return 0;  
}
