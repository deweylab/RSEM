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
#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "utils.h"
#include "my_assert.h"
#include "GTFItem.h"
#include "Transcript.hpp"
#include "Transcripts.hpp"
#include "RefSeq.hpp"
#include "Refs.hpp"

using namespace std;

bool verbose = true; // define verbose


int M;
Transcripts transcripts;
Refs refs;

// hasGTF: if has GTF file ; hasMapping: if has mapping file; isAlleleSpecific: if allele-specific; rmdup : if remove duplicate sequences; appendPolyA : if append polyA tails; n2g_idx : if generate n2g index
bool hasGTF, hasMapping, isAlleleSpecific, rmdup, appendPolyA, n2g_idx; 

int num_files; // Number of reference input files
int mappingPos, file_pos; // position in argv vector for mapping file and reference files
int polyALen; // length of poly(A)s appended

char gtfF[STRLEN], polyAExcludeF[STRLEN];
char tiF[STRLEN], dupF[STRLEN], refFastaF[STRLEN];
char transListF[STRLEN], chromListF[STRLEN];
char groupF[STRLEN], gtF[STRLEN], taF[STRLEN]; 

// Mapping file related data structures
map<string, string> mi_table; // mapping info table
map<string, string>::iterator mi_iter; //mapping info table's iterator
map<string, string> mi_table2; // allele_id to transcript_id
map<string, string>::iterator mi_iter2; // corresponding iterator

set<string> sources, no_polyA_set; // sources, trusted sources; no_polyA_set, transcripts without polyA tails

vector<GTFItem> items;
map<string, vector<int> > sn2tr; // map from chromosome name to transcripts
map<string, vector<int> >::iterator sn2tr_iter;
vector<Interval> vec;

// chromosome info
struct ChrInfo {
	string name;
	size_t len;

	ChrInfo(const string& name, size_t len) : name(name), len(len) {}

	bool operator< (const ChrInfo& o) const {
		return name < o.name;
	}
};

vector<ChrInfo> chrvec;

struct CursorPos {
	char *filename;
	int line_no, pos;
};

CursorPos cursor;


int nDup;
map<string, int> hasSeen; // map from sequence to tid
pair<map<string, int>::iterator, bool> dup_ret;
Transcripts dups; // duplicated transcripts
vector<int> dup_tids; // tid of the original transcripts


// test if we should skip this line
inline bool skip(const string& line) {
	size_t pos = 0, len = line.length();
	while (pos < len && isspace(line[pos])) ++pos;
	return pos >= len || line[pos] == '#'; // skipping if empty line or commented line
}

// check nucleotides
inline char check(char c) {
	general_assert(isalpha(c), "FASTA file " + cstrtos(cursor.filename) + " contains an unknown character, " + \
					ctos(c) + " (ASCII code " + itos(c) + "), at line " + itos(cursor.line_no) + ", position " + itos(cursor.pos + 1) + "!");
	c = toupper(c);
	if (c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'N';

	return c;
}

// GTF related, trust sources
void parseSources(char* sstr) {
	char* p = strtok(sstr, ",");
	while (p != NULL) {
		sources.insert(p);
		p = strtok(NULL, ",");
	}
}

inline bool isTrusted(const string& source) {
	return sources.size() == 0 || sources.find(source) != sources.end();
}

// load polyA exclusion ids
void load_polyA_exclusion_set(char* polyAExcludeF) {
	ifstream fin(polyAExcludeF);
	string line;

	general_assert(fin.is_open(), "Cannot open " + cstrtos(polyAExcludeF) + "! It may not exist.");

	while (getline(fin, line)) {
		no_polyA_set.insert(line);
	}

	if (verbose) { printf("%s is loaded.\n", polyAExcludeF); }
}



// load mapping info
void loadMappingInfo(char* mappingF, bool isAlleleSpecific) {
	ifstream fin(mappingF);
	string line, key, value, value2;

	general_assert(fin.is_open(), "Cannot open " + cstrtos(mappingF) + "! It may not exist.");

	if (!isAlleleSpecific) {
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

	transcripts.add(Transcript(transcript_id, gene_id, seqname, strand, vec, left, transcript_name, gene_name));
}

void parse_gtf_file(char* gtfF) {
	ifstream fin(gtfF);
	string line, tid, gid;
	GTFItem item;

	general_assert(fin.is_open(), "Cannot open " + cstrtos(gtfF) + "! It may not exist.");

	int cnt = 0;
	int n_warns = 0;

	items.clear();
	while (getline(fin, line)) {
		if (skip(line)) continue;
		item.parse(line);
  		if (item.getFeature() == "exon" && isTrusted(item.getSource())) {
 			if (item.getStart() > item.getEnd()) {
				if (++n_warns <= MAX_WARNS) {
					fprintf(stderr, "Warning: exon's start position is larger than its end position! This exon is discarded.\n");
					fprintf(stderr, "\t%s\n\n", line.c_str());
				}
 			}
 			else if (item.getStart() < 1) {
				if (++n_warns <= MAX_WARNS) {
					fprintf(stderr, "Warning: exon's start position is less than 1! This exon is discarded.\n");
					fprintf(stderr, "\t%s\n\n", line.c_str());
				}
 			}
 			else {
 				item.parseAttributes(line);
 		 		if (hasMapping) {
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

	if (n_warns > 0) fprintf(stderr, "Warning: In total, %d exons are discarded.", n_warns);

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
	
	M = transcripts.getM();
	general_assert(M > 0, "The reference contains no transcripts!");

	if (verbose) { printf("Parsing GTF File is done!\n"); }
}

void extract_reference_sequences(char* argv[]) {
	ifstream fin;
	string line, gseq, tseq, seqname;
	int len;
	size_t seqlen;

	transcripts.setType(0);
	general_assert(!isAlleleSpecific, "RSEM could not extract allele-specific transcript sequences from a genome!");
	parse_gtf_file(gtfF);

	chrvec.clear();
	refs.setM(M);

	for (int i = 0; i < num_files; ++i, ++file_pos) {
		fin.open(argv[file_pos]);
		general_assert(fin.is_open(), "Cannot open " + cstrtos(argv[file_pos]) + "! It may not exist.");
		cursor.filename = argv[file_pos]; cursor.line_no = cursor.pos = 0;

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

			chrvec.push_back(ChrInfo(seqname, seqlen));
			
			sn2tr_iter = sn2tr.find(seqname);
			if (sn2tr_iter == sn2tr.end()) continue;
		
			vector<int>& vec = sn2tr_iter->second;
			int s = vec.size();
			for (int j = 0; j < s; ++j) {
				assert(vec[j] > 0 && vec[j] <= M);
				const Transcript& transcript = transcripts.getTranscriptAt(vec[j]);
				transcript.extractSeq(gseq, tseq);
				refs.setRef(vec[j], transcript.getTranscriptID(), tseq);
			}
		}
		fin.close();

		if (verbose) { printf("%s is processed!\n", argv[file_pos]); } 
	}
		
	sort(chrvec.begin(), chrvec.end());
	if (verbose) { printf("Extracting sequences is done!\n"); }
}

void load_reference_sequences(char* argv[]) {
	ifstream fin;
	string line, tseq;
	string seqname, gene_id, transcript_id;
	int len;
	size_t seqlen;

	map<string, string> name2seq;
	map<string, string>::iterator n2s_iter;


	transcripts.setType(isAlleleSpecific ? 2 : 1);
		
	M = 0;
	name2seq.clear();
	for (int i = 0; i < num_files; ++i, ++file_pos) {
		fin.open(argv[file_pos]);
		general_assert(fin.is_open(), "Cannot open " + cstrtos(argv[file_pos]) + "! It may not exist."); 
		cursor.filename = argv[file_pos]; cursor.line_no = cursor.pos = 0;

		getline(fin, line);
		while ((fin) && (line[0] == '>')) {
			istringstream strin(line.substr(1));
			strin>> seqname;
			++cursor.line_no;

			tseq = ""; seqlen = 0;
			while((getline(fin, line)) && (line[0] != '>')) {
				++cursor.line_no;
				len = line.length();
				for (cursor.pos = 0; cursor.pos < len; ++cursor.pos) line[cursor.pos] = check(line[cursor.pos]);
				seqlen += len;
				tseq += line;
			}
			assert(seqlen > 0);
			name2seq[seqname] = tseq;
			transcript_id = seqname;
			gene_id = seqname;
	
			if (hasMapping) {
				mi_iter = mi_table.find(seqname);
				general_assert(mi_iter != mi_table.end(), "Mapping Info is not correct, cannot find " + seqname + "'s gene_id!");
				gene_id = mi_iter->second;
				if (isAlleleSpecific) {
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

	refs.setM(M);
	for (int i = 1; i <= M; ++i) {
		seqname = transcripts.getTranscriptAt(i).getSeqName();		
		n2s_iter = name2seq.find(seqname);
		assert(n2s_iter != name2seq.end());
		refs.setRef(i, seqname, n2s_iter->second);
	}
}

void shrink() {
	int curp = 0;
	int n_warns = 0;

	nDup = 0; 
	if (rmdup) hasSeen.clear(), dup_tids.clear();
	for (int i = 1; i <= M; ++i) {
		const RefSeq& ref = refs.getRef(i);
		if (ref.getSeq() == "") {
			if (++n_warns <= MAX_WARNS) {
				const Transcript& transcript = transcripts.getTranscriptAt(i);
				fprintf(stderr, "Warning: Cannot extract transcript %s's sequence since the chromosome it locates, %s, is absent!\n", transcript.getTranscriptID().c_str(), transcript.getSeqName().c_str());
			}
		}
		else {
			++curp;
			if (rmdup) {
				dup_ret = hasSeen.insert(pair<string, int>(ref.getSeq(), curp));
				if (!dup_ret.second) {
					++nDup; --curp;
					dups.add(transcripts.getTranscriptAt(i));
					dup_tids.push_back(dup_ret.first->second);
					continue;
				}
			}
			transcripts.move(i, curp);
			refs.move(i, curp);
		} 
	}

	transcripts.setM(curp);
	refs.setM(curp);
	M = transcripts.getM();
	general_assert(M > 0, "The reference contains no transcripts!");

	if (n_warns > 0) fprintf(stderr, "Warning: %d transcripts are failed to extract because their chromosome sequences are absent.\n", n_warns);
	if (nDup > 0) fprintf(stderr, "%d duplicated transcripts are removed.\n", nDup);
	if (verbose) printf("%d transcripts are extracted.\n", curp);
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
	if (verbose) printf("Transcript information file is generated.\n");

	if (rmdup && nDup > 0) {
		sprintf(tiF, "%s.dup.ti", refName);
		dups.writeTo(tiF);

		sprintf(dupF, "%s.dup.ids", refName);
		fout.open(dupF);
		for (int i = 0; i < (int)dup_tids.size(); ++i) fout<< dup_tids[i]<< endl;
		fout.close();

		if (verbose) printf("Duplicated transcript information files is generated.\n");
	}
	


	sprintf(refFastaF, "%s.transcripts.fa", refName);
	refs.writeTo(refFastaF);

	if (appendPolyA) {
		for (int i = 1; i <= M; ++i)
			if (no_polyA_set.size() == 0 || no_polyA_set.find(refs.getRef(i).getName()) == no_polyA_set.end()) 
				refs.appendPolyATail(i, polyALen);
	}

	sprintf(refFastaF, "%s.idx.fa", refName);
	refs.writeTo(refFastaF);

	sprintf(transListF, "%s.translist", refName);
	refs.writeTransListTo(transListF);

	if (n2g_idx) {
		sprintf(refFastaF, "%s.n2g.idx.fa", refName);
		fout.open(refFastaF);
		for (int i = 1; i <= M; ++i) {
			const RefSeq& ref = refs.getRef(i);
			fout<< ">"<< ref.getName()<< endl<< n2g(ref.getSeq())<< endl;
		}
		fout.close();
		if (verbose) printf("%s is generated.\n", refFastaF);
	}



	if (hasGTF) {
		sprintf(chromListF, "%s.chrlist", refName);
		fout.open(chromListF);
		for (int i = 0; i < (int)chrvec.size(); ++i)
			fout<< chrvec[i].name<< '\t'<< chrvec[i].len<< endl;
		fout.close();
		if (verbose) printf("Chromosome list file is generated.\n");
	}


	
	string cur_gene_id, cur_transcript_id, name;
	vector<int> gi, gt, ta;

	cur_gene_id = ""; gi.clear(); 
	if (isAlleleSpecific) { cur_transcript_id = ""; gt.clear(); ta.clear(); }
	for (int i = 1; i <= M; ++i) {
		const Transcript& transcript = transcripts.getTranscriptAt(i);
		if (cur_gene_id != transcript.getGeneID()) {
			gi.push_back(i);
			if (isAlleleSpecific) gt.push_back((int)ta.size());
			cur_gene_id = transcript.getGeneID();
		}
		if (isAlleleSpecific && (cur_transcript_id != transcript.getTranscriptID())) {
			ta.push_back(i);
			cur_transcript_id = transcript.getTranscriptID();
		}
	}
	
	gi.push_back(M + 1);
	if (isAlleleSpecific) { gt.push_back((int)ta.size()); ta.push_back(M + 1); }

	sprintf(groupF, "%s.grp", refName);
	fout.open(groupF);
	for (int i = 0; i < (int)gi.size(); ++i) fout<< gi[i]<< endl;
	fout.close();
	if (verbose) printf("Group file is generated.\n");

	if (isAlleleSpecific) {
		sprintf(gtF, "%s.gt", refName);
		fout.open(gtF);
		for (int i = 0; i < (int)gt.size(); ++i) fout<< gt[i]<< endl;
		fout.close();
		sprintf(taF, "%s.ta", refName);
		fout.open(taF);
		for (int i = 0; i < (int)ta.size(); ++i) fout<< ta[i]<< endl;
		fout.close();
		if (verbose) printf("Allele-specific group files are generated.\n");
	}
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		printf("Usage: rsem-build-reference refName [--gtf gtfF] [--trusted-sources sources] [--mapping mappingF] [--allele-specific] [--remove-duplicates] [--polyA-length length] [--no-polyA-subset polyAExcludeF][--n2g-index] [-q] [--files num_of_files file_1 file_2 ...]\n");
		exit(-1);
	}

	hasGTF = false;
	hasMapping = false;
	isAlleleSpecific = false;
	rmdup = false;
	appendPolyA = false;
	n2g_idx = false;

	sources.clear();
	no_polyA_set.clear();

	int argpos = 2;
	while (argpos < argc) {
		if (!strcmp(argv[argpos], "--gtf")) {
			hasGTF = true;
			strcpy(gtfF, argv[++argpos]);
		}
		if (!strcmp(argv[argpos], "--trusted-sources")) parseSources(argv[++argpos]);
		if (!strcmp(argv[argpos], "--mapping")) {
			hasMapping = true;
			mappingPos = ++argpos;
		}
		if (!strcmp(argv[argpos], "--allele-specific")) isAlleleSpecific = true;
		if (!strcmp(argv[argpos], "--remove-duplicates")) rmdup = true;
		if (!strcmp(argv[argpos], "--polyA-length")) {
			appendPolyA = true;
			polyALen = atoi(argv[++argpos]);
		}
		if (!strcmp(argv[argpos], "--no-polyA-subset")) load_polyA_exclusion_set(argv[++argpos]);
		if (!strcmp(argv[argpos], "--n2g-index")) n2g_idx = true;
		if (!strcmp(argv[argpos], "-q")) verbose = false;
		if (!strcmp(argv[argpos], "--files")) {
			num_files = atoi(argv[++argpos]);
			file_pos = argpos + 1; // the position in argv for the first file
			argpos += num_files;
		}
		++argpos;
	}

	if (hasMapping) loadMappingInfo(argv[mappingPos], isAlleleSpecific);

	if (hasGTF) extract_reference_sequences(argv);
	else load_reference_sequences(argv);

	shrink(); // remove duplicates 

	writeToDisk(argv[1]); // write out generated indices
	
	return 0;  
}
