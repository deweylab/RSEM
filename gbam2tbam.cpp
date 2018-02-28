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

#include <sys/stat.h>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <fstream>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "utils.h"
#include "my_assert.h"
#include "Transcript.hpp"
#include "Transcripts.hpp"
#include "GenomeMap.hpp"
#include "SamParser.hpp"
#include "BamWriter.hpp"
#include "CIGARstring.hpp"
#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"
#include "ConversionGroup.hpp"

using namespace std;

bool verbose = true;



int num_threads;
htsThreadPool p = {NULL, 0};

int M;
char tiF[STRLEN], dupF[STRLEN], dupidF[STRLEN], genomeMapF[STRLEN], transListF[STRLEN];
Transcripts transcripts, dups;
vector<int> dup_ids;
GenomeMap genomeMap;

SamParser *parser;
BamWriter *writer;
AlignmentGroup ag;
ConversionGroup cg;


struct cigar_type {
	CIGARstring ci;
	uint32_t n_cigar, max_size, *cigar;
	int ref_len; // aligned length indicated by cigar

	cigar_type() {
		n_cigar = 0;
		max_size = STRLEN;
		cigar = (uint32_t*)calloc(max_size, sizeof(uint32_t));
		ref_len = 0;
	}

	~cigar_type() { free(cigar); }

	void set_ci(BamAlignment* ba, int mate) {
		ba->getCIGAR(ci, mate);
		ci.setCurrent();
	}

	// calculate cigar info for transcripts
	void calc_trans_cigar() {
		int len = ci.getLen();
		uint32_t value;
		int op, oplen;

		if (len > max_size) { max_size = len; kroundup32(max_size); cigar = (uint32_t*)realloc(cigar, sizeof(uint32_t) * max_size); }

		ci.setDir('+');
		n_cigar = 0; ref_len = 0;
		for (int i = 0; i < len; ++i) {
			value = ci.valueAt(i);
			op = bam_cigar_op(value);
			oplen = bam_cigar_oplen(value);

			if (op != BAM_CREF_SKIP) {
				if (n_cigar > 0 && bam_cigar_op(cigar[n_cigar - 1]) == op)
					cigar[n_cigar - 1] = bam_cigar_gen(bam_cigar_oplen(cigar[n_cigar - 1]) + oplen, op);
				else cigar[n_cigar++] = value;
				ref_len += (bam_cigar_type(op) & 2) ? oplen : 0;
			}
		}
	}
};


inline bool exists(const char* file) {
	struct stat buffer;
	return (stat(file, &buffer) == 0);
}

void load_references(const char* reference_name) {
	sprintf(tiF, "%s.ti", reference_name);
	transcripts.readFrom(tiF);
	M = transcripts.getM();
	transcripts.updateCLens();
	printf("Transcripts are loaded.\n");

	sprintf(dupF, "%s.dup.ti", reference_name);
	if (exists(dupF)) {
		dups.readFrom(dupF);
		dups.updateCLens();
		sprintf(dupidF, "%s.dup.ids", reference_name);

		int dup_id;
		ifstream fin(dupidF);

		assert(fin.is_open());
		dup_ids.clear();
		while (fin>> dup_id) dup_ids.push_back(dup_id);

		printf("Duplicated transcripts are loaded.\n");
	}

	sprintf(genomeMapF, "%s.g2t.map", reference_name);
	genomeMap.readFrom(genomeMapF);
	printf("Genome map is loaded.\n");
}

inline bool check_consistency(int pos, const Exon& exon, CIGARstring& ci) {
	const Transcript& transcript = (exon.tid <= M ? transcripts.getTranscriptAt(exon.tid) : dups.getTranscriptAt(exon.tid - M));
	const vector<Interval>& structures = transcript.getStructure();
	int s = structures.size(), spos = exon.eid;
	int cl = ci.getLen(), oplen;

	for (int cpos = 0; cpos < cl; ++cpos) {
		oplen = ci.oplenAt(cpos);

		if (ci.opAt(cpos) == BAM_CREF_SKIP) {
			if (pos - 1 != structures[spos].end) return false;
			pos += oplen; ++spos;
			if (spos >= s || pos != structures[spos].start) return false;
		}
		else {
			if (ci.optypeAt(cpos) & 2) {
				pos += oplen;
				if (pos - 1 > structures[spos].end) return false;				
			}
		}
	}

	return true;
}

void merge_candidates(vector<Exon>& list1, vector<Exon>& list2) {
	int s1 = list1.size(), s2 = list2.size();
	int p1 = 0, p2 = 0, np1 = 0, np2 = 0;

	while (p1 < s1 && p2 < s2) {
		if (list1[p1].tid < list2[p2].tid) ++p1;
		else if (list1[p1].tid > list2[p2].tid) ++p2;
		else { // equal
			if (np1 < p1) list1[np1] = list1[p1];
			if (np2 < p2) list2[np2] = list2[p2];
			++np1; ++np2; ++p1; ++p2;
		}
	}

	list1.resize(np1);
	list2.resize(np2);
}

// tid is 0-based here
void convertCoord(int gpos, const Exon& exon, int ref_len, char dir, int& tid, int& pos, bool& is_rev) {
	const Transcript& transcript = (exon.tid <= M ? transcripts.getTranscriptAt(exon.tid) : dups.getTranscriptAt(exon.tid - M));
	const vector<Interval>& structures = transcript.getStructure();

	tid = (exon.tid <= M ? exon.tid : dup_ids[exon.tid - M - 1]) - 1;
	pos = structures[exon.eid].clen + (gpos - structures[exon.eid].start); // 0-based
	if (transcript.getStrand() == '-') pos = transcript.getLength() - pos - ref_len;
	is_rev = transcript.getStrand() != dir;
}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		printf("Usage: rsem-gbam2tbam reference_name input.bam output.bam [-p num_threads]\n");
		exit(-1);
	}

	num_threads = 1;
	for (int i = 4; i < argc; ++i) {
		if (!strcmp(argv[i], "-p")) num_threads = atoi(argv[i + 1]);
	}

	load_references(argv[1]);

	if (num_threads > 1) assert(p.pool = hts_tpool_init(num_threads));
	parser = new SamParser(argv[2], num_threads > 1 ? &p : NULL);
	writer = new BamWriter(argv[3], parser->getHeader(), num_threads > 1 ? &p : NULL);
	sprintf(transListF, "%s.translist", argv[1]);
	writer->replaceReferences(transListF);
	writer->addProgram("rsem-gbam2tbam", VERSION, generateCommand(argc, argv));
	writer->writeHeader();

	int s, pos1, pos2;
	vector<Exon> list1, list2;
	cigar_type c1, c2;
	string iv_type1, iv_type2;
	BamAlignment *ba, *nba;

	int tid, tpos;
	bool is_rev;
	double frac;

	while (ag.read(parser)) {
		if (ag.isAligned()) {
			for (int i = 0; i < ag.size(); ++i) {
				ba = ag.getAlignment(i);
				if (ba->isAligned() & 1) {
					c1.set_ci(ba, 1);
					pos1 = ba->getPos(1) + 1;
					list1.clear();
					const vector<Exon>& exon_list = genomeMap.getExonList(parser->getSeqName(ba->get_tid()), pos1, iv_type1);
					for (int j = 0; j < (int)exon_list.size(); ++j) 
						if (check_consistency(pos1, exon_list[j], c1.ci)) list1.push_back(exon_list[j]);
				}
				if (ba->isAligned() & 2) {
					c2.set_ci(ba, 2);
					pos2 = ba->getPos(2) + 1;
					list2.clear();
					const vector<Exon>& exon_list = genomeMap.getExonList(parser->getSeqName(ba->get_tid()), pos2, iv_type2);
					for (int j = 0; j < (int)exon_list.size(); ++j) 
						if (check_consistency(pos2, exon_list[j], c2.ci)) list2.push_back(exon_list[j]);
				}
				if (ba->isAligned() == 3) merge_candidates(list1, list2);

				s = (ba->isAligned() & 1) ? list1.size() : list2.size();
				if (s > 0) {
					if (ba->isAligned() & 1) c1.calc_trans_cigar();
					if (ba->isAligned() & 2) c2.calc_trans_cigar();

					frac = 1.0; // ba->getFrac(); 
					if (frac > 0.0) frac /= s;

					for (int j = 0; j < s; ++j) {
						nba = cg.getNewBA(ba);
						if (ba->isAligned() & 1) {
							convertCoord(pos1, list1[j], c1.ref_len, ba->getMateDir(1), tid, tpos, is_rev);
							nba->setConvertedInfo(1, tid, tpos, is_rev, c1.n_cigar, c1.cigar);
						}
						if (ba->isAligned() & 2) {
							convertCoord(pos2, list2[j], c2.ref_len, ba->getMateDir(2), tid, tpos, is_rev);
							nba->setConvertedInfo(2, tid, tpos, is_rev, c2.n_cigar, c2.cigar);
						}
						cg.pushBackBA(ba, frac);
					}
				}
			}
			if (cg.size() == 0) cg.asUnmap(ba, iv_type1, iv_type2);
			cg.write(writer); // write and then clear the content for next alignment group
		}
		else ag.write(writer);
	}

	delete parser;
	delete writer;
	if (num_threads > 1) hts_tpool_destroy(p.pool);

	return 0;
}
