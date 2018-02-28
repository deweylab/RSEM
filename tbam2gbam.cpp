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
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <fstream>
#include <utility>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "utils.h"
#include "my_assert.h"
#include "Transcript.hpp"
#include "Transcripts.hpp"
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
char tiF[STRLEN], chrListF[STRLEN];
Transcripts transcripts;

SamParser *parser;
BamWriter *writer;
AlignmentGroup ag;
ConversionGroup cg;


struct cigar_array {
	uint32_t n_cigar, max_size, *cigar;

	cigar_array() {
		n_cigar = 0;
		max_size = STRLEN;
		cigar = (uint32_t*)calloc(max_size, sizeof(uint32_t));
	}

	~cigar_type() { free(cigar); }

	void expand_size(uint32_t size) {
		if (size > max_size) {
			max_size = size;
			kroundup32(max_size);
			cigar = (uint32_t*)realloc(cigar, sizeof(uint32_t) * max_size);
		}
	}

	void flip() {
		for (uint32_t i = 0; i < (n_cigar >> 1); ++i) swap(cigar[i], cigar[n_cigar - i - 1]);
	}
};

// return the largest index r such as structures[r].clen <= pos
// intervals are [structures[i].clen, struxtures[i + 1].clen)
// pos is 0-based
int binary_search(const std::vector<Interval>& structures, int pos) {
	int l, r, mid;

	l = 0; r = structures.size() - 1;
	while (l <= r) {
		mid = (l + r) >> 1;
		if (pos >= structures[mid].clen) l = mid + 1;
		else r = mid - 1;
	}

	return r;
}

// pos : 0-based genomic coordinate
void convert_t2g(BamAlignment* ba, int mate, const Transcript& transcript, int& pos, bool& is_rev, cigar_array& ca) {
	int tpos, s, idx, endpos;
	const vector<Interval>& structures = transcript.getStructure();
	s = structures.size();

	int len;
	uint32_t value;
	int op, oplen, residue;
	CIGARstring ci;

	is_rev = transcript.getStrand() != ba->getMateDir(mate);

	tpos = ba->getPos(mate);
	idx = binary_search(structures, tpos);

	ba->getCIGAR(ci, mate);
	ci.setCurrent();
	len = ci.getLen();
	ca.expand_size(len + s);
	ca.n_cigar = 0;

	if (transcript.getStrand() == '+') {
		pos = structures[idx].start + (tpos - structures[idx].clen);
		endpos = pos - 1;
		for (int i = 0; i < len; ++i) {
			value = ci.valueAt(i);
			op = bam_cigar_op(value);
			oplen = bam_cigar_oplen(value);
			if (bam_cigar_type(op) & 2) {
				if (idx < s - 1 && endpos == structures[idx].end) endpos = structures[++idx].start - 1;
				endpos += oplen;
				while (idx < s - 1 && endpos > structures[idx].end) {
					residue = endpos - structures[idx].end;
					ca.cigar[n_cigar++] = bam_cigar_gen(oplen - residue, op);
					ca.cigar[n_cigar++] = bam_cigar_gen((structures[idx + 1].start - 1) - structures[idx].end, BAM_CREF_SKIP);
					oplen = residue; 
					endpos = structures[++idx].start + oplen - 1;
				}
				assert(endpos <= structures[idx].end);
				ca.cigar[n_cigar++] = bam_cigar_gen(oplen, op);
			}
			else ca.cigar[n_cigar++] = value;
		}
	}
	else {
		endpos = structures[idx].start = (tpos - structures[idx].clen);
		pos = endpos + 1;
		for (int i = 0; i < len; ++i) {
			value = ci.valueAt(i);
			op = bam_cigar_op(value);
			oplen = bam_cigar_oplen(value);
			if (bam_cigar_type(op) & 2) {
				if (idx > 0 && pos == structures[idx].start) pos = structures[--idx].end + 1;
				pos -= oplen;
				while (idx > 0 && pos < structures[idx].start) {
					residue = structures[idx].start - pos;
					ca.cigar[n_cigar++] = bam_cigar_gen(oplen - residue, op);
					ca.cigar[n_cigar++] = bam_cigar_gen((structures[idx].start - 1) - structures[idx - 1].end, BAM_CREF_SKIP);
					oplen = residue;
					pos = structures[--idx].end - oplen + 1;
				}
				assert(pos >= structures[idx].start);
				ca.cigar[n_cigar++] = bam_cigar_gen(oplen, op);		
			}
			else ca.cigar[n_cigar++] = value;
		}
	}

	--pos;
	if (ba->getMateDir(mate) == '-') ca.flip();
}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		printf("Usage: rsem-tbam2tbam reference_name input.bam output.bam [-p num_threads]\n");
		exit(-1);
	}

	num_threads = 1;
	for (int i = 4; i < argc; ++i) {
		if (!strcmp(argv[i], "-p")) num_threads = atoi(argv[i + 1]);
	}
	if (num_threads > 1) assert(p.pool = hts_tpool_init(num_threads));

	sprintf(tiF, "%s.ti", reference_name);
	transcripts.readFrom(tiF);
	M = transcripts.getM();
	transcripts.updateCLens();
	printf("Transcripts are loaded.\n");

	parser = new SamParser(argv[2], num_threads > 1 ? &p : NULL);
	writer = new BamWriter(argv[3], parser->getHeader(), num_threads > 1 ? &p : NULL);
	sprintf(chrListF, "%s.chrlist", argv[1]);
	writer->replaceReferences(chrListF);
	writer->addProgram("rsem-tbam2gbam", VERSION, generateCommand(argc, argv));
	writer->writeHeader();

	BamAlignment *ba, *nba;
	int cid, cpos; // cid, chromosome id; cpos, chromosome position
	bool is_rev;
	cigar_array ca;
	double frac;

	while (ag.read(parser)) {
		if (ag.isAligned()) {
			for (int i = 0; i < ag.size(); ++i) {
				ba = ag.getAlignment(i);
				const Transcript& transcript = transcripts.getTranscriptAt(ba->getTid());
				cid = writer->name2id(transcript.getSeqName().c_str());

				frac = ba->getFrac();

				nba = cg.getNewBA(ba);
				if (ba->isAligned() & 1) {
					convert_t2g(ba, 1, transcript, cpos, is_rev, ca);
					nba->setConvertedInfo(1, cid, cpos, is_rev, ca.n_cigar, ca.cigar);
				}
				if (ba->isAligned() & 2) {
					convert_t2g(ba, 2, transcript, cpos, is_rev, ca);
					nba->setConvertedInfo(2, cid, cpos, is_rev, ca.n_cigar, ca.cigar);
				}
				cg.pushBackBA(ba, frac, transcript.getStrand());
			}
			assert(cg.size() > 0);
			cg.write(writer); 
			cg.wrapUp(); // wrap up for this alignment and prepare for the next one
		}
		else ag.write(writer);
	}

	delete parser;
	delete writer;
	if (num_threads > 1) hts_tpool_destroy(p.pool);

	return 0;
}
