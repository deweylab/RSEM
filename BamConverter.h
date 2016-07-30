#ifndef BAMCONVERTER_H_
#define BAMCONVERTER_H_

#include<cstdio>
#include<cstring>
#include<cassert>
#include<string>
#include<map>

#include <stdint.h>
#include "htslib/sam.h"
#include "sam_utils.h"
#include "SamHeader.hpp"

#include "utils.h"
#include "my_assert.h"
#include "bc_aux.h"
#include "Transcript.h"
#include "Transcripts.h"

class BamConverter {
public:
  BamConverter(const char* inpF, const char* outF, const char* chr_list, Transcripts& transcripts, int nThreads, const std::string& command);
	~BamConverter();

	void process();
private:
	samFile *in, *out;
	bam_hdr_t *in_header, *out_header;
	Transcripts& transcripts;

	std::map<std::string, int> refmap;
	std::map<std::string, int>::iterator iter;

	CollapseMap collapseMap;

	void convert(bam1_t*, const Transcript&);

	void writeCollapsedLines();
	void flipSeq(uint8_t*, int);
	void flipQual(uint8_t*, int);
	void modifyTags(bam1_t*, const Transcript&); // modify MD tag and XS tag if needed
};

BamConverter::BamConverter(const char* inpF, const char* outF, const char* chr_list, Transcripts& transcripts, int nThreads, const std::string& command)
	: transcripts(transcripts)
{
	general_assert(transcripts.getType() == 0, "Genome information is not provided! RSEM cannot convert the transcript bam file!");

	in = sam_open(inpF, "r");
	assert(in != 0);
	in_header = sam_hdr_read(in);
	assert(in_header != 0);

	transcripts.buildMappings(in_header->n_targets, in_header->target_name);

	SamHeader hdr(in_header->text);
	hdr.replaceSQ(chr_list);
	hdr.insertPG("rsem-tbam2gbam", command);
	//	hdr.addComment("This BAM file is processed by rsem-tbam2gam to convert from transcript coordinates into genomic coordinates.");
	out_header = hdr.create_header();
	
	refmap.clear();
	for (int i = 0; i < out_header->n_targets; ++i) {
		refmap[out_header->target_name[i]] = i;
	}

	out = sam_open(outF, "wb");
	assert(out != 0);
	sam_hdr_write(out, out_header);

	if (nThreads > 1) general_assert(hts_set_threads(out, nThreads) == 0, "Fail to create threads for writing the BAM file!");
}

BamConverter::~BamConverter() {
	bam_hdr_destroy(in_header);
	sam_close(in);
	bam_hdr_destroy(out_header);
	sam_close(out);
}

void BamConverter::process() {
	bam1_t *b, *b2;
	std::string cqname, qname;
	bool isPaired = false;

	HIT_INT_TYPE cnt = 0;

	cqname = "";
	b = bam_init1(); b2 = bam_init1();

	while (sam_read1(in, in_header, b) >= 0) {
		++cnt;
		isPaired = bam_is_paired(b);
		if (isPaired) {
		  assert(sam_read1(in, in_header, b2) >= 0 && bam_is_paired(b2));
		  if (!bam_is_read1(b)) { bam1_t *tmp = b; b = b2; b2 = tmp; }
		  assert(bam_is_read1(b) && bam_is_read2(b2));
		  general_assert((bam_is_mapped(b) && bam_is_mapped(b2)) || (bam_is_unmapped(b) && bam_is_unmapped(b2)), \
				 "Detected partial alignments for read " + bam_get_canonical_name(b) + ", which RSEM currently does not support!"); 
		  ++cnt;
		}

		if (cnt % 1000000 == 0) { printf("."); fflush(stdout); }

		qname = bam_get_canonical_name(b);
		if (bam_is_mapped(b)) {
		  // for collapsing
		  if (isPaired) general_assert(b->core.tid == b2->core.tid, qname + "'s two mates are aligned to two different transcripts!");

		  const Transcript& transcript = transcripts.getTranscriptViaEid(b->core.tid + 1);

		  convert(b, transcript);
		  if (isPaired) {
		    convert(b2, transcript);
		    b->core.mpos = b2->core.pos;
		    b2->core.mpos = b->core.pos;
		  }

		  if (cqname != qname) {
		    writeCollapsedLines();
		    cqname = qname;
		    collapseMap.init(isPaired);
		  }
		  
		  uint8_t *p = bam_aux_get(b, "ZW");
		  float prb = (p != NULL? bam_aux2f(p) : 1.0);
		  collapseMap.insert(b, b2, prb);
		}
		else {
		  assert(cqname != qname);
		  
		  writeCollapsedLines();
		  cqname = qname;
		  collapseMap.init(isPaired);
		  
		  sam_write1(out, out_header, b);
		  if (isPaired) sam_write1(out, out_header, b2);
		}
	}

	writeCollapsedLines();

	bam_destroy1(b);
	bam_destroy1(b2);

	if (cnt >= 1000000) printf("\n");
}

void BamConverter::convert(bam1_t* b, const Transcript& transcript) {
	int pos = b->core.pos;
	int readlen = b->core.l_qseq;

	general_assert(readlen > 0, "One alignment line has SEQ field as *. RSEM does not support this currently!");

	iter = refmap.find(transcript.getSeqName());
	assert(iter != refmap.end());
	b->core.tid = iter->second;
	if (bam_is_paired(b)) { b->core.mtid = b->core.tid; }
	b->core.qual = 255; // set to not available temporarily

	if (transcript.getStrand() == '-') {
		b->core.flag ^= BAM_FREVERSE;
		if (bam_is_paired(b)) {
			b->core.flag ^= BAM_FMREVERSE;
			b->core.isize = -b->core.isize;
		}
		flipSeq(bam_get_seq(b), readlen);
		flipQual(bam_get_qual(b), readlen);
	}

	std::vector<uint32_t> data;
	data.clear();

	int core_pos, core_n_cigar;
	tr2chr(transcript, pos + 1, pos + readlen, core_pos, core_n_cigar, data);
	assert(core_pos >= 0);

	int rest_len = b->l_data - b->core.l_qname - b->core.n_cigar * 4;
	b->l_data = b->core.l_qname + core_n_cigar * 4 + rest_len;
	expand_data_size(b);
	uint8_t* pt = b->data + b->core.l_qname;
	memmove(pt + core_n_cigar * 4, pt + b->core.n_cigar * 4, rest_len);
	for (int i = 0; i < core_n_cigar; ++i) { memmove(pt, &data[i], 4); pt += 4; }

	b->core.pos = core_pos;
	b->core.n_cigar = core_n_cigar;
	b->core.bin = hts_reg2bin(b->core.pos, bam_endpos(b), 14, 5);

	modifyTags(b, transcript); // check if need to add XS tag, if need, add it
}

inline void BamConverter::writeCollapsedLines() {
	bam1_t *tmp_b = NULL,*tmp_b2 = NULL;
	float prb;
	bool isPaired;
	uint8_t *p;

	if (!collapseMap.empty(isPaired)) {
		while (collapseMap.next(tmp_b, tmp_b2, prb)) {
			p = bam_aux_get(tmp_b, "ZW");
			if (p != NULL) {
				memcpy(bam_aux_get(tmp_b, "ZW") + 1, (uint8_t*)&(prb), bam_aux_type2size('f'));
				tmp_b->core.qual = bam_prb_to_mapq(prb);
			}
			// otherwise, just use the MAPQ score of the orignal alignment

			sam_write1(out, out_header, tmp_b);
			if (isPaired) {
				if (p != NULL) memcpy(bam_aux_get(tmp_b2, "ZW") + 1, (uint8_t*)&(prb), bam_aux_type2size('f'));
				tmp_b2->core.qual = tmp_b->core.qual;
				sam_write1(out, out_header, tmp_b2);
			}
			bam_destroy1(tmp_b);
			if (isPaired) bam_destroy1(tmp_b2);
		}
	}
}

inline void BamConverter::flipSeq(uint8_t* s, int readlen) {
	uint8_t code, base;
	std::vector<uint8_t> seq;

	code = 0; base = 0;
	seq.clear();
	for (int i = 0; i < readlen; ++i) {
		switch (bam_seqi(s, readlen - i - 1)) {
		case 1: base = 8; break;
		case 2: base = 4; break;
		case 4: base = 2; break;
		case 8: base = 1; break;
		case 15: base = 15; break;
		default: assert(false);
		}
		code |=  base << (4 * (1 - i % 2));
		if (i % 2 == 1) { seq.push_back(code); code = 0; }
	}
	if (readlen % 2 == 1) { seq.push_back(code); }

	for (int i = 0; i < (int)seq.size(); ++i) s[i] = seq[i];
}

inline void BamConverter::flipQual(uint8_t* q, int readlen) {
	int32_t mid = readlen / 2;
	uint8_t tmp;
	for (int i = 0; i < mid; ++i) {
		tmp = q[i]; q[i] = q[readlen - i - 1]; q[readlen - i - 1] = tmp;
	}
}

inline void BamConverter::modifyTags(bam1_t* b, const Transcript& transcript) {
	char strand = transcript.getStrand();
	uint8_t *s = NULL;

	if (strand == '-') {
	  s = bam_aux_get(b, "MD");
	  if ((s != NULL) && (*(s) == 'Z') && (bam_aux2Z(s) != NULL)) {
	    char *mis = bam_aux2Z(s);
	    int len = strlen(mis);
	    char *tmp = new char[len];
	    int cur_type = -1, fr = -1, type, base;
	    for (int i = 0; i < len; i++) {
	      type = (mis[i] >= '0' && mis[i] <= '9');
	      if (cur_type != type) {
		switch(cur_type) {
		case 0:
		  base = len - 1;
		  if (mis[fr] == '^') { tmp[len - i] = mis[fr]; ++fr; ++base; }
		  for (int j = fr; j < i; j++) tmp[base - j] = ((mis[j] == 'A' || mis[j] == 'C' || mis[j] == 'G' || mis[j] == 'T') ? getOpp(mis[j]) : mis[j]);
		  break;
		case 1: 
		  base = len - i - fr;
		  for (int j = fr; j < i; j++) tmp[base + j] = mis[j]; 
		  break; 
		}
		cur_type = type;
		fr = i;
	      }
	    }
	    switch(cur_type) {
	    case 0:
	      base = len - 1;
	      if (mis[fr] == '^') { tmp[0] = mis[fr]; ++fr; ++base; }
	      for (int j = fr; j < len; j++) tmp[base - j] = ((mis[j] == 'A' || mis[j] == 'C' || mis[j] == 'G' || mis[j] == 'T') ? getOpp(mis[j]) : mis[j]);
	      break;
	    case 1: 
	      for (int j = fr; j < len; j++) tmp[j - fr] = mis[j]; 
	      break; 
	    }
	    strncpy(mis, tmp, len);
	    delete[] tmp;
	  }
	}

	// append XS:A field if necessary
	s = bam_aux_get(b, "XS");
	if (s != NULL) bam_aux_del(b, s);
	bool hasN = false;
	uint32_t* p = bam_get_cigar(b);
	for (int i = 0; i < (int)b->core.n_cigar; i++)
		if ((*(p + i) & BAM_CIGAR_MASK) == BAM_CREF_SKIP) { hasN = true; break; }
	if (hasN) bam_aux_append(b, "XS", 'A', 1, (uint8_t*)&strand);
}

#endif /* BAMCONVERTER_H_ */
