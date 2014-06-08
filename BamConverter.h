#ifndef BAMCONVERTER_H_
#define BAMCONVERTER_H_

#include<cstdio>
#include<cstring>
#include<cassert>
#include<string>
#include<map>

#include <stdint.h>
#include "sam/bam.h"
#include "sam/sam.h"
#include "sam_rsem_aux.h"
#include "sam_rsem_cvt.h"

#include "utils.h"
#include "my_assert.h"
#include "bc_aux.h"
#include "Transcript.h"
#include "Transcripts.h"

class BamConverter {
public:
	BamConverter(const char*, const char*, const char*, Transcripts&);
	~BamConverter();

	void process();
private:
	samfile_t *in, *out;
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

BamConverter::BamConverter(const char* inpF, const char* outF, const char* chr_list, Transcripts& transcripts)
	: transcripts(transcripts)
{
	general_assert(transcripts.getType() == 0, "Genome information is not provided! RSEM cannot convert the transcript bam file!");

	in = samopen(inpF, "rb", NULL);
	assert(in != 0);

	transcripts.buildMappings(in->header->n_targets, in->header->target_name);

	bam_header_t *out_header = sam_header_read2(chr_list);

	refmap.clear();
	for (int i = 0; i < out_header->n_targets; i++) {
		refmap[out_header->target_name[i]] = i;
	}

	if (in->header->l_text > 0) {
		char comment[] = "@CO\tThis BAM file is processed by rsem-tbam2gam to convert from transcript coordinates into genomic coordinates.\n";
		int comment_len = strlen(comment);

		//Filter @SQ fields if the BAM file is user provided
		char *text = in->header->text;
		int l_text = in->header->l_text;
		char *new_text = new char[l_text + comment_len];
		int pos = 0, s = 0;
		while (pos < l_text) {
			if ((pos + 2 < l_text) && (text[pos] == '@') && (text[pos + 1] == 'S') && (text[pos + 2] == 'Q')) {
				pos += 3;
				while (pos < l_text && text[pos] != '\n') ++pos;
			}
			else new_text[s++] = text[pos];
			++pos;
		}
		strncpy(new_text + s, comment, comment_len);
		s += comment_len;

		append_header_text(out_header, new_text, s);
		delete[] new_text;
	}

	out = samopen(outF, "wb", out_header);
	assert(out != 0);

	bam_header_destroy(out_header);
}

BamConverter::~BamConverter() {
	samclose(in);
	samclose(out);
}

void BamConverter::process() {
	bam1_t *b, *b2;
	std::string cqname;
	bool isPaired = false;

	HIT_INT_TYPE cnt = 0;

	cqname = "";
	b = bam_init1(); b2 = bam_init1();

	while (samread(in, b) >= 0) {
		++cnt;
		isPaired = (b->core.flag & 0x0001) > 0;
		if (isPaired) {
			assert(samread(in, b2) >= 0 && (b2->core.flag & 0x0001));
			assert((b->core.flag & 0x0001) && (b2->core.flag & 0x0001));
			assert(((b->core.flag & 0x0040) && (b2->core.flag & 0x0080)) || ((b->core.flag & 0x0080) && (b2->core.flag & 0x0040)));
 			++cnt;
		}

		if (cnt % 1000000 == 0) { printf("."); fflush(stdout); }

		// at least one segment is not properly mapped
		bool notgood = (b->core.flag & 0x0004) || (isPaired && (b2->core.flag & 0x0004));
		
		if (isPaired && notgood) general_assert((b->core.flag & 0x0004) && (b2->core.flag & 0x0004), cstrtos(bam1_qname(b)) + "'s two mates are not all marked as unalignable!");

		if (!notgood) {
			// for collapsing
			if (isPaired) {
				assert(b->core.tid == b2->core.tid);
				general_assert(b->core.tid == b2->core.tid, cstrtos(bam1_qname(b)) + "'s two mates are aligned to two different transcripts!");
				if ((b->core.flag & 0x0080) && (b2->core.flag & 0x0040)) {
					bam1_t *tmp = b; b = b2; b2 = tmp;
				}
			}

			const Transcript& transcript = transcripts.getTranscriptViaEid(b->core.tid + 1);

			convert(b, transcript);
			if (isPaired) {
				convert(b2, transcript);
				b->core.mpos = b2->core.pos;
				b2->core.mpos = b->core.pos;
			}

			if (cqname != bam1_qname(b)) {
				writeCollapsedLines();
				cqname = bam1_qname(b);
				collapseMap.init(isPaired);
			}

			uint8_t *p = bam_aux_get(b, "ZW");
			float prb = (p != NULL? bam_aux2f(p) : 1.0);
			collapseMap.insert(b, b2, prb);
		}
		else {
			assert(cqname != bam1_qname(b));

			writeCollapsedLines();
			cqname = bam1_qname(b);
			collapseMap.init(isPaired);

			samwrite(out, b);
			if (isPaired) samwrite(out, b2);
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
	if (b->core.flag & 0x0001) { b->core.mtid = b->core.tid; }
	b->core.qual = 255; // set to not available temporarily

	if (transcript.getStrand() == '-') {
		b->core.flag ^= 0x0010;
		if (b->core.flag & 0x0001) {
			b->core.flag ^= 0x0020;
			b->core.isize = -b->core.isize;
		}
		flipSeq(bam1_seq(b), readlen);
		flipQual(bam1_qual(b), readlen);
	}

	std::vector<uint32_t> data;
	data.clear();

	int core_pos, core_n_cigar;
	tr2chr(transcript, pos + 1, pos + readlen, core_pos, core_n_cigar, data);
	assert(core_pos >= 0);

	int rest_len = b->data_len - b->core.l_qname - b->core.n_cigar * 4;
	b->data_len = b->core.l_qname + core_n_cigar * 4 + rest_len;
	expand_data_size(b);
	uint8_t* pt = b->data + b->core.l_qname;
	memmove(pt + core_n_cigar * 4, pt + b->core.n_cigar * 4, rest_len);
	for (int i = 0; i < core_n_cigar; i++) { memmove(pt, &data[i], 4); pt += 4; }

	b->core.pos = core_pos;
	b->core.n_cigar = core_n_cigar;
	b->core.bin = bam_reg2bin(b->core.pos, bam_calend(&(b->core), bam1_cigar(b)));

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
				tmp_b->core.qual = getMAPQ(prb);
			}
			// otherwise, just use the MAPQ score of the orignal alignment

			samwrite(out, tmp_b);
			if (isPaired) {
				if (p != NULL) memcpy(bam_aux_get(tmp_b2, "ZW") + 1, (uint8_t*)&(prb), bam_aux_type2size('f'));
				tmp_b2->core.qual = tmp_b->core.qual;
				samwrite(out, tmp_b2);
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
	for (int i = 0; i < readlen; i++) {
		switch (bam1_seqi(s, readlen - i - 1)) {
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

	for (int i = 0; i < (int)seq.size(); i++) s[i] = seq[i];
}

inline void BamConverter::flipQual(uint8_t* q, int readlen) {
	int32_t mid = readlen / 2;
	uint8_t tmp;
	for (int i = 0; i < mid; i++) {
		tmp = q[i]; q[i] = q[readlen - i - 1]; q[readlen -i -1] = tmp;
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
	uint32_t* p = bam1_cigar(b);
	for (int i = 0; i < (int)b->core.n_cigar; i++)
		if ((*(p + i) & BAM_CIGAR_MASK) == BAM_CREF_SKIP) { hasN = true; break; }
	if (hasN) bam_aux_append(b, "XS", 'A', 1, (uint8_t*)&strand);
}

#endif /* BAMCONVERTER_H_ */
