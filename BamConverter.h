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
	void addXSTag(bam1_t*, const Transcript&);
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

	append_header_text(out_header, in->header->text, in->header->l_text);

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
			assert(samread(in, b2) >= 0 && (b2->core.flag & 0x0001) && b->core.tid == b2->core.tid);
			assert((b->core.flag & 0x0040) && (b2->core.flag & 0x0080)); // for collapsing
			++cnt;
		}

		if (cnt % 1000000 == 0) { printf("."); fflush(stdout); }

		// at least one segment is not properly mapped
		if ((b->core.flag & 0x0004) || (isPaired && (b2->core.flag & 0x0004))) continue;

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

		collapseMap.insert(b, b2, bam_aux2f(bam_aux_get(b, "ZW")));
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

	addXSTag(b, transcript); // check if need to add XS tag, if need, add it
}

inline void BamConverter::writeCollapsedLines() {
	bam1_t *tmp_b = NULL,*tmp_b2 = NULL;
	float prb;
	bool isPaired;

	if (!collapseMap.empty(isPaired)) {
		while (collapseMap.next(tmp_b, tmp_b2, prb)) {
			memcpy(bam_aux_get(tmp_b, "ZW") + 1, (uint8_t*)&(prb), bam_aux_type2size('f'));
			tmp_b->core.qual = getMAPQ(prb);
			if (tmp_b->core.qual > 0) {
				samwrite(out, tmp_b);
				if (isPaired) {
					memcpy(bam_aux_get(tmp_b2, "ZW") + 1, (uint8_t*)&(prb), bam_aux_type2size('f'));
					tmp_b2->core.qual = tmp_b->core.qual;
					samwrite(out, tmp_b2);
				}
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

inline void BamConverter::addXSTag(bam1_t* b, const Transcript& transcript) {
	uint32_t* p = bam1_cigar(b);
	bool hasN = false;
	for (int i = 0; i < (int)b->core.n_cigar; i++)
		if ((*(p + i) & BAM_CIGAR_MASK) == BAM_CREF_SKIP) { hasN = true; break; }
	if (!hasN) return;
	char strand = transcript.getStrand();
	bam_aux_append(b, "XS", 'A', 1, (uint8_t*)&strand);
}

#endif /* BAMCONVERTER_H_ */
