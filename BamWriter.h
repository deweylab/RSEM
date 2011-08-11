#ifndef BAMWRITER_H_
#define BAMWRITER_H_

#include<cmath>
#include<cstdio>
#include<cstring>
#include<cassert>
#include<string>
#include<map>
#include<sstream>

#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"
#include "SingleHit.h"
#include "PairedEndHit.h"

#include "HitWrapper.h"
#include "Transcript.h"
#include "Transcripts.h"

class BamWriter {
public:
	BamWriter(char, const char*, const char*, const char*, const char*);
	~BamWriter();

	void work(HitWrapper<SingleHit>, Transcripts&);
	void work(HitWrapper<PairedEndHit>, Transcripts&);
private:
	samfile_t *in, *out;

	std::map<std::string, int> refmap;
	std::map<std::string, int>::iterator iter;

	struct SingleEndT {
		bam1_t *b;

		SingleEndT(bam1_t *b = NULL) {
			this->b = b;
		}

		bool operator< (const SingleEndT& o) const {
			int strand1, strand2;
			uint32_t *p1, *p2;

			if (b->core.tid != o.b->core.tid) return b->core.tid < o.b->core.tid;
			if (b->core.pos != o.b->core.pos) return b->core.pos < o.b->core.pos;
			strand1 = b->core.flag & 0x0010; strand2 = o.b->core.flag & 0x0010;
			if (strand1 != strand2) return strand1 < strand2;
			if (b->core.n_cigar != o.b->core.n_cigar) return b->core.n_cigar < o.b->core.n_cigar;
			p1 = bam1_cigar(b); p2 = bam1_cigar(o.b);
			for (int i = 0; i < (int)b->core.n_cigar; i++) {
				if (*p1 != *p2) return *p1 < *p2;
				++p1; ++p2;
			}
			return false;
		}
	};

	//b is mate 1, b2 is mate 2
	struct PairedEndT {
		bam1_t *b, *b2;

		PairedEndT() { b = NULL; b2 = NULL;}

		PairedEndT(bam1_t *b, bam1_t *b2) {
			this->b = b;
			this->b2 = b2;
		}

		bool operator< (const PairedEndT& o) const {
			int strand1, strand2;
			uint32_t *p1, *p2;

			//compare b
			if (b->core.tid != o.b->core.tid) return b->core.tid < o.b->core.tid;
			if (b->core.pos != o.b->core.pos) return b->core.pos < o.b->core.pos;
			strand1 = b->core.flag & 0x0010; strand2 = o.b->core.flag & 0x0010;
			if (strand1 != strand2) return strand1 < strand2;
			if (b->core.n_cigar != o.b->core.n_cigar) return b->core.n_cigar < o.b->core.n_cigar;
			p1 = bam1_cigar(b); p2 = bam1_cigar(o.b);
			for (int i = 0; i < (int)b->core.n_cigar; i++) {
				if (*p1 != *p2) return *p1 < *p2;
				++p1; ++p2;
			}

			//compare b2
			if (b2->core.tid != o.b2->core.tid) return b2->core.tid < o.b2->core.tid;
			if (b2->core.pos != o.b2->core.pos) return b2->core.pos < o.b2->core.pos;
			strand1 = b2->core.flag & 0x0010; strand2 = o.b2->core.flag & 0x0010;
			if (strand1 != strand2) return strand1 < strand2;
			if (b2->core.n_cigar != o.b2->core.n_cigar) return b2->core.n_cigar < o.b2->core.n_cigar;
			p1 = bam1_cigar(b2); p2 = bam1_cigar(o.b2);
			for (int i = 0; i < (int)b2->core.n_cigar; i++) {
				if (*p1 != *p2) return *p1 < *p2;
				++p1; ++p2;
			}

			return false;
		}
	};

	uint8_t getMAPQ(double val) {
		double err = 1.0 - val;
		if (err <= 1e-10) return 100;
		return (uint8_t)(-10 * log10(err) + .5); // round it
	}

	void push_qname(const uint8_t* qname, int l_qname, std::vector<uint8_t>& data) {
		for (int i = 0; i < l_qname; i++) data.push_back(*(qname + i));
	}

	void push_seq(const uint8_t* seq, int readlen, char strand, std::vector<uint8_t>& data) {
		int seq_len = (readlen + 1) / 2;

		switch (strand) {
		case '+': for (int i = 0; i < seq_len; i++) data.push_back(*(seq + i)); break;
		case '-':
			uint8_t code, base;
			code = 0; base = 0;
			for (int i = 0; i < readlen; i++) {
				switch (bam1_seqi(seq, readlen - i - 1)) {
				case 1: base = 8; break;
				case 2: base = 4; break;
				case 4: base = 2; break;
				case 8: base = 1; break;
				case 15: base = 15; break;
				default: assert(false);
				}
				code |=  base << (4 * (1 - i % 2));
				if (i % 2 == 1) { data.push_back(code); code = 0; }
			}

			if (readlen % 2 == 1) { data.push_back(code); }
			break;
		default: assert(false);
		}
	}

	void push_qual(const uint8_t* qual, int readlen, char strand, std::vector<uint8_t>& data) {
		switch (strand) {
		case '+': for (int i = 0; i < readlen; i++) data.push_back(*(qual + i)); break;
		case '-': for (int i = readlen - 1; i >= 0; i--) data.push_back(*(qual + i)); break;
		default: assert(false);
		}
	}

	//convert transcript coordinate to chromosome coordinate and generate CIGAR string
	void tr2chr(const Transcript&, int, int, int&, int&, std::vector<uint8_t>&);
};

//fn_list can be NULL
BamWriter::BamWriter(char inpType, const char* inpF, const char* fn_list, const char* outF, const char* chr_list) {
	switch(inpType) {
	case 's': in = samopen(inpF, "r", fn_list); break;
	case 'b': in = samopen(inpF, "rb", fn_list); break;
	default: assert(false);
	}
	assert(in != 0);

	//generate output's header
	bam_header_t *out_header = NULL;
	refmap.clear();

	if (chr_list == NULL) {
		out_header = in->header;
	}
	else {
		out_header = sam_header_read2(chr_list);

		for (int i = 0; i < out_header->n_targets; i++) {
			refmap[out_header->target_name[i]] = i;
		}
	}


	out = samopen(outF, "wb", out_header);
	assert(out != 0);

	if (chr_list != NULL) { bam_header_destroy(out_header); }
}

BamWriter::~BamWriter() {
	samclose(in);
	samclose(out);
}

void BamWriter::work(HitWrapper<SingleHit> wrapper, Transcripts& transcripts) {
	bam1_t *b;
	std::string cqname; // cqname : current query name
	std::map<SingleEndT, double> hmap;
	std::map<SingleEndT, double>::iterator hmapIter;
	SingleHit *hit;

	int cnt = 0;

	cqname = "";
	b = bam_init1();
	hmap.clear();

	while (samread(in, b) >= 0) {

		if (verbose && cnt > 0 && cnt % 1000000 == 0) { printf("%d entries are finished!\n", cnt); }
		++cnt;

		if (b->core.flag & 0x0004) continue;

		hit = wrapper.getNextHit();
		assert(hit != NULL);

		int sid = b->core.tid + 1;
		assert(sid == hit->getSid());
		const Transcript& transcript = transcripts.getTranscriptAt(sid);

		if (transcripts.getType() == 0) {
			int pos = b->core.pos;
			int readlen = b->core.l_qseq;
			uint8_t *qname = b->data, *seq = bam1_seq(b), *qual = bam1_qual(b);
			std::vector<uint8_t> data;
			data.clear();

			iter = refmap.find(transcript.getSeqName());
			assert(iter != refmap.end());
			b->core.tid = iter->second;
			b->core.qual = 255;

			uint16_t rstrand = b->core.flag & 0x0010; // read strand
			b->core.flag -= rstrand;
			rstrand = (((!rstrand && transcript.getStrand() == '+') || (rstrand && transcript.getStrand() == '-')) ? 0 : 0x0010);
			b->core.flag += rstrand;

			push_qname(qname, b->core.l_qname, data);
			int core_pos, core_n_cigar;
			tr2chr(transcript, pos + 1, pos + readlen, core_pos, core_n_cigar, data);
			if (core_pos < 0) b->core.tid = -1;
			b->core.pos = core_pos;
			b->core.n_cigar = core_n_cigar;
			push_seq(seq, readlen, transcript.getStrand(), data);
			push_qual(qual, readlen, transcript.getStrand(), data);

			free(b->data);
			b->m_data = b->data_len = data.size() + 7; // 7 extra bytes for ZW tag
			b->l_aux = 7;
			b->data = (uint8_t*)malloc(b->m_data);
			for (int i = 0; i < b->data_len; i++) b->data[i] = data[i];

			b->core.bin = bam_reg2bin(b->core.pos, bam_calend(&(b->core), bam1_cigar(b)));
		}
		else {
			b->m_data = b->data_len = b->data_len - b->l_aux + 7; // 7 extra bytes for ZW tag
			b->l_aux = 7;
			b->data = (uint8_t*)realloc(b->data, b->m_data);
		}


		if (cqname != bam1_qname(b)) {
			if (!hmap.empty()) {
				for (hmapIter = hmap.begin(); hmapIter != hmap.end(); hmapIter++) {
					bam1_t *tmp_b = hmapIter->first.b;
					tmp_b->core.qual = getMAPQ(hmapIter->second);
					uint8_t *p = bam1_aux(tmp_b);
					*p = 'Z'; ++p; *p = 'W'; ++p; *p = 'f'; ++p;
					float val = (float)hmapIter->second;
					memcpy(p, &val, 4);
					samwrite(out, tmp_b);
					bam_destroy1(tmp_b); // now hmapIter->b makes no sense
				}
				hmap.clear();
			}
			cqname = bam1_qname(b);
		}

		hmapIter = hmap.find(SingleEndT(b));
		if (hmapIter == hmap.end()) {
			hmap[SingleEndT(bam_dup1(b))] = hit->getConPrb();
		}
		else {
			hmapIter->second += hit->getConPrb();
		}
	}

	assert(wrapper.getNextHit() == NULL);

	if (!hmap.empty()) {
		for (hmapIter = hmap.begin(); hmapIter != hmap.end(); hmapIter++) {
			bam1_t *tmp_b = hmapIter->first.b;
			tmp_b->core.qual = getMAPQ(hmapIter->second);
			uint8_t *p = bam1_aux(tmp_b);
			*p = 'Z'; ++p; *p = 'W'; ++p; *p = 'f'; ++p;
			float val = (float)hmapIter->second;
			memcpy(p, &val, 4);
			samwrite(out, tmp_b);
			bam_destroy1(tmp_b); // now hmapIter->b makes no sense
		}
		hmap.clear();
	}

	bam_destroy1(b);
	if (verbose) { printf("Bam output file is generated!\n"); }
}

void BamWriter::work(HitWrapper<PairedEndHit> wrapper, Transcripts& transcripts) {
	bam1_t *b, *b2;
	std::string cqname; // cqname : current query name
	std::map<PairedEndT, double> hmap;
	std::map<PairedEndT, double>::iterator hmapIter;
	PairedEndHit *hit;

	int cnt = 0;

	cqname = "";
	b = bam_init1();
	b2 = bam_init1();
	hmap.clear();

	while (samread(in, b) >= 0 && samread(in, b2) >= 0) {

		if (verbose && cnt > 0 && cnt % 1000000 == 0) { printf("%d entries are finished!\n", cnt); }
		++cnt;

		if (!((b->core.flag & 0x0002) && (b2->core.flag & 0x0002))) continue;

		//swap if b is mate 2
		if (b->core.flag & 0x0080) {
			assert(b2->core.flag & 0x0040);
			bam1_t *tmp = b;
			b = b2; b2 = tmp;
		}

		hit = wrapper.getNextHit();
		assert(hit != NULL);

		int sid = b->core.tid + 1;
		assert(sid == hit->getSid());
		assert(sid == b2->core.tid + 1);
		const Transcript& transcript = transcripts.getTranscriptAt(sid);

		if (transcripts.getType() == 0) {
			int pos = b->core.pos, pos2 = b2->core.pos;
			int readlen = b->core.l_qseq, readlen2 = b2->core.l_qseq;
			uint8_t *qname = b->data, *seq = bam1_seq(b), *qual = bam1_qual(b);
			uint8_t *qname2 = b2->data, *seq2 = bam1_seq(b2), *qual2 = bam1_qual(b2);
			std::vector<uint8_t> data, data2;

			data.clear();
			data2.clear();

			iter = refmap.find(transcript.getSeqName());
			assert(iter != refmap.end());
			b->core.tid = iter->second; b->core.mtid = iter->second;
			b2->core.tid = iter->second; b2->core.mtid = iter->second;

			uint16_t rstrand = b->core.flag & 0x0010;
			b->core.flag = b->core.flag - (b->core.flag & 0x0010) - (b->core.flag & 0x0020);
			b2->core.flag = b2->core.flag - (b2->core.flag & 0x0010) - (b2->core.flag & 0x0020);

			uint16_t add, add2;
			if ((!rstrand && transcript.getStrand() == '+') || (rstrand && transcript.getStrand() == '-')) {
				add = 0x0020; add2 = 0x0010;
			}
			else {
				add = 0x0010; add2 = 0x0020;
			}
			b->core.flag += add;
			b2->core.flag += add2;

			b->core.qual = b2->core.qual = 255;

			//Do I really need this? The insert size uses transcript coordinates
			if (transcript.getStrand() == '-') {
				b->core.isize = -b->core.isize;
				b2->core.isize = -b2->core.isize;
			}

			push_qname(qname, b->core.l_qname, data);
			push_qname(qname2, b2->core.l_qname, data2);
			int core_pos, core_n_cigar;
			tr2chr(transcript, pos + 1, pos + readlen, core_pos, core_n_cigar, data);
			if (core_pos < 0) b->core.tid = -1;
			b->core.pos = core_pos; b->core.n_cigar = core_n_cigar;
			tr2chr(transcript, pos2 + 1, pos2 + readlen2, core_pos, core_n_cigar, data2);
			if (core_pos < 0) b2->core.tid = -1;
			b2->core.pos = core_pos; b2->core.n_cigar = core_n_cigar;
			b->core.mpos = b2->core.pos;
			b2->core.mpos = b->core.pos;
			push_seq(seq, readlen, transcript.getStrand(), data);
			push_seq(seq2, readlen2, transcript.getStrand(), data2);
			push_qual(qual, readlen, transcript.getStrand(), data);
			push_qual(qual2, readlen2, transcript.getStrand(), data2);

			free(b->data);
			b->m_data = b->data_len = data.size() + 7; // 7 extra bytes for ZW tag
			b->l_aux = 7;
			b->data = (uint8_t*)malloc(b->m_data);
			for (int i = 0; i < b->data_len; i++) b->data[i] = data[i];

			free(b2->data);
			b2->m_data = b2->data_len = data2.size() + 7; // 7 extra bytes for ZW tag
			b2->l_aux = 7;
			b2->data = (uint8_t*)malloc(b2->m_data);
			for (int i = 0; i < b2->data_len; i++) b2->data[i] = data2[i];

			b->core.bin = bam_reg2bin(b->core.pos, bam_calend(&(b->core), bam1_cigar(b)));
			b2->core.bin = bam_reg2bin(b2->core.pos, bam_calend(&(b2->core), bam1_cigar(b2)));
		}
		else {
			b->m_data = b->data_len = b->data_len - b->l_aux + 7; // 7 extra bytes for ZW tag
			b->l_aux = 7;
			b->data = (uint8_t*)realloc(b->data, b->m_data);

			b2->m_data = b2->data_len = b2->data_len - b2->l_aux + 7; // 7 extra bytes for ZW tag
			b2->l_aux = 7;
			b2->data = (uint8_t*)realloc(b2->data, b2->m_data);
		}

		if (cqname != bam1_qname(b)) {
			if (!hmap.empty()) {
				for (hmapIter = hmap.begin(); hmapIter != hmap.end(); hmapIter++) {
					bam1_t *tmp_b = hmapIter->first.b;
					bam1_t *tmp_b2 = hmapIter->first.b2;

					tmp_b->core.qual = tmp_b2->core.qual = getMAPQ(hmapIter->second);

					uint8_t *p = bam1_aux(tmp_b), *p2 = bam1_aux(tmp_b2);
					*p = 'Z'; ++p; *p = 'W'; ++p; *p = 'f'; ++p;
					*p2 = 'Z'; ++p2; *p2 = 'W'; ++p2; *p2 = 'f'; ++p2;

					float val = (float)hmapIter->second;
					memcpy(p, &val, 4);
					memcpy(p2, &val, 4);

					samwrite(out, tmp_b);
					samwrite(out, tmp_b2);

					bam_destroy1(tmp_b);
					bam_destroy1(tmp_b2);
				}
				hmap.clear();
			}
			cqname = bam1_qname(b);
		}

		hmapIter = hmap.find(PairedEndT(b, b2));
		if (hmapIter == hmap.end()) {
			hmap[PairedEndT(bam_dup1(b), bam_dup1(b2))] = hit->getConPrb();
		}
		else {
			hmapIter->second += hit->getConPrb();
		}
	}

	assert(wrapper.getNextHit() == NULL);

	if (!hmap.empty()) {
		for (hmapIter = hmap.begin(); hmapIter != hmap.end(); hmapIter++) {
			bam1_t *tmp_b = hmapIter->first.b;
			bam1_t *tmp_b2 = hmapIter->first.b2;

			tmp_b->core.qual = tmp_b2->core.qual = getMAPQ(hmapIter->second);

			uint8_t *p = bam1_aux(tmp_b), *p2 = bam1_aux(tmp_b2);
			*p = 'Z'; ++p; *p = 'W'; ++p; *p = 'f'; ++p;
			*p2 = 'Z'; ++p2; *p2 = 'W'; ++p2; *p2 = 'f'; ++p2;

			float val = (float)hmapIter->second;
			memcpy(p, &val, 4);
			memcpy(p2, &val, 4);

			samwrite(out, tmp_b);
			samwrite(out, tmp_b2);

			bam_destroy1(tmp_b);
			bam_destroy1(tmp_b2);
		}
		hmap.clear();
	}

	bam_destroy1(b);
	bam_destroy1(b2);

	if (verbose) { printf("Bam output file is generated!"); }
}

void BamWriter::tr2chr(const Transcript& transcript, int sp, int ep, int& pos, int& n_cigar, std::vector<uint8_t>& data) {
	int length = transcript.getLength();
	char strand = transcript.getStrand();
	const std::vector<Interval>& structure = transcript.getStructure();

	int s, i;
	int oldlen, curlen;

	uint32_t operation;
	uint8_t *p;

	n_cigar = 0;
	s = structure.size();

	if (strand == '-') {
		int tmp = sp;
		sp = length - ep + 1;
		ep = length - tmp + 1;
	}

	if (ep < 1 || sp > length) { // a read which align to polyA tails totally! 
	  pos = (sp > length ? structure[s - 1].end : structure[0].start - 1); // 0 based

	  n_cigar = 1;
	  operation = (ep - sp + 1) << BAM_CIGAR_SHIFT | BAM_CINS; //BAM_CSOFT_CLIP;
	  p = (uint8_t*)(&operation);
	  for (int j = 0; j < 4; j++) data.push_back(*(p + j));

	  return;
	}

	if (sp < 1) {
		n_cigar++;
		operation = (1 - sp) << BAM_CIGAR_SHIFT | BAM_CINS; //BAM_CSOFT_CLIP;
		p = (uint8_t*)(&operation);
		for (int j = 0; j < 4; j++) data.push_back(*(p + j));
		sp = 1;
	}

	oldlen = curlen = 0;

	for (i = 0; i < s; i++) {
		oldlen = curlen;
		curlen += structure[i].end - structure[i].start + 1;
		if (curlen >= sp) break;
	}
	assert(i < s);
	pos = structure[i].start + (sp - oldlen - 1) - 1; // 0 based

	while (curlen < ep && i < s) {
		n_cigar++;
		operation = (curlen - sp + 1) << BAM_CIGAR_SHIFT | BAM_CMATCH;
		p = (uint8_t*)(&operation);
		for (int j = 0; j < 4; j++) data.push_back(*(p + j));

		++i;
		if (i >= s) continue;
		n_cigar++;
		operation = (structure[i].start - structure[i - 1].end - 1) << BAM_CIGAR_SHIFT | BAM_CREF_SKIP;
		p = (uint8_t*)(&operation);
		for (int j = 0; j < 4; j++) data.push_back(*(p + j));

		oldlen = curlen;
		sp = oldlen + 1;
		curlen += structure[i].end - structure[i].start + 1;
	}

	if (i >= s) {
		n_cigar++;
		operation = (ep - length) << BAM_CIGAR_SHIFT | BAM_CINS; //BAM_CSOFT_CLIP;
		p = (uint8_t*)(&operation);
		for (int j = 0; j < 4; j++) data.push_back(*(p + j));
	}
	else {
		n_cigar++;
		operation = (ep - sp + 1) << BAM_CIGAR_SHIFT | BAM_CMATCH;
		p = (uint8_t*)(&operation);
		for (int j = 0; j < 4; j++) data.push_back(*(p + j));
	}
}

#endif /* BAMWRITER_H_ */
