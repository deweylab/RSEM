#ifndef SAM_UTILS_H_
#define SAM_UTILS_H_

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<vector>
#include<string>
#include<fstream>

#include<stdint.h>

#include "htslib/sam.h"

#include "my_assert.h"
#include "Transcript.h"
#include "Transcripts.h"


/******************************************************/

// These functions are adopted/modified from samtools source codes because the original codes are not visible from sam.h/bam.h

inline int bam_aux_type2size(char x) {
  if (x == 'C' || x == 'c' || x == 'A') return 1;
  else if (x == 'S' || x == 's') return 2;
  else if (x == 'I' || x == 'i' || x == 'f') return 4;
  else if (x == 'd') return 8;
  else return 0;
}

inline void expand_data_size(bam1_t *b) {
  if (b->m_data < b->l_data) {
    b->m_data = b->l_data;
    kroundup32(b->m_data);
    b->data = (uint8_t*)realloc(b->data, b->m_data);
  }
}

/******************************************************/

// These functions are specially designed for RSEM

const char* whitespaces = " \t\n\r\f\v";

inline bool bam_is_paired(const bam1_t* b) { return (b->core.flag & BAM_FPAIRED); }
inline bool bam_is_proper(const bam1_t* b) { return (b->core.flag & BAM_FPROPER_PAIR); }
inline bool bam_is_mapped(const bam1_t* b) { return !(b->core.flag & BAM_FUNMAP); }
inline bool bam_is_unmapped(const bam1_t* b) { return (b->core.flag & BAM_FUNMAP); }
inline bool bam_is_read1(const bam1_t* b) { return (b->core.flag & BAM_FREAD1); }
inline bool bam_is_read2(const bam1_t* b) { return (b->core.flag & BAM_FREAD2); }

inline std::string bam_get_canonical_name(const bam1_t* b) {
  // Retain only the first whitespace-delimited word as the read name
  // This prevents issues of mismatching names when aligners do not
  // strip off extra words in read name strings
  const char* raw_query_name = bam_get_qname(b);
  const char* whitespace_pos = std::strpbrk(raw_query_name, whitespaces);
  return (whitespace_pos == NULL ? std::string(raw_query_name) : std::string(raw_query_name, whitespace_pos - raw_query_name));
}

// Current RSEM only accept matches
inline bool bam_check_cigar(bam1_t *b) {
  uint32_t *cigar = bam_get_cigar(b);
  char op = bam_cigar_op(*cigar);
  int32_t oplen = bam_cigar_oplen(*cigar);
  
  return (b->core.n_cigar == 1) && (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) && (b->core.l_qseq == oplen);
}

uint8_t bam_prb_to_mapq(double val) {
  double err = 1.0 - val;
  if (err <= 1e-10) return 100;
  return (uint8_t)(-10 * log10(err) + .5); // round it
}

inline std::string bam_get_read_seq(const bam1_t* b) {
  uint8_t *p = bam_get_seq(b);
  std::string readseq = "";
  char base = 0;

  if (bam_is_rev(b)) {
    for (int i = b->core.l_qseq - 1; i >= 0; i--) {
      switch(bam_seqi(p, i)) {
	//case 0 : base = '='; break;
      case 1 : base = 'T'; break;
      case 2 : base = 'G'; break;
      case 4 : base = 'C'; break;
      case 8 : base = 'A'; break;
      case 15 : base = 'N'; break;
      default : assert(false);
      }
      readseq.append(1, base);
    }
  }
  else {
    for (int i = 0; i < b->core.l_qseq; ++i) {
      switch(bam_seqi(p, i)) {
	//case 0 : base = '='; break;
      case 1 : base = 'A'; break;
      case 2 : base = 'C'; break;
      case 4 : base = 'G'; break;
      case 8 : base = 'T'; break;
      case 15 : base = 'N'; break;
      default : assert(false);
      }
      readseq.append(1, base);
    }
  }
  
  return readseq;
}

inline std::string bam_get_qscore(const bam1_t* b) {
  uint8_t *p = bam_get_qual(b);
  std::string qscore = "";
  
  if (bam_is_rev(b)) {
    p = p + b->core.l_qseq - 1;
    for (int i = 0; i < b->core.l_qseq; ++i) {
      qscore.append(1, (char)(*p + 33));
      --p;
    }
  }
  else {
    for (int i = 0; i < b->core.l_qseq; ++i) {
      qscore.append(1, (char)(*p + 33));
      ++p;
    }
  }
  
  return qscore;
}

//convert transcript coordinate to chromosome coordinate and generate CIGAR string
void tr2chr(const Transcript& transcript, int sp, int ep, int& pos, int& n_cigar, std::vector<uint32_t>& data) {
	int length = transcript.getLength();
	char strand = transcript.getStrand();
	const std::vector<Interval>& structure = transcript.getStructure();

	int s, i;
	int oldlen, curlen;

	uint32_t operation;

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
	  data.push_back(operation);

	  return;
	}

	if (sp < 1) {
		n_cigar++;
		operation = (1 - sp) << BAM_CIGAR_SHIFT | BAM_CINS; //BAM_CSOFT_CLIP;
		data.push_back(operation);
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
		data.push_back(operation);
		++i;
		if (i >= s) continue;
		n_cigar++;
		operation = (structure[i].start - structure[i - 1].end - 1) << BAM_CIGAR_SHIFT | BAM_CREF_SKIP;
		data.push_back(operation);

		oldlen = curlen;
		sp = oldlen + 1;
		curlen += structure[i].end - structure[i].start + 1;
	}

	if (i >= s) {
		n_cigar++;
		operation = (ep - length) << BAM_CIGAR_SHIFT | BAM_CINS; //BAM_CSOFT_CLIP;
		data.push_back(operation);
	}
	else {
		n_cigar++;
		operation = (ep - sp + 1) << BAM_CIGAR_SHIFT | BAM_CMATCH;
		data.push_back(operation);
	}
}

#endif /* SAM_RSEM_AUX_H_ */
