#ifndef SAM_RSEM_CVT_H_
#define SAM_RSEM_CVT_H_

#include<vector>

#include "stdint.h"
#include "sam/bam.h"

#include "Transcript.h"
#include "Transcripts.h"

uint8_t getMAPQ(double val) {
	double err = 1.0 - val;
	if (err <= 1e-10) return 100;
	return (uint8_t)(-10 * log10(err) + .5); // round it
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

#endif /* SAM_RSEM_CVT_H_ */
