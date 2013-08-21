#ifndef BAMWRITER_H_
#define BAMWRITER_H_

#include<cmath>
#include<cstdio>
#include<cstring>
#include<cassert>
#include<string>
#include<sstream>
#include<iostream>

#include <stdint.h>
#include "sam/bam.h"
#include "sam/sam.h"
#include "sam_rsem_aux.h"
#include "sam_rsem_cvt.h"

#include "utils.h"

#include "SingleHit.h"
#include "PairedEndHit.h"

#include "HitWrapper.h"
#include "Transcript.h"
#include "Transcripts.h"

class BamWriter {
public:
	BamWriter(char, const char*, const char*, const char*, Transcripts&);
	~BamWriter();

	void work(HitWrapper<SingleHit>);
	void work(HitWrapper<PairedEndHit>);
private:
	samfile_t *in, *out;
	Transcripts& transcripts;

	//convert bam1_t
	void convert(bam1_t*, double);
};

//fn_list can be NULL
BamWriter::BamWriter(char inpType, const char* inpF, const char* fn_list, const char* outF, Transcripts& transcripts)
	: transcripts(transcripts)
{
	switch(inpType) {
	case 's': in = samopen(inpF, "r", fn_list); break;
	case 'b': in = samopen(inpF, "rb", fn_list); break;
	default: assert(false);
	}
	assert(in != 0);

	//build mappings from external sid to internal sid
	transcripts.buildMappings(in->header->n_targets, in->header->target_name);

	//generate output's header
	bam_header_t *out_header = bam_header_dwt(in->header);

	std::ostringstream strout;
	strout<<"@HD\tVN:1.4\tSO:unknown\n@PG\tID:RSEM\n";
	std::string content = strout.str();
	append_header_text(out_header, content.c_str(), content.length());

	out = samopen(outF, "wb", out_header);
	assert(out != 0);

	bam_header_destroy(out_header);
}

BamWriter::~BamWriter() {
	samclose(in);
	samclose(out);
}

void BamWriter::work(HitWrapper<SingleHit> wrapper) {
	bam1_t *b;
	SingleHit *hit;

	HIT_INT_TYPE cnt = 0;

	b = bam_init1();

	while (samread(in, b) >= 0) {
		++cnt;
		if (verbose && cnt % 1000000 == 0) { std::cout<< cnt<< " alignment lines are loaded!"<< std::endl; }

		if (!(b->core.flag & 0x0004)) {
		  hit = wrapper.getNextHit();
		  assert(hit != NULL);

		  assert(transcripts.getInternalSid(b->core.tid + 1) == hit->getSid());
		  convert(b, hit->getConPrb());
		}

		//if (b->core.qual > 0) samwrite(out, b); // output only when MAPQ > 0
		samwrite(out, b);
	}

	assert(wrapper.getNextHit() == NULL);

	bam_destroy1(b);
	if (verbose) { std::cout<< "Bam output file is generated!"<< std::endl; }
}

void BamWriter::work(HitWrapper<PairedEndHit> wrapper) {
	bam1_t *b, *b2;
	PairedEndHit *hit;

	HIT_INT_TYPE cnt = 0;

	b = bam_init1();
	b2 = bam_init1();

	while (samread(in, b) >= 0 && samread(in, b2) >= 0) {
		cnt += 2;
		if (verbose && cnt % 1000000 == 0) { std::cout<< cnt<< " alignment lines are loaded!"<< std::endl; }

		assert((b->core.flag & 0x0001) && (b2->core.flag & 0x0001));
		assert(((b->core.flag & 0x0040) && (b2->core.flag & 0x0080)) || ((b->core.flag & 0x0080) && (b2->core.flag & 0x0040)));

		//unalignable reads, skip		
		bool notgood = (b->core.flag & 0x0004) || (b2->core.flag & 0x0004);

		if (!notgood) {
			//swap if b is mate 2
			if (b->core.flag & 0x0080) {
				assert(b2->core.flag & 0x0040);
				bam1_t *tmp = b;
				b = b2; b2 = tmp;
			}

			hit = wrapper.getNextHit();
			assert(hit != NULL);

			assert(transcripts.getInternalSid(b->core.tid + 1) == hit->getSid());
			assert(transcripts.getInternalSid(b2->core.tid + 1) == hit->getSid());

			convert(b, hit->getConPrb());
			convert(b2, hit->getConPrb());

			// omit "b->core.mpos = b2->core.pos; b2->core.mpos = b->core.pos" for the reason that it is possible that only one mate is alignable 
		}

		samwrite(out, b);
		samwrite(out, b2);
	}

	assert(wrapper.getNextHit() == NULL);

	bam_destroy1(b);
	bam_destroy1(b2);

	if (verbose) { std::cout<< "Bam output file is generated!"<< std::endl; }
}

void BamWriter::convert(bam1_t *b, double prb) {
	b->core.qual = getMAPQ(prb);
	float val = (float)prb;
	bam_aux_append(b, "ZW", 'f', bam_aux_type2size('f'), (uint8_t*)&val);
}

#endif /* BAMWRITER_H_ */
