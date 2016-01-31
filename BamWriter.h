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
#include "bam.h"
#include "sam.h"
#include "sam_utils.h"

#include "utils.h"

#include "SingleHit.h"
#include "PairedEndHit.h"

#include "HitWrapper.h"
#include "Transcript.h"
#include "Transcripts.h"

class BamWriter {
public:
	BamWriter(const char* inpF, const char* aux, const char* outF, Transcripts& transcripts);
	~BamWriter();

	void work(HitWrapper<SingleHit> wrapper);
	void work(HitWrapper<PairedEndHit> wrapper);
private:
	samfile_t *in, *out;
	Transcripts& transcripts;
	
	void set_alignment_weight(bam1_t *b, double prb) {
	  b->core.qual = bam_prb_to_mapq(prb);
	  float val = (float)prb;
	  bam_aux_append(b, "ZW", 'f', bam_aux_type2size('f'), (uint8_t*)&val);
	}
};

//fn_list can be NULL
BamWriter::BamWriter(const char* inpF, const char* aux, const char* outF, Transcripts& transcripts) : transcripts(transcripts) {
  in = samopen(inpF, "r", aux);
  assert(in != 0);

  //build mappings from external sid to internal sid
  transcripts.buildMappings(in->header->n_targets, in->header->target_name);
  
  //generate output's header
  bam_hdr_t *out_header = bam_header_dwt(in->header);
  
  std::ostringstream strout;
  strout<<"@HD\tVN:1.4\tSO:unknown\n@PG\tID:RSEM\n";
  std::string content = strout.str();
  append_header_text(out_header, content.c_str(), content.length());
  
  out = samopen(outF, "wb", out_header); // If CRAM format is desired, use "wc"
  assert(out != 0);
  
  bam_hdr_destroy(out_header);
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

		if (bam_is_mapped(b)) {
		  hit = wrapper.getNextHit();
		  assert(hit != NULL);

		  assert(transcripts.getInternalSid(b->core.tid + 1) == hit->getSid());
		  set_alignment_weight(b, hit->getConPrb());
		}
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

		if (!bam_is_read1(b)) { bam1_t * tmp = b; b = b2; b2 = tmp; }
		
		if (bam_is_mapped(b) && bam_is_mapped(b2)) {
			hit = wrapper.getNextHit();
			assert(hit != NULL);

			assert(transcripts.getInternalSid(b->core.tid + 1) == hit->getSid());
			assert(transcripts.getInternalSid(b2->core.tid + 1) == hit->getSid());

			set_alignment_weight(b, hit->getConPrb());
			set_alignment_weight(b2, hit->getConPrb());
		}

		samwrite(out, b);
		samwrite(out, b2);
	}

	assert(wrapper.getNextHit() == NULL);

	bam_destroy1(b);
	bam_destroy1(b2);

	if (verbose) { std::cout<< "Bam output file is generated!"<< std::endl; }
}

#endif /* BAMWRITER_H_ */
