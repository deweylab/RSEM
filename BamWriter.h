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
#include "htslib/sam.h"
#include "sam_utils.h"
#include "SamHeader.hpp"

#include "utils.h"
#include "my_assert.h"

#include "SingleHit.h"
#include "PairedEndHit.h"

#include "HitWrapper.h"
#include "Transcript.h"
#include "Transcripts.h"

class BamWriter {
public:
	BamWriter(const char* inpF, const char* aux, const char* outF, Transcripts& transcripts, int nThreads);
	~BamWriter();

	void work(HitWrapper<SingleHit> wrapper);
	void work(HitWrapper<PairedEndHit> wrapper);
private:
	samFile *in, *out;
	bam_hdr_t *in_header, *out_header;
	Transcripts& transcripts;
	
	void set_alignment_weight(bam1_t *b, double prb) {
	  b->core.qual = bam_prb_to_mapq(prb);
	  float val = (float)prb;
          uint8_t *p = bam_aux_get(b, "ZW");
          if (p != NULL) {
            memcpy(p + 1, (uint8_t*)&(val), bam_aux_type2size('f'));
          } else {
            bam_aux_append(b, "ZW", 'f', bam_aux_type2size('f'), (uint8_t*)&val);
          }
	}
};

//aux can be NULL
BamWriter::BamWriter(const char* inpF, const char* aux, const char* outF, Transcripts& transcripts, int nThreads) : transcripts(transcripts) {
  in = sam_open(inpF, "r");
  assert(in != 0);

  if (aux == NULL) hts_set_fai_filename(in, aux);
  in_header = sam_hdr_read(in);
  assert(in_header != 0);

  //build mappings from external sid to internal sid
  transcripts.buildMappings(in_header->n_targets, in_header->target_name);
  
  //generate output's header
  SamHeader hdr(in_header->text);
  hdr.insertPG("RSEM");
  out_header = hdr.create_header();
  
  out = sam_open(outF, "wb"); // If CRAM format is desired, use "wc"
  assert(out != 0);
  sam_hdr_write(out, out_header);
    
  if (nThreads > 1) general_assert(hts_set_threads(out, nThreads) == 0, "Fail to create threads for writing the BAM file!");
}

BamWriter::~BamWriter() {
	bam_hdr_destroy(in_header);
	sam_close(in);
	bam_hdr_destroy(out_header);
	sam_close(out);
}

void BamWriter::work(HitWrapper<SingleHit> wrapper) {
	bam1_t *b;
	SingleHit *hit;

	HIT_INT_TYPE cnt = 0;

	b = bam_init1();

	while (sam_read1(in, in_header, b) >= 0) {
		++cnt;
		if (verbose && cnt % 1000000 == 0) { std::cout<< cnt<< " alignment lines are loaded!"<< std::endl; }

		if (bam_is_mapped(b)) {
		  hit = wrapper.getNextHit();
		  assert(hit != NULL);

		  assert(transcripts.getInternalSid(b->core.tid + 1) == hit->getSid());
		  set_alignment_weight(b, hit->getConPrb());
		}
		sam_write1(out, out_header, b);
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

	while (sam_read1(in, in_header, b) >= 0 && sam_read1(in, in_header, b2) >= 0) {
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

		sam_write1(out, out_header, b);
		sam_write1(out, out_header, b2);
	}

	assert(wrapper.getNextHit() == NULL);

	bam_destroy1(b);
	bam_destroy1(b2);

	if (verbose) { std::cout<< "Bam output file is generated!"<< std::endl; }
}

#endif /* BAMWRITER_H_ */
