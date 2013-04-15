#include <cstring>
#include <cstdlib>
#include <cassert>
#include <iostream>

#include <stdint.h>
#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"
#include "wiggle.h"

bool no_fractional_weight = false;

void add_bam_record_to_wiggle(const bam1_t *b, Wiggle& wiggle) {
    double w;

    if (no_fractional_weight) w = 1.0;
    else {
      uint8_t *p_tag = bam_aux_get(b, "ZW");
      if (p_tag == NULL) return;
      w = bam_aux2f(p_tag);
    }

    int pos = b->core.pos;
    uint32_t *p = bam1_cigar(b);
    
    for (int i = 0; i < (int)b->core.n_cigar; i++, ++p) {
        int op = *p & BAM_CIGAR_MASK;
        int op_len = *p >> BAM_CIGAR_SHIFT;
        
        switch (op) {
            //case BAM_CSOFT_CLIP : pos += op_len; break;
        case BAM_CINS : pos += op_len; break;
        case BAM_CMATCH :
            for (int j = 0; j < op_len; j++, ++pos) {
                wiggle.read_depth[pos] += w;
            }
            break;
        case BAM_CREF_SKIP : pos += op_len; break;
        default : assert(false);
        }
    }
}

void build_wiggles(const std::string& bam_filename,
                   WiggleProcessor& processor) {
    samfile_t *bam_in = samopen(bam_filename.c_str(), "rb", NULL);
	if (bam_in == 0) { fprintf(stderr, "Cannot open %s!\n", bam_filename.c_str()); exit(-1); }

	bam_header_t *header = bam_in->header;
	bool *used = new bool[header->n_targets];
	memset(used, 0, sizeof(bool) * header->n_targets);

	int cur_tid = -1; //current tid;
	HIT_INT_TYPE cnt = 0;
	bam1_t *b = bam_init1();
	Wiggle wiggle;
	while (samread(bam_in, b) >= 0) {
		if (b->core.flag & 0x0004) continue;

		if (b->core.tid != cur_tid) {
			if (cur_tid >= 0) { used[cur_tid] = true; processor.process(wiggle); }
			cur_tid = b->core.tid;
			wiggle.name = header->target_name[cur_tid];
			wiggle.length = header->target_len[cur_tid];
			wiggle.read_depth.assign(wiggle.length, 0.0);
		}
		add_bam_record_to_wiggle(b, wiggle);
		++cnt;
		if (cnt % 1000000 == 0) std::cout<< cnt<< std::endl;
	}
	if (cur_tid >= 0) { used[cur_tid] = true; processor.process(wiggle); }

	for (int32_t i = 0; i < header->n_targets; i++)
		if (!used[i]) {
			wiggle.name = header->target_name[i];
			wiggle.length = header->target_len[i];
			wiggle.read_depth.clear();
			processor.process(wiggle);
		}

	samclose(bam_in);
	bam_destroy1(b);
	delete[] used;
}

UCSCWiggleTrackWriter::UCSCWiggleTrackWriter(const std::string& output_filename,
                                             const std::string& track_name) {
    fo = fopen(output_filename.c_str(), "w");
    fprintf(fo, "track type=wiggle_0 name=\"%s\" description=\"%s\" visibility=full\n",
            track_name.c_str(),
            track_name.c_str());
}

UCSCWiggleTrackWriter::~UCSCWiggleTrackWriter() {
    fclose(fo);
}

void UCSCWiggleTrackWriter::process(const Wiggle& wiggle) {
    int sp, ep;

    if (wiggle.read_depth.empty()) return;
    
    sp = ep = -1;
    for (size_t i = 0; i < wiggle.length; i++) {
        if (wiggle.read_depth[i] >= 0.0095) {
            ep = i;
        }
        else {
            if (sp < ep) {
                ++sp;
                fprintf(fo, "fixedStep chrom=%s start=%d step=1\n", wiggle.name.c_str(), sp + 1);
                for (int j = sp; j <= ep; j++) fprintf(fo, "%.2f\n", wiggle.read_depth[j]);
            }
            sp = i;
        }
    }
    if (sp < ep) {
        ++sp;
        fprintf(fo, "fixedStep chrom=%s start=%d step=1\n", wiggle.name.c_str(), sp + 1);
        for (int j = sp; j <= ep; j++) fprintf(fo, "%.2f\n", wiggle.read_depth[j]);
    }
}

ReadDepthWriter::ReadDepthWriter(std::ostream& stream) 
    : stream_(stream) {
}

void ReadDepthWriter::process(const Wiggle& wiggle) {

    stream_ << wiggle.name << '\t'
            << wiggle.length << '\t';

    if (wiggle.read_depth.empty()) { stream_ << "NA\n"; return; }

    for (size_t i = 0; i < wiggle.length; ++i) {
        if (i > 0) stream_ << ' ';
        stream_ << wiggle.read_depth[i];
    }
    stream_ << '\n';
}
