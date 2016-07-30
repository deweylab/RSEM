#include <cstring>
#include <cstdlib>
#include <cassert>
#include <iostream>

#include <stdint.h>
#include "htslib/sam.h"
#include "sam_utils.h"

#include "utils.h"
#include "my_assert.h"
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
    uint32_t *p = bam_get_cigar(b);
    
    for (int i = 0; i < (int)b->core.n_cigar; ++i, ++p) {
      char op = bam_cigar_op(*p);
      int op_len = bam_cigar_oplen(*p);

      if (op == BAM_CMATCH)
	for (int j = 0; j < op_len; ++j, ++pos) wiggle.read_depth[pos] += w;
      else pos += ((bam_cigar_type(op) & 2) ? op_len : 0);
    }
}

void build_wiggles(const std::string& bam_filename,
                   WiggleProcessor& processor) {
  
    samFile *bam_in = sam_open(bam_filename.c_str(), "r");
    general_assert(bam_in != NULL, "Cannot open " + bam_filename + "!");

    bam_hdr_t *header = sam_hdr_read(bam_in);
    general_assert(header != 0, "Cannot load SAM header!");
    bool *used = new bool[header->n_targets];
    memset(used, 0, sizeof(bool) * header->n_targets);

    int cur_tid = -1; //current tid;
    HIT_INT_TYPE cnt = 0;
    bam1_t *b = bam_init1();
    Wiggle wiggle;
    while (sam_read1(bam_in, header, b) >= 0) {
      if (bam_is_unmapped(b)) continue;
      
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

    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(bam_in);

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
