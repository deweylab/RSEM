#include <cstring>
#include <cstdlib>
#include <cassert>

#include "wiggle.h"

#include "sam/bam.h"
#include "sam/sam.h"

void add_bam_record_to_wiggle(const bam1_t *b, Wiggle& wiggle) {
    float w = bam_aux2f(bam_aux_get(b, "ZW"));
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
	//assert(bam_in != 0);

    int cur_tid = -1; //current tid;
	int cnt = 0;
    bam1_t *b = bam_init1();
    Wiggle wiggle;
	while (samread(bam_in, b) >= 0) {
		if (b->core.tid != cur_tid) {
			if (cur_tid >= 0) processor.process(wiggle);
			cur_tid = b->core.tid;
            wiggle.name = bam_in->header->target_name[cur_tid];
            wiggle.read_depth.assign(bam_in->header->target_len[cur_tid], 0.0);
		}
        add_bam_record_to_wiggle(b, wiggle);
		++cnt;
		if (cnt % 1000000 == 0) fprintf(stderr, "%d FIN\n", cnt);
	}
	if (cur_tid >= 0) processor.process(wiggle);

	samclose(bam_in);
	bam_destroy1(b);
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
    
    sp = ep = -1;
    for (size_t i = 0; i < wiggle.read_depth.size(); i++) {
        if (wiggle.read_depth[i] > 0) {
            ep = i;
        }
        else {
            if (sp < ep) {
                ++sp;
                fprintf(fo, "fixedStep chrom=%s start=%d step=1\n", wiggle.name.c_str(), sp + 1);
                for (int j = sp; j <= ep; j++) fprintf(fo, "%.7g\n", wiggle.read_depth[j]);
            }
            sp = i;
        }
    }
    if (sp < ep) {
        ++sp;
        fprintf(fo, "fixedStep chrom=%s start=%d step=1\n", wiggle.name.c_str(), sp + 1);
        for (int j = sp; j <= ep; j++) fprintf(fo, "%.7g\n", wiggle.read_depth[j]);
    }
}

ReadDepthWriter::ReadDepthWriter(std::ostream& stream) 
    : stream_(stream) {
}

void ReadDepthWriter::process(const Wiggle& wiggle) {
    stream_ << wiggle.name << '\t'
            << wiggle.read_depth.size() << '\t';
    for (size_t i = 0; i < wiggle.read_depth.size(); ++i) {
        if (i > 0) stream_ << ' ';
        stream_ << wiggle.read_depth[i];
    }
    stream_ << '\n';
}
