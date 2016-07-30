#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<vector>
#include<algorithm>

#include <stdint.h>
#include "htslib/sam.h"
#include "sam_utils.h"

#include "utils.h"
#include "my_assert.h"

using namespace std;

int nThreads;
samFile *in, *out;
bam_hdr_t *header;
bam1_t *b, *b2;
vector<bam1_t*> arr_both, arr_partial_1, arr_partial_2, arr_partial_unknown;

inline void add_to_appropriate_arr(bam1_t *b) {
  if (bam_is_mapped(b) && bam_is_proper(b)) {
    arr_both.push_back(bam_dup1(b)); return;
  }

  if (bam_is_read1(b)) arr_partial_1.push_back(bam_dup1(b));
  else if (bam_is_read2(b)) arr_partial_2.push_back(bam_dup1(b));
  else arr_partial_unknown.push_back(bam_dup1(b));
}

char get_pattern_code(uint32_t flag) {
  if (flag & BAM_FREAD1) return ((flag & BAM_FREVERSE) ? 1 : 0);
  else return ((flag & BAM_FREVERSE) ? 0 : 1);
}

bool less_than(bam1_t *a, bam1_t *b) {
	int32_t ap1 = min(a->core.pos, a->core.mpos);
	int32_t ap2 = max(a->core.pos, a->core.mpos);
	int32_t bp1 = min(b->core.pos, b->core.mpos);
	int32_t bp2 = max(b->core.pos, b->core.mpos);
	char apat = get_pattern_code(a->core.flag); // apt: a's pattern of strand and mate information
	char bpat = get_pattern_code(b->core.flag);

	if (a->core.tid != b->core.tid) return a->core.tid < b->core.tid;
	if (ap1 != bp1) return ap1 < bp1;
	if (ap2 != bp2) return ap2 < bp2;
	return apat < bpat;
}

int main(int argc, char* argv[]) {
	if (argc != 4) {
		printf("Usage: rsem-scan-for-paired-end-reads number_of_threads input.[sam/bam/cram] output.bam\n");
		exit(-1);
	}

        nThreads = atoi(argv[1]);
	in = sam_open(argv[2], "r");
	general_assert(in != 0, "Cannot open " + cstrtos(argv[2]) + " !");
	header = sam_hdr_read(in);
	general_assert(header != 0, "Cannot load SAM header!");
	out = sam_open(argv[3], "wb");
	general_assert(out != 0, "Cannot open " + cstrtos(argv[3]) + " !");
	sam_hdr_write(out, header);
	if (nThreads > 1) general_assert(hts_set_threads(out, nThreads) == 0, "Fail to create threads for writing the BAM file!");

	b = bam_init1(); b2 = bam_init1();

	string qname;
	bool go_on = (sam_read1(in, header, b) >= 0);
	bool isPaired;
	HIT_INT_TYPE cnt = 0;

	printf("."); fflush(stdout);

	while (go_on) {
	  qname = bam_get_canonical_name(b);
	  isPaired = bam_is_paired(b);

	  if (isPaired) {
	    add_to_appropriate_arr(b);
	    while ((go_on = (sam_read1(in, header, b) >= 0)) && (qname == bam_get_canonical_name(b))) {
	      general_assert_1(bam_is_paired(b), "Read " + qname + " is detected as both single-end and paired-end read!");
	      add_to_appropriate_arr(b);
	    }

	    general_assert_1(arr_both.size() % 2 == 0, "Number of first and second mates in read " + qname + "'s full alignments (both mates are aligned) are not matched!");
	    general_assert_1((arr_partial_1.size() + arr_partial_2.size() + arr_partial_unknown.size()) % 2 == 0, "Number of first and second mates in read " + qname + "'s partial alignments (at most one mate is aligned) are not matched!");

	    if (!arr_both.empty()) {
	      sort(arr_both.begin(), arr_both.end(), less_than);
	      for (size_t i = 0; i < arr_both.size(); i++) { sam_write1(out, header, arr_both[i]); bam_destroy1(arr_both[i]); }
	      arr_both.clear();
	    }
	    
	    while (!arr_partial_1.empty() || !arr_partial_2.empty()) {
	      if (!arr_partial_1.empty() && !arr_partial_2.empty()) {
		sam_write1(out, header, arr_partial_1.back()); bam_destroy1(arr_partial_1.back()); arr_partial_1.pop_back();
		sam_write1(out, header, arr_partial_2.back()); bam_destroy1(arr_partial_2.back()); arr_partial_2.pop_back();
	      }
	      else if (!arr_partial_1.empty()) {
		sam_write1(out, header, arr_partial_1.back()); bam_destroy1(arr_partial_1.back()); arr_partial_1.pop_back();
		sam_write1(out, header, arr_partial_unknown.back()); bam_destroy1(arr_partial_unknown.back()); arr_partial_unknown.pop_back();
	      }
	      else {
		sam_write1(out, header, arr_partial_2.back()); bam_destroy1(arr_partial_2.back()); arr_partial_2.pop_back();
		sam_write1(out, header, arr_partial_unknown.back()); bam_destroy1(arr_partial_unknown.back()); arr_partial_unknown.pop_back();
	      }
	    }
	    
	    while (!arr_partial_unknown.empty()) {
	      sam_write1(out, header, arr_partial_unknown.back()); bam_destroy1(arr_partial_unknown.back()); arr_partial_unknown.pop_back();
	    }
	  }
	  else {
	    sam_write1(out, header, b);
	    while ((go_on = (sam_read1(in, header, b) >= 0)) && (qname == bam_get_qname(b))) {
	      sam_write1(out, header, b);
	    }
	  }
	  
	  ++cnt;
	  if (cnt % 1000000 == 0) { printf("."); fflush(stdout); }
	}
	
	printf("\nFinished!\n");
	
	bam_destroy1(b); bam_destroy1(b2);
	bam_hdr_destroy(header);
	
	sam_close(in);
	sam_close(out);
	
	return 0;
}
