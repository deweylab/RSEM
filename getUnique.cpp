#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<vector>

#include <stdint.h>
#include "htslib/sam.h"
#include "sam_utils.h"

#include "utils.h"
#include "my_assert.h"

using namespace std;

int nThreads;
string cqname;
samFile *in, *out;
bam_hdr_t *header;
bam1_t *b;
vector<bam1_t*> arr;
bool unaligned;

void output() {
	if (unaligned || arr.size() == 0) return;
	bool isPaired = bam_is_paired(arr[0]);
	if ((isPaired && arr.size() != 2) || (!isPaired && arr.size() != 1)) return;
	for (size_t i = 0; i < arr.size(); ++i) sam_write1(out, header, arr[i]);
}

int main(int argc, char* argv[]) {
	if (argc != 4) {
		printf("Usage: rsem-get-unique number_of_threads unsorted_transcript_bam_input bam_output\n");
		exit(-1);
	}

        nThreads = atoi(argv[1]);
	in = sam_open(argv[2], "r");
	assert(in != 0);
	header = sam_hdr_read(in);
	assert(header != 0);
	out = sam_open(argv[3], "wb");
	assert(out != 0);
	sam_hdr_write(out, header);
	if (nThreads > 1) general_assert(hts_set_threads(out, nThreads) == 0, "Fail to create threads for writing the BAM file!");

	HIT_INT_TYPE cnt = 0;

	cqname = "";
	arr.clear();
	b = bam_init1();
	unaligned = false;

	while (sam_read1(in, header, b) >= 0) {
		if (cqname != bam_get_qname(b)) {
			output();
			cqname = bam_get_qname(b);
			for (size_t i = 0; i < arr.size(); ++i) bam_destroy1(arr[i]);
			arr.clear();
			unaligned = false;
		}

		unaligned = unaligned || bam_is_unmapped(b);
		arr.push_back(bam_dup1(b));

		++cnt;
		if (cnt % 1000000 == 0) { printf("."); fflush(stdout); }
	}

	if (cnt >= 1000000) printf("\n");

	output();

	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(in);
	sam_close(out);

	printf("done!\n");

	return 0;
}
