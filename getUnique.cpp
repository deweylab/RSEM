#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<vector>

#include <stdint.h>
#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"

using namespace std;

string cqname;
samfile_t *in, *out;
bam1_t *b;
vector<bam1_t*> arr;
bool unaligned;

void output() {
	if (unaligned || arr.size() == 0) return;
	bool isPaired = (arr[0]->core.flag & 0x0001);
	if ((isPaired && arr.size() != 2) || (!isPaired && arr.size() != 1)) return;
	for (size_t i = 0; i < arr.size(); i++) samwrite(out, arr[i]);
}

int main(int argc, char* argv[]) {
	if (argc != 3) {
		printf("Usage: rsem-get-unique unsorted_transcript_bam_input bam_output\n");
		exit(-1);
	}

	in = samopen(argv[1], "rb", NULL);
	assert(in != 0);
	out = samopen(argv[2], "wb", in->header);
	assert(out != 0);

	HIT_INT_TYPE cnt = 0;

	cqname = "";
	arr.clear();
	b = bam_init1();
	unaligned = false;

	while (samread(in, b) >= 0) {
		if (cqname != bam1_qname(b)) {
			output();
			cqname = bam1_qname(b);
			for (size_t i = 0; i < arr.size(); i++) bam_destroy1(arr[i]);
			arr.clear();
			unaligned = false;
		}

		unaligned = unaligned || (b->core.flag & 0x0004);
		arr.push_back(bam_dup1(b));

		++cnt;
		if (cnt % 1000000 == 0) { printf("."); fflush(stdout); }
	}

	if (cnt >= 1000000) printf("\n");

	output();

	bam_destroy1(b);
	samclose(in);
	samclose(out);

	printf("done!\n");

	return 0;
}
