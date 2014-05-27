#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<vector>
#include<algorithm>

#include <stdint.h>
#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"
#include "my_assert.h"

using namespace std;

samfile_t *in, *out;
bam1_t *b, *b2;
vector<bam1_t*> arr_both, arr_partial_1, arr_partial_2, arr_partial_unknown;

inline void add_to_appropriate_arr(bam1_t *b) {
	if (!(b->core.flag & 0x0004) && (b->core.flag & 0x0002)) {
		arr_both.push_back(bam_dup1(b)); return;
	}

	if (b->core.flag & 0x0040) arr_partial_1.push_back(bam_dup1(b));
	else if (b->core.flag & 0x0080) arr_partial_2.push_back(bam_dup1(b));
	else arr_partial_unknown.push_back(bam_dup1(b));
}

char get_pattern_code(uint32_t flag) {
  if (flag & 0x0040) return (flag & 0x0010 ? 1 : 0);
  else return (flag & 0x0010 ? 0 : 1);
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
	if (argc != 3) {
		printf("Usage: rsem-scan-for-paired-end-reads input.sam output.bam\n");
		exit(-1);
	}

	in = samopen(argv[1], "r", NULL);
	general_assert(in != 0, "Cannot open " + cstrtos(argv[1]) + " !");
	out = samopen(argv[2], "wb", in->header);
	general_assert(out != 0, "Cannot open " + cstrtos(argv[2]) + " !");

	b = bam_init1(); b2 = bam_init1();

	string qname;
	bool go_on = (samread(in, b) >= 0);
	bool isPaired;
	HIT_INT_TYPE cnt = 0;

	printf("."); fflush(stdout);

	while (go_on) {
		qname.assign(bam1_qname(b));
		isPaired = (b->core.flag & 0x0001);

		if (isPaired) {
			add_to_appropriate_arr(b);
			while ((go_on = (samread(in, b) >= 0)) && (qname == bam1_qname(b))) {
				general_assert_1(b->core.flag & 0x0001, "Read " + qname + " is detected as both single-end and paired-end read!");
				add_to_appropriate_arr(b);
			}

			general_assert_1(arr_both.size() % 2 == 0, "Number of first and second mates in read " + qname + "'s full alignments (both mates are aligned) are not matched!");
			general_assert_1((arr_partial_1.size() + arr_partial_2.size() + arr_partial_unknown.size()) % 2 == 0, "Number of first and second mates in read " + qname + "'s partial alignments (at most one mate is aligned) are not matched!");

			if (!arr_both.empty()) {
				sort(arr_both.begin(), arr_both.end(), less_than);
				for (size_t i = 0; i < arr_both.size(); i++) { samwrite(out, arr_both[i]); bam_destroy1(arr_both[i]); }
				arr_both.clear();
			}

			while (!arr_partial_1.empty() || !arr_partial_2.empty()) {
				if (!arr_partial_1.empty() && !arr_partial_2.empty()) {
					samwrite(out, arr_partial_1.back()); bam_destroy1(arr_partial_1.back()); arr_partial_1.pop_back();
					samwrite(out, arr_partial_2.back()); bam_destroy1(arr_partial_2.back()); arr_partial_2.pop_back();
				}
				else if (!arr_partial_1.empty()) {
					samwrite(out, arr_partial_1.back()); bam_destroy1(arr_partial_1.back()); arr_partial_1.pop_back();
					samwrite(out, arr_partial_unknown.back()); bam_destroy1(arr_partial_unknown.back()); arr_partial_unknown.pop_back();
				}
				else {
					samwrite(out, arr_partial_2.back()); bam_destroy1(arr_partial_2.back()); arr_partial_2.pop_back();
					samwrite(out, arr_partial_unknown.back()); bam_destroy1(arr_partial_unknown.back()); arr_partial_unknown.pop_back();
				}
			}

			while (!arr_partial_unknown.empty()) {
				samwrite(out, arr_partial_unknown.back()); bam_destroy1(arr_partial_unknown.back()); arr_partial_unknown.pop_back();
			}
		}
		else {
			samwrite(out, b);
			while ((go_on = (samread(in, b) >= 0)) && (qname == bam1_qname(b))) {
				samwrite(out, b);
			}
		}

		++cnt;
		if (cnt % 1000000 == 0) { printf("."); fflush(stdout); }
	}

	printf("\nFinished!\n");

	bam_destroy1(b); bam_destroy1(b2);

	samclose(in);
	samclose(out);

	return 0;
}
