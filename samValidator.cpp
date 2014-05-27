#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<set>

#include <stdint.h>
#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"
#include "my_assert.h"

using namespace std;

samfile_t *in;
bam1_t *b, *b2;

set<string> used;

bool isValid;

int main(int argc, char* argv[]) {
	if (argc != 2) {
		printf("Usage: rsem-sam-validator <input.sam/input.bam>\n");
		exit(-1);
	}

	string input_file(argv[1]);
	size_t pos = input_file.find_last_of('.');
	general_assert(pos != string::npos, "Input file does not have a suffix!");
	++pos;

	string suffix = input_file.substr(pos);
	size_t len = suffix.length();
	for (size_t i = 0; i < len; i++) suffix[i] = tolower(suffix[i]);

	general_assert(suffix == "sam" || suffix == "bam", "Cannot recognize input file's file type! The file suffix is neither sam nor bam.");

	in = (suffix == "sam" ? samopen(argv[1], "r", NULL) : samopen(argv[1], "rb", NULL));
	general_assert(in != 0, "Cannot open input file!");

	used.clear();
	b = bam_init1(); b2 = bam_init1();

	isValid = true;

	HIT_INT_TYPE cnt = 0;
	string cqname(""), qname;
	uint64_t creadlen = 0, readlen;
	bool cispaired = false, ispaired;

	printf("."); fflush(stdout);
	do {
		if (samread(in, b) < 0) break;
		assert(b->core.l_qseq > 0);

		qname.assign(bam1_qname(b));

		// if this is a paired-end read
		ispaired = b->core.flag & 0x0001;
		if (ispaired) {

			isValid = (samread(in, b2) >= 0) && (qname == bam1_qname(b2)) && (b2->core.flag & 0x0001);
			if (!isValid) { printf("\nOnly find one mate for paired-end read %s!\n", qname.c_str()); continue; }
			assert(b2->core.l_qseq > 0);
			isValid = !((b->core.flag & 0x0040) && (b->core.flag & 0x0080)) && !((b2->core.flag & 0x0040) & (b2->core.flag & 0x0080));
			if (!isValid) { printf("\nRead %s has more than 2 segments (e.g. tripled or more ended reads)!\n", qname.c_str()); continue; }
			isValid = !(((b->core.flag & 0x0040) && (b2->core.flag & 0x0040)) || ((b->core.flag & 0x0080) && (b2->core.flag & 0x0080)));
			if (!isValid) { printf("\nThe two mates of paired-end read %s are not adjacent!\n", qname.c_str()); continue; }


			// both mates are mapped
			if (!(b->core.flag & 0x0004) && !(b2->core.flag & 0x0004)) {
				isValid = (b->core.tid == b2->core.tid) && (b->core.pos == b2->core.mpos) && (b2->core.pos == b->core.mpos);
				if (!isValid) { printf("\nOne of paired-end read %s's alignment does not have two mates adjacent to each other! If you're running covert-sam-for-rsem now, this might mean the read contains duplicate alignments.\n", qname.c_str()); continue; }
			}

			readlen = ((b->core.flag & 0x0040) ? (uint64_t(b->core.l_qseq) << 32) + b2->core.l_qseq : (uint64_t(b2->core.l_qseq) << 32) + b->core.l_qseq);
		}
		else readlen = b->core.l_qseq;

		if (cqname != qname) {
			isValid = used.find(qname) == used.end();
			if (!isValid) { printf("\nThe alignments of read %s are not grouped together!", qname.c_str()); continue; }
			used.insert(cqname);
			cqname = qname;
			creadlen = readlen;
			cispaired = ispaired;
		}
		else {
			assert(cqname != "");
			isValid = (creadlen == readlen);
			if (!isValid) { printf("\nRead %s have different read/mate lengths!\n", qname.c_str()); continue; }
			isValid = (cispaired == ispaired);
			if (!isValid) { printf("\nRead %s is detected as both single-end read and paired-end read!\n", qname.c_str()); continue; }
		}

		++cnt;
		if (cnt % 1000000 == 0) { printf("."); fflush(stdout); }

	} while(isValid);

	bam_destroy1(b); bam_destroy1(b2);
	samclose(in);

	if (isValid) printf("\nThe input file is valid!\n");
	else printf("The input file is not valid!\n");

	return 0;
}
