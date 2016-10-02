#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <string>
#include <set>
#include <algorithm>

#include <stdint.h>
#include "htslib/sam.h"
#include "sam_utils.h"

#include "utils.h"
#include "my_assert.h"

using namespace std;

samFile *in;
bam_hdr_t *header;
bam1_t *b, *b2;

set<string> used;

bool isValid;

bool check_read(bam1_t *b, bam_hdr_t *header) {
	uint32_t* cigar = bam_get_cigar(b);
	for (int i = 0; i < b->core.n_cigar; ++i) {
		char op = bam_cigar_opchr(*cigar);
		if (op == 'N') {
			printf("\nSkipped region is detected (cigar N) for read %s!\nTo use RSEM, please align your reads to a set of transcript sequences instead of a genome.\n", bam_get_qname(b));
			return false;
		}
		else if (op == 'I' || op == 'D') {
			printf("\nIndel alignment is detected (cigar %c) for read %s!\nRSEM currently does not support indel alignments.\n", op, bam_get_qname(b));
			return false;
		}
		else if (op == 'S' || op == 'H' || op == 'P') {
			printf("\nClipping or padding is detected (cigar %c) for read %s!\nRSEM currently doest not support clipping or padding.\n", op, bam_get_qname(b));
			return false;
		}
		++cigar;
	}

	if (b->core.pos < 0 || bam_endpos(b) > header->target_len[b->core.tid]) {
		printf("\n");
		if (bam_is_paired(b)) {
			printf("Mate %d of paired-end read %s", (bam_is_read1(b) ? 1 : 2), bam_get_qname(b));
		}
		else {
			printf("Read %s", bam_get_qname(b));
		}
		printf(" aligns to [%d, %d) of transcript %s, which exceeds the transcript's boundary [0, %d)!\n",
				b->core.pos, bam_endpos(b), header->target_name[b->core.tid], header->target_len[b->core.tid]);
		return false;
	}

	return true;
}

int main(int argc, char* argv[]) {
	if (argc != 2) {
		printf("Usage: rsem-sam-validator <input.sam/input.bam/input.cram>\n");
		exit(-1);
	}
	
	in = sam_open(argv[1], "r");
	general_assert(in != 0, "Cannot open input file!");
	header = sam_hdr_read(in);
	general_assert(header != 0, "Cannot load SAM header!");
	used.clear();
	b = bam_init1(); b2 = bam_init1();
	
	isValid = true;
	
	HIT_INT_TYPE cnt = 0;
	string cqname(""), qname;
	int creadlen = 0, readlen, creadlen2, readlen2;
	char ispaired = -1;
	
	printf("."); fflush(stdout);
	do {
		int ret = sam_read1(in, header, b);
		if (ret == -1) break;
		else if (ret < 0) { isValid = false; continue; }
		assert(b->core.l_qseq > 0);
		
		qname = bam_get_canonical_name(b);
		
		// if this is a paired-end read
		if (ispaired == -1) ispaired = bam_is_paired(b);
		else {
			isValid = (ispaired == bam_is_paired(b));
			if (!isValid) {
				printf("\nWe detected both single-end and paired-end reads in the data!\nRSEM currently does not support a mixture of single-end/paired-end reads.\n");
				continue;
			}
		}

		if (ispaired) {
			isValid = (sam_read1(in, header, b2) >= 0) && (qname == bam_get_canonical_name(b2)) && bam_is_paired(b2);
			if (!isValid) { 
				printf("\nOnly find one mate for paired-end read %s!\nPlease make sure that the two mates of a paired-end read are adjacent to each other.\n", bam_get_qname(b)); 
				continue; 
			}

			assert(b2->core.l_qseq > 0);

			isValid = (bam_is_read1(b) && bam_is_read2(b2)) || (bam_is_read1(b2) && bam_is_read2(b));
			if (!isValid) { 
				printf("\nThe two mates of paired-end read %s are marked as both mate1 or both mate2!\n", bam_get_qname(b)); 
				continue; 
			}

			int value = int(bam_is_mapped(b)) + int(bam_is_mapped(b2));
			isValid = (value != 1);
			if (!isValid) {
				printf("\nPaired-end read %s has an alignment with only one mate aligned!\n", bam_get_qname(b));
				printf("Currently RSEM does not handle mixed alignments for paired-end reads.\n");
				continue;
			}

			if (!bam_is_read1(b)) { bam1_t *tmp = b; b = b2; b2 = tmp; }

			if (value == 2) {
				isValid = (b->core.tid == b2->core.tid);
				if (!isValid) {
					printf("\nPaired-end read %s has a discordant alignment (two mates aligned to different reference sequences)!\n", bam_get_qname(b));
					printf("Mate 1 aligns to %s and mate 2 aligns to %s\n", header->target_name[b->core.tid], header->target_name[b2->core.tid]);
					printf("Currently RSEM does not handle discordant alignments.\n");
					continue;
				}

				int strandedness = (int(bam_is_rev(b)) << 1) + int(bam_is_rev(b2));
				isValid = (strandedness == 1 || strandedness == 2);
				if (!isValid) {
					printf("\nPaired-end read %s has an alignment in which two mates aligned to the same strand!\n", bam_get_qname(b));
					printf("Its two mates aligned to %s in %s direction.\n", header->target_name[b->core.tid], (strandedness == 0 ? "forward" : "reverse"));
					continue;
				}

				bam1_t *tb = (b->core.pos < b2->core.pos ? b : b2);
				isValid = tb->core.pos >= 0 && tb->core.pos + abs(tb->core.isize) <= header->target_len[tb->core.tid];
				if (!isValid) {
					printf("\nPaired-end read %s aligns to [%d, %d) of transcript %s, which exceeds the transcript's boundary [0, %d)!\n", 
						bam_get_qname(b), tb->core.pos, tb->core.pos + abs(tb->core.isize), header->target_name[tb->core.tid], header->target_len[tb->core.tid]);
					continue;
				}
				isValid = check_read(b, header);
				if (!isValid) continue;
				isValid = check_read(b2, header);
				if (!isValid) continue;
			}
			
			readlen = b->core.l_qseq;
			readlen2 = b2->core.l_qseq;
		}
		else {
			if (bam_is_mapped(b)) {
				isValid = check_read(b, header);
				if (!isValid) continue;
			}
			readlen = b->core.l_qseq;
		}

		if (cqname != qname) {
			isValid = used.find(qname) == used.end();
			if (!isValid) { printf("\nThe alignments of read %s are not grouped together!\n", qname.c_str()); continue; }
			used.insert(cqname);
			cqname = qname;
			creadlen = readlen;
			if (ispaired) creadlen2 = readlen2;			
		}
		else {
			assert(cqname != "");
			isValid = (creadlen == readlen && (!ispaired || creadlen2 == readlen2));
			if (!isValid) { printf("\nRead %s have alignments showing different read/mate lengths!\n", qname.c_str()); continue; }
		}
		
		++cnt;
		if (cnt % 1000000 == 0) { printf("."); fflush(stdout); }
		
	} while(isValid);
	
	bam_destroy1(b); bam_destroy1(b2);
	bam_hdr_destroy(header);
	sam_close(in);
	
	if (isValid) printf("\nThe input file is valid!\n");
	else printf("The input file is not valid!\n");
	
	return 0;
}
