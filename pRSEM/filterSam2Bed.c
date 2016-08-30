/*
 *  pliu 20150621
 *
 *  filter Sam file by flag 1548 and output alignment in Bed format
 *
 *  this code is modified from sam/examples/bam2bed.c
 *
*/

#include <stdio.h>
#include "sam.h"
//#include "../samtools-1.3/htslib-1.3/htslib/sam.h"

static int fetch_func(const bam1_t *b, void *data) {
	samfile_t *fp = (samfile_t*)data;
	uint32_t *cigar = bam1_cigar(b);
	const bam1_core_t *c = &b->core;
	int i, l;
	if (b->core.tid < 0) return 0;
	if ( (b->core.flag & 0x4) || (b->core.flag & 0x8) || (b->core.flag & 0x200) ||
		   (b->core.flag & 0x400) ) return 0;
	for (i = l = 0; i < c->n_cigar; ++i) {
		int op = cigar[i]&0xf;
		if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
			l += cigar[i]>>4;
	}
	printf("%s\t%d\t%d\tN\t%d\t%c\n", fp->header->target_name[c->tid],
		   c->pos, c->pos + l, c->qual, (c->flag&BAM_FREVERSE)? '-' : '+');
	return 0;
}


int main(int argc, char *argv[]) {
	samfile_t *fp;
	if (argc != 2) {
		fprintf(stderr, "\nUsage: filterSam2Bed <in.sam>\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Filter SAM file by flag 1548 (0x4, 0x8, 0x200, 0x400)\n"); 
		fprintf(stderr, "and write to STDOUT in BED format\n");
		fprintf(stderr, "<in.sam>: input SAM file name, '-' for STDIN\n\n");
		return 1;
	}
	if ((fp = samopen(argv[1], "r", 0)) == 0) {
		fprintf(stderr, "filterSam2Bed: Fail to open SAM file %s\n", argv[1]);
		return 1;
	}
	bam1_t *b = bam_init1();
	while (samread(fp, b) >= 0) fetch_func(b, fp);
	bam_destroy1(b);
	samclose(fp);
	return 0;
}
