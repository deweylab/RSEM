#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>

#include "sam/bam.h"
#include "sam/sam.h"

using namespace std;

samfile_t *bam_in;
bam1_t *b;

int cur_tid; //current tid;
float *wig_arr; // wiggle array
FILE *fo;

void generateWiggle(int tid) {
	int chr_len = bam_in->header->target_len[tid];
	char *chr_name = bam_in->header->target_name[tid];
	int sp, ep;

	sp = ep = -1;
	for (int i = 0; i < chr_len; i++) {
		if (wig_arr[i] > 0) {
			ep = i;
		}
		else {
			if (sp < ep) {
				++sp;
				fprintf(fo, "fixedStep chrom=%s start=%d step=1\n", chr_name, sp + 1);
				for (int j = sp; j <= ep; j++) fprintf(fo, "%.7g\n", wig_arr[j]);
			}
			sp = i;
		}
	}
	if (sp < ep) {
		++sp;
		fprintf(fo, "fixedStep chrom=%s start=%d step=1\n", chr_name, sp + 1);
		for (int j = sp; j <= ep; j++) fprintf(fo, "%.7g\n", wig_arr[j]);
	}
}

int main(int argc, char* argv[]) {
	int cnt = 0;

	if (argc != 4) {
		printf("Usage : rsem-bam2wig sorted_bam_input wig_output wiggle_name\n");
		exit(-1);
	}

	bam_in = samopen(argv[1], "rb", NULL);
	if (bam_in == 0) { fprintf(stderr, "Cannot open %s!\n", argv[1]); exit(-1); }
	//assert(bam_in != 0);
	b = bam_init1();

	fo = fopen(argv[2], "w");
	fprintf(fo, "track type=wiggle_0 name=\"%s\" description=\"%s\" visibility=full\n", argv[3], argv[3]);

	cur_tid = -1;
	wig_arr = NULL;
	while (samread(bam_in, b) >= 0) {
		if (b->core.tid != cur_tid) {
			if (cur_tid >= 0) generateWiggle(cur_tid);
			cur_tid = b->core.tid;
			size_t len = sizeof(float) * bam_in->header->target_len[cur_tid];
			wig_arr = (float*)realloc(wig_arr, len);
			memset(wig_arr, 0, len);
		}

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
				for (int j = 0; j < op_len; j++, ++pos) wig_arr[pos] += w;
				break;
			case BAM_CREF_SKIP : pos += op_len; break;
			default : assert(false);
			}
		}

		++cnt;
		if (cnt % 1000000 == 0) printf("%d FIN\n", cnt);
	}
	if (cur_tid >= 0) generateWiggle(cur_tid);
	free(wig_arr);

	samclose(bam_in);
	bam_destroy1(b);

	fclose(fo);

	return 0;
}
