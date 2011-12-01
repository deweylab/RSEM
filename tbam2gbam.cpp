#include<cstdio>
#include<cstring>
#include<cstdlib>

#include "utils.h"
#include "Transcripts.h"
#include "BamConverter.h"

using namespace std;

char tiF[STRLEN], chr_list[STRLEN];
Transcripts transcripts;

int main(int argc, char* argv[]) {
	if (argc != 4) {
		printf("Usage: rsem-tbam2gbam reference_name unsorted_transcript_bam_input genome_bam_output\n");
		exit(-1);
	}

	sprintf(tiF, "%s.ti", argv[1]);
	sprintf(chr_list, "%s.chrlist", argv[1]);
	transcripts.readFrom(tiF);

	printf("Start converting:\n");
	BamConverter bc(argv[2], argv[3], chr_list, transcripts);
	bc.process();
	printf("Genome bam file is generated!\n");

	return 0;
}

