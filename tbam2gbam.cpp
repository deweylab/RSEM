#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>

#include "utils.h"
#include "Transcripts.h"
#include "BamConverter.h"

using namespace std;

int nThreads;
char tiF[STRLEN], chr_list[STRLEN];
Transcripts transcripts;

int main(int argc, char* argv[]) {
	if (argc != 4 && argc != 6) {
		printf("Usage: rsem-tbam2gbam reference_name unsorted_transcript_bam_input genome_bam_output [-p number_of_threads]\n");
		exit(-1);
	}

        nThreads = 1; // default is 1
        if (argc == 6) { assert(strcmp(argv[4], "-p") == 0); nThreads = atoi(argv[5]); }

	sprintf(tiF, "%s.ti", argv[1]);
	sprintf(chr_list, "%s.chrlist", argv[1]);
	transcripts.readFrom(tiF);

	printf("Start converting:\n");
	BamConverter bc(argv[2], argv[3], chr_list, transcripts, nThreads, assemble_command(argc, argv));
	bc.process();
	printf("Genome bam file is generated!\n");

	return 0;
}

