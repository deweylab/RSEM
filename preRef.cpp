#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cctype>
#include<string>
#include<fstream>
#include<sstream>
#include<cassert>

#include "utils.h"
#include "Refs.h"
#include "PolyARules.h"
#include "RefSeqPolicy.h"
#include "AlignerRefSeqPolicy.h"

using namespace std;

bool verbose = true;

int M;

RefSeqPolicy refp;
AlignerRefSeqPolicy aligner_refp;
PolyARules rules;
Refs refs;

ofstream fout;
char refF[STRLEN], idxF[STRLEN], n2g_idxF[STRLEN];

int polyAChoice, polyALen;
char exceptionF[STRLEN];
bool quiet; // verbose = !quiet;

// always generate references for aligners, default convert all N into G
int main(int argc, char* argv[]) {

	if (argc < 4) {
		printf("USAGE : rsem-preref refFastaF polyAChoice refName [-l polyALen] [-f exceptionF] [-q]\n\n");
		printf("  refFastaF: a FASTA format file contains all reference transcripts\n");
		printf("  polyAChoice: choice for polyA tail padding.It is a number from {0,1,2}\n");
		printf("    0: pad polyA tail\n");
		printf("    1: do not pad polyA tail at all\n");
		printf("    2: pad polyA tail for all references but those in exceptionF\n");
		printf("  -l: polyALen: specify the length of polyA tail you want to pad. Default is 100\n");
		printf("  -f: exceptionF: file contains a list of exception reference ids. IDs starts from 1. Must set if polyAChoice = 2\n");
		printf("  -q: quiet\n");
		exit(-1);
	}

	polyAChoice = atoi(argv[2]);

	polyALen = 125;
	quiet = false;
	memset(exceptionF, 0, sizeof(exceptionF));

	for (int i = 4; i < argc; i++) {
		if (!strcmp(argv[i], "-l")) { polyALen = atoi(argv[i + 1]); }
		if (!strcmp(argv[i], "-f")) { strcpy(exceptionF, argv[i + 1]); }
		if (!strcmp(argv[i], "-q")) { quiet = true; }
	}

	verbose = !quiet;

	//make references
	rules = PolyARules(polyAChoice, polyALen, exceptionF);
	refs.makeRefs(argv[1], refp, rules);
	M = refs.getM();

	//save references
	sprintf(refF, "%s.seq", argv[3]);
	refs.saveRefs(refF);

	sprintf(idxF, "%s.idx.fa", argv[3]);
	fout.open(idxF);
	for (int i = 1; i <= M; i++) {
		fout<< ">"<< refs.getRef(i).getName()<< endl<< refs.getRef(i).getSeq()<< endl;
	}
	fout.close();
	if (verbose) printf("%s is generated!\n", idxF);

	sprintf(n2g_idxF, "%s.n2g.idx.fa", argv[3]);
	fout.open(n2g_idxF);
	for (int i = 1; i <= M; i++) {
	  fout<< ">"<< refs.getRef(i).getName()<< endl<< aligner_refp.convert(refs.getRef(i).getSeq())<< endl;
	}
	fout.close();
	if (verbose) printf("%s is generated!\n", n2g_idxF);

	return 0;
}
