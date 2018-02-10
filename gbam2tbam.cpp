#include <sys/stat.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <fstream>
#include <iostream>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "utils.h"
#include "my_assert.h"
#include "Transcripts.hpp"
#include "GenomeMap.hpp"
#include "SamParser.hpp"
#include "BamWriter.hpp"
#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

using namespace std;

bool verbose = true;

int num_threads;
htsThreadPool p = {NULL, 0};

char tiF[STRLEN], dupF[STRLEN], dupidF[STRLEN], genomeMapF[STRLEN];
Transcripts transcripts, dups;
vector<int> dup_ids;
GenomeMap genomeMap;

SamParser *parser;
BamWriter *writer;
AlignmentGroup ag;



inline bool exists(const char* file) {
	struct stat buffer;
	return (stat(file, &buffer) == 0);
}

void load_references(const char* reference_name) {
	sprintf(tiF, "%s.ti", reference_name);
	transcripts.readFrom(tiF);
	printf("Transcripts are loaded.\n");

	sprintf(dupF, "%s.dup.ti", reference_name);
	if (exists(dupF)) {
		dups.readFrom(dupF);
		sprintf(dupidF, "%s.dup.ids", reference_name);

		int dup_id;
		ifstream fin(dupidF);

		assert(fin.is_open());
		dup_ids.clear();
		while (fin>> dup_id) dup_ids.push_back(dup_id);

		printf("Duplicated transcripts are loaded.\n");
	}

	sprintf(genomeMapF, "%s.g2t.map", reference_name);
	genomeMap.readFrom(genomeMapF);
	printf("Genome map is loaded.\n");
}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		printf("Usage: gbam2tbam reference_name input.bam output.bam [-p num_threads]\n");
		exit(-1);
	}

	num_threads = 1;
	for (int i = 4; i < argc; ++i) {
		if (!strcmp(argv[i], "-p")) num_threads = atoi(argv[i + 1]);
	}

	load_references(argv[1]);

	if (num_threads > 1) assert(p.pool = hts_tpool_init(num_threads));

	parser = new SamParser(argv[2], &p);

	while (ag.read(parser)) {
		if (ag.isAligned()) {
			for (int i = 0; i < ag.size(); ++i) {
				BamAlignment* ba = ag.getAlignment(i);
				if (ba->isAligned() & 1) {
					const vector<Exon>& exon_list = genomeMap.getExonList(parser->getSeqName(ba->get_tid()), ba->getPos(1));
					cout<< ba->getName()<< "::1";
					for (int j = 0; j < (int)exon_list.size(); ++j) cout<< " ("<< exon_list[i].tid<< ", "<< exon_list[i].eid<< ")";
					cout<< endl;
				}
				if (ba->isPaired() && (ba->isAligned() & 2)) {
					const vector<Exon>& exon_list = genomeMap.getExonList(parser->getSeqName(ba->get_tid()), ba->getPos(2));
					cout<< ba->getName()<< "::2";
					for (int j = 0; j < (int)exon_list.size(); ++j) cout<< " ("<< exon_list[i].tid<< ", "<< exon_list[i].eid<< ")";
					cout<< endl;
				}
			}
		}
	}
		
	delete parser;

	if (num_threads > 1) hts_tpool_destroy(p.pool);

	return 0;
}
