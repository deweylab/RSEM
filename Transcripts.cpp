/* Copyright (c) 2017
   Bo Li (The Broad Institute of MIT and Harvard)
   libo@broadinstitute.org

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.   

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA
*/

#include <cstdio>
#include <cassert>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>

#include "utils.h"
#include "my_assert.h"

#include "Transcript.hpp"
#include "Transcripts.hpp"

void Transcripts::readFrom(const char* inpF) {
	std::string line;
	std::ifstream fin(inpF);

	general_assert(fin.is_open(), "Cannot open " + cstrtos(inpF) + "! It may not exist.");

	fin>> M>> type;
	std::getline(fin, line);
	transcripts.assign(M + 1, Transcript());
	for (int i = 1; i <= M; ++i) {
		transcripts[i].read(fin);
	}
	fin.close();
}

void Transcripts::writeTo(const char* outF) {
	std::ofstream fout(outF);
	fout<< M<< " "<< type<< std::endl;
	for (int i = 1; i <= M; ++i) {
		transcripts[i].write(fout);
	}
	fout.close();
}

void Transcripts::buildMappings(const char* imdName, int n_targets, char** target_name) {
	char file[STRLEN];

	if (imdName != NULL && target_name == NULL) {
		// Load mapping file for outputting BAM
		sprintf(file, "%s.mappings", imdName);
		FILE *fi = fopen(file, "r");
		assert(fi != NULL);

		e2i = new int[M];
		i2e = new int[M + 1];
		// the first entry, 0 is omitted here
		for (int i = 0; i < M; i++) assert(fscanf(fi, "%d", &e2i[i]) == 1);
		for (int i = 1; i <= M; i++) assert(fscanf(fi, "%d", &i2e[i]) == 1);
		fclose(fi);

		return;
	}

	std::map<std::string, int> dict;
	std::map<std::string, int>::iterator iter;

	general_assert(n_targets > 0, "The SAM/BAM file declares less than one reference sequence!");
	general_assert(n_targets <=  M, "The SAM/BAM file declares more reference sequences (" + itos(n_targets) + ") than RSEM knows (" + itos(M) + ")!");
	if (n_targets < M) printf("Warning: The SAM/BAM file declares less reference sequences (%d) than RSEM knows (%d)!\n", n_targets, M);

	dict.clear();
	for (int i = 1; i <= M; ++i) {
		const std::string& tid = isAlleleSpecific() ? transcripts[i].getSeqName() : transcripts[i].getTranscriptID();
		iter = dict.find(tid);
		general_assert(iter == dict.end(), "RSEM's indices might be corrupted, " + tid + " appears more than once!");
		dict[tid] = i;
	}

	e2i = new int[M];
	for (int i = 0; i < M; ++i) e2i[i] = -1;
	i2e = new int[M + 1];
	for (int i = 0; i <= M; ++i) i2e[i] = -1;

	for (int i = 0; i < n_targets; ++i) {
		iter = dict.find(std::string(target_name[i]));
		general_assert(iter != dict.end(), "RSEM can not recognize reference sequence name " + cstrtos(target_name[i]) + "!");
		general_assert(iter->second > 0, "Reference sequence name " + cstrtos(target_name[i]) + " appears more than once in the SAM/BAM file!");
		e2i[i] = iter->second;
		i2e[iter->second] = i;
		iter->second = -1;
	}

	if (imdName != NULL) {
		sprintf(file, "%s.mappings", imdName);
		FILE *fo = fopen(file, "w");
		for (int i = 0; i < M; ++i) fprintf(fo, "%d%c", e2i[i], (i == M - 1 ? '\n' : '\t'));
		for (int i = 1; i <= M; ++i) fprintf(fo, "%d%c", i2e[i], (i == M ? '\n' : '\t'));
		fclose(fo);

		sprintf(file, "%s.omit", imdName);
		fo = fopen(file, "w");

		for (int i = 1; i <= M; ++i) 
		if (i2e[i] < 0) fprintf(fo, "%d\n", i);
		fclose(fo);
	}
}
