#ifndef READINDEX_H_
#define READINDEX_H_

#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<fstream>

#include "utils.h"

struct ReadIndex {
	READ_INT_TYPE nReads;
	int gap, nPos;
	std::streampos *index;

	ReadIndex () {
		nReads = 0; gap = nPos = 0;
		index = NULL;
	}

	ReadIndex(const char* readF) {
		char indexF[STRLEN];
		std::ifstream fin;

		sprintf(indexF, "%s.ridx", readF);
		fin.open(indexF, std::ios::binary);
		if (!fin.is_open()) { fprintf(stderr, "Cannot open %s! It may not exist.\n", indexF); exit(-1); }

		nReads = 0; gap = nPos = 0;
		index = NULL;
		if (fin.is_open()) {
			fin.read((char*)&nReads, sizeof(nReads));
			fin.read((char*)&gap, sizeof(gap));
			fin.read((char*)&nPos, sizeof(nPos));
			index = new std::streampos[nPos];
			for (int i = 0; i < nPos; i++) {
				fin.read((char*)&index[i], sizeof(std::streampos));
			}
		}
	}

	~ReadIndex() {
		nReads = 0; gap = nPos = 0;
		if (index != NULL) delete[] index;
	}

	//rid  0-based , return crid : current seeked rid
	READ_INT_TYPE locate(READ_INT_TYPE rid, std::ifstream& out) {
		if (index == NULL) {
			out.seekg(0, std::ios::beg);
			return 0;
		}
		assert(rid >= 0 && rid < nReads);
		out.seekg(index[rid / gap]);
		return (rid / gap) * gap;
	}
};

#endif /* READINDEX_H_ */
