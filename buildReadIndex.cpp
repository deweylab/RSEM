#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<string>
#include<fstream>
#include<iostream>

#include "utils.h"
using namespace std;

bool verbose = true;

int gap;
bool hasQ;

void buildIndex(char* readF, int gap, bool hasQ) {
	int nPos;
	READ_INT_TYPE nReads;
	bool success;
	string line;
	char idxF[STRLEN];
	char buf[sizeof(nReads) + sizeof(gap) + sizeof(nPos)];
	streampos startPos;

	sprintf(idxF, "%s.ridx", readF);

	ifstream fin(readF);
	if (!fin.is_open()) { fprintf(stderr, "Cannot open %s! It may not exist.\n", readF); exit(-1); }
	ofstream fout(idxF, ios::binary);

	startPos = fout.tellp();
	memset(buf, 0, sizeof(buf));
	fout.write((char*)buf, sizeof(buf));

	nReads = 0; nPos = 0;
	do {
		streampos pos = fin.tellg();
		success = true;

		success = (getline(fin, line));
		if (!success) continue;
		success = (getline(fin, line));
		if (!success) continue;

		if (hasQ) {
			success = (getline(fin, line));
			if (!success) continue;
			success = (getline(fin, line));
			if (!success) continue;
		}

		if (nReads % gap == 0) {
			++nPos;
			fout.write((char*)&pos, sizeof(pos));
		}
		++nReads;

		if (verbose && nReads % 1000000 == 0) { cout<< "FIN "<< nReads<< endl; }
	} while (success);

	fout.seekp(startPos);
	fout.write((char*)&nReads, sizeof(nReads));
	fout.write((char*)&gap, sizeof(gap));
	fout.write((char*)&nPos, sizeof(nPos));

	fin.close();
	fout.close();

	if (verbose) { cout<< "Build Index "<< readF<< " is Done!"<< endl; }
}

int main(int argc, char* argv[]) {
	if (argc < 5) {
		printf("Usage : rsem-build-read-index gap hasQ quiet readFile1, readFile2, ...\n");
		exit(-1);
	}

	gap = atoi(argv[1]);
	hasQ = atoi(argv[2]);
	verbose = !atoi(argv[3]);
	for (int i = 4; i < argc; i++) {
		buildIndex(argv[i], gap, hasQ);
	}

	return 0;
}
