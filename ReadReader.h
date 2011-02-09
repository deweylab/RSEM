#ifndef READREADER_H_
#define READREADER_H_

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<iostream>
#include<cassert>
#include<fstream>

#include "utils.h"
#include "SingleRead.h"
#include "SingleReadQ.h"
#include "PairedEndRead.h"
#include "PairedEndReadQ.h"
#include "ReadIndex.h"


template<class ReadType>
class ReadReader {
public:
	ReadReader() { s = 0; indices = NULL; arr = NULL; }
	ReadReader(int s, char readFs[][STRLEN]);
	~ReadReader();

	void setIndices(ReadIndex** indices) {
		this->indices = indices;
	}

	bool locate(long); // You should guarantee that indices exist and rid is valid, otherwise return false; If it fails, you should reset it manually!
	void reset();

	bool next(ReadType& read, int flags = 7) {
		return read.read(s, (std::istream**)arr, flags);
	}

private:
	int s; // number of files
	ReadIndex **indices;
	std::ifstream** arr;
	std::streampos *locations;
};

template<class ReadType>
ReadReader<ReadType>::ReadReader(int s, char readFs[][STRLEN]) {
	assert(s > 0);
	this->s = s;
	arr = new std::ifstream*[s];
	locations = new std::streampos[s];
	indices = NULL;
	for (int i = 0; i < s; i++) {
		arr[i] = new std::ifstream(readFs[i]);
		if (!arr[i]->is_open()) { fprintf(stderr, "Cannot open %s! It may not exist.\n", readFs[i]); exit(-1); }
		locations[i] = arr[i]->tellg();
	}
}

template<class ReadType>
ReadReader<ReadType>::~ReadReader() {
	indices = NULL;
	if (arr != NULL) {
		for (int i = 0; i < s; i++) {
			arr[i]->close();
			delete arr[i];
		}
		delete[] arr;
	}
	if (locations != NULL) {
		delete[] locations;
	}
}

template<class ReadType>
bool ReadReader<ReadType>::locate(long rid) {
	long crid = -1;
	ReadType read;

	if (indices == NULL) return false;

	//We should make sure that crid returned by each indices is the same
	for (int i = 0; i < s; i++) {
		long val = indices[i]->locate(rid, *arr[i]);
		if (i == 0) { crid = val; } else { assert(crid == val); }
	}
	assert(crid <= rid);
	while (crid < rid && read.read(s, (std::istream**)arr, 0)) ++crid;

	if (crid < rid) return false;

	std::streampos tmp[s];
	for (int i = 0; i < s; i++) { tmp[i] = arr[i]->tellg(); }

	if (!read.read(s, (std::istream**)arr, 0)) return false;

	for (int i = 0; i < s; i++) {
		locations[i] = tmp[i];
		arr[i]->seekg(locations[i]);
	}

	return true;
}

template<class ReadType>
void ReadReader<ReadType>::reset() {
	for (int i = 0; i < s; i++) {
		arr[i]->seekg(locations[i]);
	}
}

#endif /* READREADER_H_ */
