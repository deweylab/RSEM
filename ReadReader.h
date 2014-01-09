#ifndef READREADER_H_
#define READREADER_H_

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<iostream>
#include<cassert>
#include<fstream>
#include<vector>

#include "utils.h"
#include "SingleRead.h"
#include "SingleReadQ.h"
#include "PairedEndRead.h"
#include "PairedEndReadQ.h"
#include "ReadIndex.h"


template<class ReadType>
class ReadReader {
public:
	ReadReader() { s = 0; indices = NULL; arr = NULL; hasPolyA = false; seedLen = -1; }
	ReadReader(int s, char readFs[][STRLEN], bool hasPolyA = false, int seedLen = -1);
	~ReadReader();

	void setIndices(ReadIndex** indices) {
		this->indices = indices;
	}

	bool locate(READ_INT_TYPE); // You should guarantee that indices exist and rid is valid, otherwise return false; If it fails, you should reset it manually!
	void reset();

	bool next(ReadType& read, int flags = 7) {
		bool success = read.read(s, (std::istream**)arr, flags);
		if (success && seedLen > 0) { read.calc_lq(hasPolyA, seedLen); }
		return success;
	}

private:
	int s; // number of files
	ReadIndex **indices;
	std::ifstream** arr;
	std::streampos *locations;

	bool hasPolyA;
	int seedLen;
};

template<class ReadType>
ReadReader<ReadType>::ReadReader(int s, char readFs[][STRLEN], bool hasPolyA, int seedLen) {
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
	this->hasPolyA = hasPolyA;
	this->seedLen = seedLen;
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
bool ReadReader<ReadType>::locate(READ_INT_TYPE rid) {
	READ_INT_TYPE crid = -1;
	ReadType read;

	if (indices == NULL) return false;

	//We should make sure that crid returned by each indices is the same
	for (int i = 0; i < s; i++) {
		READ_INT_TYPE val = indices[i]->locate(rid, *arr[i]);
		if (i == 0) { crid = val; } else { assert(crid == val); }
	}
	assert(crid <= rid);
	while (crid < rid && read.read(s, (std::istream**)arr, 0)) ++crid;

	if (crid < rid) return false;

	std::vector<std::streampos> tmp(s);
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
