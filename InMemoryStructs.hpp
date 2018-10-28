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

#ifndef INMEMORYSTRUCTS_H_
#define INMEMORYSTRUCTS_H_

#include <vector>

#include "utils.h"

// In memory alignment
struct InMemAlign {
	int tid; // transcript id
	double conprb, frac; // conprb, the conditional probability of generating the read based on read_model; frac the expected weight

	InMemAlign() : tid(0), conprb(0.0), frac(0.0) {}

	InMemAlign(int tid, double conprb, double frac = 0.0) : tid(tid), conprb(conprb), frac(frac) {}
};

// In memory alignment group
struct InMemAlignG {
	int size; 
	double noise_conprb, noise_frac;

	InMemAlignG() : size(0), noise_conprb(0.0) {}
};

// Store in memory information for all alignments of a thread
struct InMemChunk {
	READ_INT_TYPE pos, nreads;
	InMemAlign *aligns, *pointer;
	InMemAlignG *reads;

	InMemChunk(READ_INT_TYPE nreads, HIT_INT_TYPE nlines) {
		this->nreads = nreads;
		reads = new InMemAlignG[nreads];
		--reads; // because reads[1] is the start position
		aligns = new InMemAlign[nlines];

		reset();
	}

	/*
		@func   reset InMemAlignG.aligns to NULL
	 */
	void reset() {
		pos = 0;
		pointer = aligns;
	}

	/*
		@param   aread     In memory read group information
		@param   alignArr  In memory alignment pointer
		@return  true if has more reads, false otherwise
	 */
	bool next(InMemAlignG*& aread, InMemAlign*& alignArr) {
		if (pos >= nreads) return false;
		if (pos > 0) pointer += reads[pos].size;
		++pos;

		aread = reads + pos;
		alignArr = pointer;

		return true;
	}

	~InMemChunk() {
		++reads;
		delete[] reads;
		delete[] aligns;
	}  
};

#endif
