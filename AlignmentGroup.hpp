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

#ifndef ALIGNMENTGROUP_H_
#define ALIGNMENTGROUP_H_

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <vector>

#include "SEQstring.hpp"
#include "QUALstring.hpp"
#include "SamParser.hpp"
#include "BamWriter.hpp"
#include "BamAlignment.hpp"

// One alignment group associates with one file. 
class AlignmentGroup {
public:
	AlignmentGroup() { 
		leftover = -1;
		s = max_size = 0;
		alignments.clear(); 
	}

	~AlignmentGroup() {
		for (int i = 0; i < max_size; ++i) delete alignments[i];
	}

	/*
		@func   clear the alignment group so that we can use it for next SAM/BAM file
	 */
	void clear() {
		leftover = -1;
		s = 0;
	}

	bool isFiltered() const { return (s > 0) && alignments[0]->isFiltered(); }

	void markAsFiltered() {
		for (int i = 0; i < s; ++i) alignments[i]->markAsFiltered();
	}

	bool read(SamParser* in);
	bool write(BamWriter* out, int choice = 0);
	
	bool isPaired() const { return (s > 0) && alignments[0]->isPaired(); }

	bool isAligned() const { return (s > 0) && (alignments[0]->isAligned() > 0); }

	int size() const { return s; }

	std::string getName(int mate = 0) const { 
		assert(s > 0);
		return alignments[0]->getName(mate);
	}
	
	int getSeqLength(int mate = 1) const { 
		assert(s > 0);
		return alignments[0]->getSeqLength(mate);
	}

	bool getSEQ(SEQstring& si, int mate = 1) {
		assert(s > 0); 
		return alignments[0]->getSEQ(si, mate);
	}

	bool getQUAL(QUALstring& qi, int mate = 1) {
		assert(s > 0);
		return alignments[0]->getQUAL(qi, mate);
	}

	BamAlignment* getAlignment(int id) { 
		assert(id >=0 && id < s);
		return alignments[id];
	}

private:
	char leftover; // if has next read. -1, initial value, standing for not called; 0, no next read; 1, has next read, which is the alignments[s]
	int s, max_size; // s, total number of alignments; max_size, max capacity
	std::vector<BamAlignment*> alignments; // pointers to BamAlignment objects

	void allocate() {
		if (s >= max_size) { 
			alignments.push_back(new BamAlignment());
			++max_size;
		}    
	}
};

inline bool AlignmentGroup::read(SamParser *in) {
	const char *cname = NULL, *name = NULL;
	BamAlignment *tmp = NULL;

	switch (leftover) {
	case 1 :
		// swap
		tmp = alignments[s]; alignments[s] = alignments[0]; alignments[0] = tmp;
		break;
	case 0 : return false;
	case -1 :
		s = 0; allocate();
		leftover = alignments[s]->read(in);
		if (leftover == 0) return false;
		break;
	default : assert(false);
	}

	cname = alignments[0]->getName();
	assert(cname[0] != 0);
	s = 1;

	while (allocate(), (leftover = alignments[s]->read(in, alignments[0]))) {
		name = alignments[s]->getName();
		if (name[0] != 0 && strcmp(cname, name)) break;
		assert(alignments[s]->isPaired() == alignments[0]->isPaired());
		++s;
	}

	return true;
}

// choice: 0, do nothing; 1, delete read sequence and qual score; 2, add read sequence and qual score.
inline bool AlignmentGroup::write(BamWriter *out, int choice) {
	assert(s > 0);
	alignments[0]->write(out);
	for (int i = 1; i < s; ++i) alignments[i]->write(out, choice, alignments[0]);
	return true;
}

#endif
