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

#ifndef CIGARSTRING_H_
#define CIGARSTRING_H_

#include <cstdint>
#include <cassert>
#include <string>
#include <sstream>

#include "htslib/sam.h"

// An iterator for BAM CIGAR strings
class CIGARstring {
public:
	CIGARstring() : cigar(NULL), is_ori(true), return_current(true), len(0) {}
	
	void setUp(const uint32_t *cigar, int len, bool is_ori) { 
		this->cigar = cigar;
		this->len = len;
		this->is_ori = is_ori;
		// Default, return the CIGAR string in original read sequence's order
		return_current = (is_ori ? true : false);
	}

	char getDir() { return ((is_ori && return_current) || (!is_ori && !return_current)) ? '+' : '-'; }
	
	// '+' returns CIGAR string in original read sequence's order; '-' returns the reverse CIGAR
	void setDir(char dir) { 
		assert(dir == '+' || dir == '-');
		switch(dir) {
		case '+' : return_current = (is_ori ? true : false); break;
		case '-' : return_current = (is_ori ? false : true); break;
		default: assert(false);
		}
	}

	// return current cigar string
	void setCurrent() { return_current = true; }

	int getLen() const { return len; }

	uint32_t valueAt(int pos) const {
		assert(pos >= 0 && pos < len);
		return return_current ? cigar[pos] : cigar[len - pos - 1];
	}

	int opAt(int pos) const { return bam_cigar_op(valueAt(pos)); }

	char opchrAt(int pos) const { return bam_cigar_opchr(valueAt(pos)); }

	int oplenAt(int pos) const { return bam_cigar_oplen(valueAt(pos)); }

	// 0: consume nothing; 1: query; 2: reference; 3: both
	int optypeAt(int pos) const { return bam_cigar_type(opAt(pos)); }
	
	// toString will reset dir
	std::string toString(char dir = '+') {
		setDir(dir);
		std::ostringstream strout;
		for (int i = 0; i < len; ++i) strout<< oplenAt(i)<< opchrAt(i); 
		return strout.str();
	}

private:
	const uint32_t *cigar;
	bool is_ori; // if the stored cigar is in original read sequence's order
	bool return_current; // if we can return the current cigar or return the reverse
	int len; // len, cigar string length
};

#endif
