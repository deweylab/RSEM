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

#ifndef SEQSTRING_H_
#define SEQSTRING_H_

#include <cassert>
#include <string>

#include <stdint.h>
#include "htslib/sam.h"

// An iterator for BAM sequence strings
class SEQstring {
public:
	SEQstring() : seq(NULL), is_ori(true), return_current(true), len(0) {}
	
	void setUp(uint8_t *seq, int len, bool is_ori) { 
		this->seq = seq;
		this->len = len;
		this->is_ori = is_ori;
		// Default, return the original read sequence
		return_current = (is_ori ? true : false);
	}

	char getDir() { return ((is_ori && return_current) || (!is_ori && !return_current)) ? '+' : '-'; }
	
	// dir: '+', return original read sequence; '-' return reverse complement sequence
	void setDir(char dir) {
		switch(dir) {
		case '+' : return_current = (is_ori ? true : false); break;
		case '-' : return_current = (is_ori ? false : true); break;
		default: assert(false);
		}
	}

	// return current seq string
	void setCurrent() { return_current = true; }

	int getLen() const { return len; }

	char baseAt(int pos) const {
		assert(pos >= 0 && pos < len);
		return (return_current ? decode[bam_seqi(seq, pos)] : decode_r[bam_seqi(seq, len - pos - 1)]);
	}

	int baseCodeAt(int pos) const {
		assert(pos >= 0 && pos < len);
		return (return_current ? codes[bam_seqi(seq, pos)] : rcodes[bam_seqi(seq, len - pos - 1)]);
	}

	// toString will reset dir
	std::string toString(char dir = '+');

private:
	uint8_t *seq;
	bool is_ori; // if the stored seq is in read's original form
	bool return_current; // if we can return the current seq or return the reverse complement
	int len; // len, sequence length

	static const char decode[17], decode_r[17];
	static const int codes[16], rcodes[16];
};

#endif
