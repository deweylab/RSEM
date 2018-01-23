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

#ifndef QUALSTRING_H_
#define QUALSTRING_H_

#include <cassert>
#include <string>
#include <sstream>

#include <stdint.h>

// An iterator for BAM quality score strings
class QUALstring {
public:
	QUALstring() : qual(NULL), is_ori(true), return_current(true), len(0) {}
	
	void setUp(uint8_t *qual, int len, bool is_ori) { 
		this->qual = qual;
		this->len = len;
		this->is_ori = is_ori;
		// Default, return the original quality score sequence
		return_current = (is_ori ? true : false);
	}

	char getDir() { return ((is_ori && return_current) || (!is_ori && !return_current)) ? '+' : '-'; }
	
	// '+' returns the original qual string; '-' returns the reverse string
	void setDir(char dir) { 
		assert(dir == '+' || dir == '-');
		switch(dir) {
		case '+' : return_current = (is_ori ? true : false); break;
		case '-' : return_current = (is_ori ? false : true); break;
		default: assert(false);
		}
	}
	
	int getLen() const { return len; }
	
	// 33 is already deducted
	int qualAt(int pos) const {
		assert(pos >= 0 && pos < len);
		return qual[return_current ? pos : len - pos - 1];
	}

	// default is the original quality score string
	std::string toString(char dir = '+') {
		setDir(dir);
		std::ostringstream strout;
		for (int i = 0; i < len; ++i) strout<< char(qualAt(i) + 33);
		return strout.str();
	}

private:
	uint8_t *qual;
	bool is_ori; // if the stored qual is in read's original form
	bool return_current; // if we can return the current seq or return the reverse
	int len; // len, quality score sequence length
};

#endif
