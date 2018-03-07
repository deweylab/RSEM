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

#ifndef MDSTRING_H_
#define MDSTRING_H_

#include <cctype>
#include <cassert>

// An iterator for MD string of a BAM alignment
class MDstring {
public:
	MDstring() : mdstr(NULL) {}

	void setUp(const char* mdstr) {
		this->mdstr = mdstr;
		counter = 0;
	}

	// return next reference base
	char next() {
		if (counter > 0) { --counter; return 0; } // MATCH
		if (isdigit(*mdstr)) {
			counter = *mdstr - '0', ++mdstr;
			while (isdigit(*mdstr)) counter = counter * 10 + (*mdstr - '0'), ++mdstr;
			if (counter > 0) { --counter; return 0; }
		}
		if (*mdstr == 0) return -1; // end of the MD string
		if (*mdstr == '^') ++mdstr; assert(isalpha(*mdstr));
		return *mdstr++;
	}
	
private:
	const char* mdstr;
	int counter;
};

#endif
