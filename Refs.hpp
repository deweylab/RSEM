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

#ifndef REFS_H_
#define REFS_H_

#include <string>
#include <vector>

#include "RefSeq.hpp"

class Refs {
 public:
	Refs();
	
	void readFrom(char* inpF);
	void writeTo(char* outF);
	void writeTransListTo(char* outF);
	
	int getM() const { return M; } // get number of isoforms
	void setM(int M) { this->M = M; seqs.resize(M + 1); } 

	const RefSeq& getRef(int sid) const { return seqs[sid]; } // get a particular reference

	void setRef(int sid, const std::string& name, const std::string& seq) {
		seqs[sid].set(name, seq);
	}
	
	void appendPolyATail(int sid, int polyALen) {
		seqs[sid].appendPolyATail(polyALen);
	}
	
	void move(int from, int to) {
		assert(from >= to);
		if (from > to) seqs[to] = seqs[from];
	}

 private:
	int M; // # of isoforms, id starts from 1
	std::vector<RefSeq> seqs;  // reference sequences, starts from 1; 0 is for noise transcript
};

#endif
