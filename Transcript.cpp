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

#include <cassert>
#include <string>
#include <vector>
#include <fstream>

#include "utils.h"
#include "my_assert.h"
#include "Transcript.hpp"

//gseq : genomic sequence
void Transcript::extractSeq(const std::string& gseq, std::string& seq) const {
	seq = "";
	int s = structure.size();

	general_assert(structure[0].start >= 1 && structure[s - 1].end <= (int)gseq.length(), "Transcript " + transcript_id + " is out of chromosome " + seqname + "'s boundary!");
 
	switch(strand) {
	case '+':
		for (int i = 0; i < s; ++i) {
			seq += gseq.substr(structure[i].start - 1, structure[i].end - structure[i].start + 1); // gseq starts from 0!
		}
		break;
	case '-':
		for (int i = s - 1; i >= 0; --i) {
			for (int j = structure[i].end; j >= structure[i].start; --j) {
				seq += base2rbase[gseq[j - 1]];
			}
		}
		break;
	default: assert(false);
	}

	assert(seq.length() > 0);
}

void Transcript::read(std::ifstream& fin) {
	int s;
	std::string tmp;
	std::istringstream strin;

	getline(fin, tmp);
	strin.str(tmp);
	getline(strin, transcript_id, '\t');
	getline(strin, transcript_name);

	getline(fin, tmp);
	strin.clear(); strin.str(tmp);
	getline(strin, gene_id, '\t');
	getline(strin, gene_name);

	std::getline(fin, seqname);

	fin>> tmp>> length;
	assert(tmp.length() == 1 && (tmp[0] == '+' || tmp[0] == '-'));
	strand = tmp[0];
	structure.clear();
	fin>> s;
	for (int i = 0; i < s; ++i) {
		int start, end;
		fin>> start>> end;
		structure.push_back(Interval(start, end));
	}
	std::getline(fin, tmp); //get the end of this line
	std::getline(fin, left);
}

void Transcript::write(std::ofstream& fout) {
	int s = structure.size();

	fout<< transcript_id;
	if (transcript_name != "") fout<< '\t'<< transcript_name;
	fout<< std::endl;

	fout<< gene_id;
	if (gene_name != "") fout<< '\t'<< gene_name;
	fout<< std::endl;

	fout<< seqname<< std::endl;
	fout<< strand<<" "<< length<< std::endl;
	fout<< s;
	for (int i = 0; i < s; ++i) fout<< " "<< structure[i].start<< " "<< structure[i].end;
	fout<< std::endl;
	fout<< left<< std::endl;
}
