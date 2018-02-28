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

#ifndef CONVERSIONGROUP_H_
#define CONVERSIONGROUP_H_

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <map>
#include <utility>

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

struct Signature {
	BamAlignment *ba;

	Signature() : ba(NULL) {}
	Signature(BamAlignment *ba) : ba(ba) {}

	bool operator< (const Signature& o) const { return (*ba) < (*o.ba); }
};

class ConversionGroup : public AlignmentGroup {
public:
	ConversionGroup() { s = 0; align2pos.clear(); fracs.clear(); }

	// get new BamAlignment
	BamAlignment* getNewBA(const BamAlignment* o) { 
		allocate();
		alignments[s]->init_with(o);
		return alignments[s];
	}

	void pushBackBA(const BamAlignment *o, double frac, char strand1 = '.', char strand2 = '.') {
		pair.first.ba = alignments[s];
		pair.second = s;
		ret = align2pos.insert(pair);
		if (ret.second) {
			alignments[s]->completeAlignment(o, s > 0);
			if (alignments[s]->isAligned() & 1) {
				if (strand1 == '.') alignments[s]->removeTag("XS", 1);
				else alignments[s]->insertTag("XS", 'A', 1, (uint8_t*)&strand1, 1);
			}
			if (alignments[s]->isAligned() & 2) { 
				if (strand2 == '.') alignments[s]->removeTag("XS", 2);
				else alignments[s]->insertTag("XS", 'A', 1, (uint8_t*)&strand2, 2); 
			}
			if (frac >= 0.0) fracs.push_back(frac);
			++s;
		}
		else {
			if (frac > 0.0) fracs[ret.first->second] += frac;
		}
	}

	void asUnmap(const BamAlignment *o, const std::string& iv_type1, const std::string& iv_type2) {
		assert(s == 0);
		allocate();
		alignments[0]->asUnmap(o, iv_type1, iv_type2);
		++s;
	}

	void write(BamWriter* out, int choice = 0) {
		if (fracs.size() > 0) {
			assert(s == fracs.size());
			for (int i = 0; i < s; ++i) alignments[i]->setFrac(fracs[i]);
		}
		AlignmentGroup::write(out, choice);
		s = 0; align2pos.clear(); fracs.clear();
	}

private:
	std::map<Signature, int> align2pos; // a map from alignment to position in the group, used to remove duplicated alignments
	std::vector<double> fracs; // stores the ZW values
	std::pair<Signature, int> pair;
	std::pair<std::map<Signature, int>::iterator, bool> ret;
};

#endif
