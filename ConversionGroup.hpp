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
#include <set>

#include "BamAlignment.hpp"
#include "AlignmentGroup.hpp"

struct Signature {
	BamAlignment *ba;

	Signature(BamAlignment *ba) : ba(ba) {}

	bool operator< (const Signature& o) const { return (*ba) < (*o.ba); }
};

class ConversionGroup : public AlignmentGroup {
public:
	ConversionGroup() { has_seen.clear(); }

	void clear() {
		s = 0;
		has_seen.clear();
	}

	// get new BamAlignment
	BamAlignment* getNewBA(const BamAlignment* o) { 
		allocate();
		alignments[s]->init_with(o);
		return alignments[s];
	}

	void pushBackBA(const BamAlignment *o) {
		if (has_seen.insert(Signature(alignments[s])).second) {
			alignments[s++]->completeAlignment(o, has_seen.size() > 1);
		}
	}

	void asUnmap(const BamAlignment *o, const std::string& iv_type1, const std::string& iv_type2) {
		assert(s == 0);
		allocate();
		alignments[0]->asUnmap(o, iv_type1, iv_type2);
		++s;
	}

private:
	std::set<Signature> has_seen; // if we have seen this alignment
};

#endif
