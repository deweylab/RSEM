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

#ifndef GENOMEMAP_H_
#define GENOMEMAP_H_

#include <cstdio>
#include <string>
#include <vector>
#include <map>

struct Point {
	int pos; // genomic coordinate, 1-based.
	int tid; // transcript id; sign + start, - end
	int eid; // exon id; -1 refers to the whole transcript
   
	Point() : pos(0), tid(0), eid(-1) {}
	Point(int pos, int tid, int eid) : pos(pos), tid(tid), eid(eid) {}
   
	bool operator< (const Point& o) const {
		return pos < o.pos;
	}
};

struct Exon {
	int tid, eid;

	Exon() : tid(0), eid(0) {}
	Exon(int tid, int eid) : tid(tid), eid(eid) {}
};

struct Chrom {
	std::vector<int> coords; // genomic coordinates
	std::vector<std::vector<Exon> > exons;
	std::vector<std::string> types; // interval types: "intergenic_region", "intron", "exon"

	void clear() {
		coords.clear();
		exons.clear();
		types.clear();
	}

	void resize(int m) {
		clear();
		coords.resize(m);
		exons.resize(m);
		types.resize(m);
	}
};

class GenomeMap {
public:
	GenomeMap() { genomeMap.clear(); }

	void constructMap(const Transcripts* transcripts, const Transcripts* duplicates = NULL); 

	void readFrom(const char* inpF);
	void writeTo(const char* outF);

private:
	std::map<std::string, Chrom> genomeMap;
};
	
#endif /* GENOMEMAP_H_ */

