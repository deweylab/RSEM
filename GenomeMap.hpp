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

	bool operator< (const Exon& o) const {
		return tid < o.tid;
	}
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

typedef std::map<std::string, Chrom> genome_map;
typedef std::map<std::string, Chrom>::iterator genome_map_iter;

class GenomeMap {
public:
	GenomeMap() { genomeMap.clear(); empty.clear(); }

	void constructMap(const Transcripts* transcripts, const Transcripts* duplicates = NULL); 

	void readFrom(const char* inpF);
	void writeTo(const char* outF);

	// pos is one-based
	const std::vector<Exon>& getExonList(const std::string& chrname, int pos, std::string& iv_type) {
		g_it = genomeMap.find(chrname);
		if (g_it == genomeMap.end()) return empty;
		Chrom& chrom = g_it->second;
		int idx = binary_search(chrom.coords, pos);
		if (idx < 0) { iv_type = "intergenic_region"; return empty; }
		else iv_type = chrom.types[idx];
		return chrom.exons[idx];
	}

private:
	genome_map genomeMap;
	genome_map_iter g_it;

	std::vector<Exon> empty;

	// return the largest index r such as coords[r] <= pos
	int binary_search(const std::vector<int>& coords, int pos) {
		int l, r, mid;

		l = 0; r = coords.size() - 1;
		while (l <= r) {
			mid = (l + r) >> 1;
			if (pos >= coords[mid]) l = mid + 1;
			else r = mid - 1;
		}

		return r;
	}
};
	
#endif /* GENOMEMAP_H_ */
