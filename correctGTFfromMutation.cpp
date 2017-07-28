#include <cstdio>
#include <cassert>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

#include "utils.h"
#include "my_assert.h"
#include "GTFItem.h"

using namespace std;

bool verbose = true;

struct GenomicCoordinate {
	GenomicCoordinate(int x, bool isStart) : x(x), isStart(isStart) {}

	bool operator<(const GenomicCoordinate& o) const {
		if (x != o.x) return x < o.x;
		return isStart > o.isStart;
	}

	bool operator==(const GenomicCoordinate& o) const {
		return x == o.x && isStart == o.isStart;
	}

	int x;
	bool isStart;
};


class Exon {
public:
	Exon(const GTFItem& item, GenomicCoordinate* start, GenomicCoordinate* end) {
		this->gtf_item = item;
		this->start = start;
		this->end = end;
	}

	bool isValid() const { return start->x <= end->x; }

	string toString() {
		gtf_item.setStart(start->x);
		gtf_item.setEnd(end->x);

		return gtf_item.toString();
	}

private:
	GTFItem gtf_item;
	GenomicCoordinate *start, *end;
};

class FusionEvent {
public:
	FusionEvent(const string& outstr1, const string& outstr2, const string& outstr3, GenomicCoordinate* left, GenomicCoordinate* right) : outstr1(outstr1), outstr2(outstr2), outstr3(outstr3), left(left), right(right) {}

	string toString() {
		return outstr1 + itos(left->x) + outstr2 + itos(right->x) + outstr3;
	}

private:
	GenomicCoordinate *left, *right;	
	string outstr1, outstr2, outstr3;
};


typedef unordered_map<int, GenomicCoordinate*> CoordMap;
typedef CoordMap::iterator CoordMapIter;
typedef unordered_map<string, CoordMap> ChrMap;
typedef ChrMap::iterator ChrMapIter;
typedef unordered_map<string, unordered_set<string> > GnameMap;
typedef GnameMap::iterator GnameMapIter;

ChrMap chrMap;
vector<Exon> exons;
vector<FusionEvent> fusions;

GnameMap gname2gid;
GnameMapIter g2g_it;


GenomicCoordinate* insert2chrMap(const string& chr, int x, bool isStart) {
	ChrMapIter it1 = chrMap.find(chr);
	if (it1 == chrMap.end()) {
		it1 = chrMap.emplace(chr, CoordMap()).first;
	}

	CoordMap &coordMap = it1->second;
	int key = (x << 1) + isStart;
	CoordMapIter it2 = coordMap.find(key);
	if (it2 == coordMap.end()) {
		it2 = coordMap.emplace(key, new GenomicCoordinate(x, isStart)).first;
	}

	return it2->second;
}

void parseGTF(const char* gtfF) {
	ifstream fin(gtfF);
	string line;
	GTFItem item;

	general_assert(fin.is_open(), "Cannot open " + cstrtos(gtfF) + "! It may not exist.");

	int cnt = 0;

	int n_warns = 0;
	
 	while (getline(fin, line)) {
 		if (line[0] == '#') continue; // if this line is comment, jump it
 		item.parse(line);
  		if (item.getFeature() == "exon") {
 			if (item.getStart() > item.getEnd()) {
			  if (++n_warns <= MAX_WARNS) {
			    fprintf(stderr, "Warning: exon's start position is larger than its end position! This exon is discarded.\n");
			    fprintf(stderr, "\t%s\n\n", line.c_str());
			  }
 			}
 			else if (item.getStart() < 1) {
			  if (++n_warns <= MAX_WARNS) {
			    fprintf(stderr, "Warning: exon's start position is less than 1! This exon is discarded.\n");
			    fprintf(stderr, "\t%s\n\n", line.c_str());
			  }
 			}
 			else {
 				exons.emplace_back(item, \
 					insert2chrMap(item.getSeqName(), item.getStart(), true), \
 					insert2chrMap(item.getSeqName(), item.getEnd(), false));
 				
 				item.parseAttributes(line);
 				g2g_it = gname2gid.find(item.getGeneName());
 				if (g2g_it == gname2gid.end()) g2g_it = gname2gid.emplace(item.getGeneName(), unordered_set<string>()).first;
 				g2g_it->second.insert(item.getGeneID());
 			}
 		}

 		++cnt;
 		if (verbose && cnt % 200000 == 0) { printf("Parsed %d lines\n", cnt); }
	}
	fin.close();

	if (n_warns > 0) fprintf(stderr, "Warning: In total, %d exons are discarded.", n_warns);

	if (verbose) { printf("GTF file is parsed.\n"); }
}

inline string get_gids(unordered_set<string>& gid_set) {
	string gidstr = "";
	for (auto&& gid : gid_set) {
		gidstr += (gidstr == "" ? gid : "," + gid);
	}
	return gidstr;
}

void parseFusionFile(const char* fusionF) {
	ifstream fin(fusionF);
	string line, outstr1, outstr2, outstr3;
	GenomicCoordinate *left, *right;
	vector<string> fields, coords;

	general_assert(fin.is_open(), "Cannot open " + cstrtos(fusionF) + "! It may not exist.");

	getline(fin, line);
	while (getline(fin, line)) {
		split(line, '\t', fields);
		
		outstr1 = fields[1] + "\t" + fields[8] + "\t" + fields[9] + "\tONLY_REF_SPLICE\t" + fields[2] + "^";
		g2g_it = gname2gid.find(fields[2]);
		general_assert(g2g_it != gname2gid.end(), "Could not find gene id for " + fields[2] + "!");
		outstr1 += get_gids(g2g_it->second) + "\t";

		split(fields[3], ':', coords);
		outstr1 += coords[0] + ":";

		left = insert2chrMap(coords[0], stoi(coords[1]), false);

		outstr2 = ":" + coords[2] + "\t" + fields[5] + "^";
		g2g_it = gname2gid.find(fields[5]);
		general_assert(g2g_it != gname2gid.end(), "Could not find gene id for " + fields[5] + "!");		
		outstr2 += get_gids(g2g_it->second) + "\t";

		split(fields[6], ':', coords);
		outstr2 += coords[0] + ":";

		right = insert2chrMap(coords[0], stoi(coords[1]), true);

		outstr3 = ":" + coords[2];

		fusions.emplace_back(outstr1, outstr2, outstr3, left, right);
	}

	fin.close();

	if (verbose) { printf("Fusion file is parsed.\n"); }
}


inline bool cmp(const GenomicCoordinate* a, const GenomicCoordinate* b) {
	return (*a) < (*b);
}

void update(const string& chr, string& cur_chr, int& last_pos, int& delta, int& ptr, vector<GenomicCoordinate*>& points) {
	while (ptr < points.size()) points[ptr++]->x += delta;
	
	points.clear();
	ChrMapIter it = chrMap.find(chr);
	general_assert(it != chrMap.end(), "Could not find " + chr + "!");
	for (auto&& kv : it->second) {
		points.push_back(kv.second);
	}
	sort(points.begin(), points.end(), cmp);

	cur_chr = chr; 
	last_pos = delta = ptr = 0;
}

void processMutations(const char* mutF) {
	ifstream fin(mutF);
	string line;
	vector<string> fields;

	string chr;
	bool isINS;
	int pos, len;

	string cur_chr;
	int last_pos, delta, ptr; // delta, the difference between new and old coordinate; ptr, current position at points
	vector<GenomicCoordinate*> points;

	general_assert(fin.is_open(), "Cannot open " + cstrtos(mutF) + "! It may not exist.");

	cur_chr = "";
	last_pos = delta = ptr = 0;

	int cnt = 0;

	getline(fin, line);
	while (getline(fin, line)) {
		split(line, '\t', fields);
		chr = "chr" + fields[4];
		pos = stoi(fields[5]);
		isINS = (fields[9] == "INS");
		assert(fields[11] == fields[12]);
		len = (isINS ? fields[11].length() : fields[10].length());

		if (cur_chr != chr) update(chr, cur_chr, last_pos, delta, ptr, points);

		assert(last_pos < pos);

		// Assume that the pos for INS will never be deleted.
		if (isINS) {
			while (ptr < points.size() && points[ptr]->x < pos) points[ptr++]->x += delta;
			while (ptr < points.size() && points[ptr]->x == pos) points[ptr]->x += delta + (points[ptr]->isStart ? 0 : len), ++ptr;
			while (ptr < points.size() && points[ptr]->x == pos + 1 && points[ptr]->isStart) points[ptr++]->x += delta;
			delta += len;
			last_pos = pos;
		}
		else {
			while (ptr < points.size() && points[ptr]->x < pos) points[ptr++]->x += delta;
			while (ptr < points.size() && points[ptr]->x < pos + len) points[ptr]->x = (points[ptr]->isStart ? pos + delta : pos + delta - 1), ++ptr;
			delta -= len;
			last_pos = pos + len - 1;
		}

		if (verbose && ++cnt % 200000 == 0) printf("%d lines parsed.\n", cnt);
	}

	while (ptr < points.size()) points[ptr++]->x += delta;

	fin.close();

	if (verbose) { printf("processMutations is finished.\n"); }
}

void output(const char* outF) {
	ofstream fout(outF);
	for (auto&& exon : exons) {
		if (exon.isValid()) fout<< exon.toString()<< endl;
	}
	fout.close();
	if (verbose) printf("%s is generated.\n", outF);
}

void outputFusion(const string& fusionF) {
	string outF = fusionF + ".sf";
	ofstream fout(outF);
	for (auto&& fusion : fusions) fout<< fusion.toString()<< endl;
	fout.close();
	if (verbose) printf("%s is generated.\n", outF.c_str());
}

void release() {
	for (auto&& kv1 : chrMap) 
		for (auto&& kv2 : kv1.second)
			delete kv2.second; 
}

int main(int argc, char* argv[]) {
	if (argc != 4 && argc != 6) {
		printf("Usage: correctGTFfromMutation input.GTF input.mutation output.GTF [--fusion fusionF]\n");
		exit(-1);
	}

	parseGTF(argv[1]);

	if (argc == 6) parseFusionFile(argv[5]);

	processMutations(argv[2]);

	output(argv[3]);
	if (argc == 6) outputFusion(argv[5]);

	release();

	return 0;
}
