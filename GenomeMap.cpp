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

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <fstream>

#include "my_assert.h"
#include "Transcript.hpp"
#include "Transcripts.hpp"
#include "GenomeMap.hpp"

typedef std::map<std::string, std::vector<Point> > chrom_type;
typedef std::map<std::string, std::vector<Point> >::iterator chrom_type_iter;
typedef std::pair<std::string, std::vector<Point> > value_type_pair;

void GenomeMap::constructMap(const Transcripts* transcripts, const Transcripts* duplicates) {
	int s;
	int M = transcripts->getM();

	chrom_type chroms;
	value_type_pair value_pair;

	value_pair.second.clear();

	for (int i = 1; i <= M; ++i) {
		const Transcript& transcript = transcripts->getTranscriptAt(i);
		value_pair.first = transcript.getSeqName();
		std::vector<Point>& vec = chroms.insert(value_pair).first->second;
		const std::vector<Interval>& structure = transcript.getStructure();
		s = structure.size();
		for (int j = 0; j < s; ++j) {
			vec.push_back(Point(structure[j].start, i, j));
			vec.push_back(Point(structure[j].end + 1, -i, j));
		}
		vec.push_back(Point(structure[0].start, i, -1));
		vec.push_back(Point(structure[s - 1].end + 1, -i, -1)); // push transcript boundary
	}

	if (duplicates != NULL) {
		int DM = duplicates->getM();
		for (int i = 1; i <= DM; ++i) {
			const Transcript& transcript = duplicates->getTranscriptAt(i);
			value_pair.first = transcript.getSeqName();
			std::vector<Point>& vec = chroms.insert(value_pair).first->second;
			const std::vector<Interval>& structure = transcript.getStructure();
			s = structure.size();
			for (int j = 0; j < s; ++j) {
				vec.push_back(Point(structure[j].start, i + M, j));
				vec.push_back(Point(structure[j].end + 1, -(i + M), j));
			}
			vec.push_back(Point(structure[0].start, i + M, -1));
			vec.push_back(Point(structure[s - 1].end + 1, -(i + M), -1)); // push transcript boundary
		}		
	}

	int transcript_count; // how many transcript cover current location
	int n_exon; // number of forward and reverse strand exons covering current position
	int tid, sign;
	std::pair<std::string, Chrom> gm_pair;
	std::list<Exon> exon_list;
	std::list<Exon>::iterator iter;
	std::vector<Exon> dummy_vec;

	genomeMap.clear();
	exon_list.clear();
	dummy_vec.clear();
	transcript_count = n_exon = 0;
	for (chrom_type_iter it = chroms.begin(); it != chroms.end(); ++it) {
		std::vector<Point>& vec = it->second;
		s = vec.size();
		std::sort(vec.begin(), vec.end());

		gm_pair.first = it->first;
		Chrom& chrom = genomeMap.insert(gm_pair).first->second;
		chrom.clear();

		for (int i = 0; i < s; ++i) {
			sign = vec[i].tid > 0;
			tid = abs(vec[i].tid);

			if (sign) { // start
				if (vec[i].eid < 0) ++transcript_count;
				else {
					exon_list.push_back(Exon(tid, vec[i].eid));
					++n_exon;
				}
			}
			else { // end
				if (vec[i].eid < 0) --transcript_count;
				else {
					for (iter = exon_list.begin(); iter != exon_list.end(); ++iter) 
						if (iter->tid == tid && iter->eid == vec[i].eid) {
							exon_list.erase(iter);
							break;
						}
					assert(iter != exon_list.end());
					--n_exon;
				}
			}

			if (i + 1 == s || vec[i].pos < vec[i + 1].pos) {
				chrom.coords.push_back(vec[i].pos);
				chrom.exons.push_back(dummy_vec);
				if (n_exon > 0) {
					chrom.exons.back().assign(exon_list.begin(), exon_list.end());
					chrom.types.push_back("exon");
				}
				else chrom.types.push_back(transcript_count > 0 ? "intron" : "intergenic_region");
			}
		}
	}
}

void GenomeMap::readFrom(const char* inpF) {
	int n, m, s;
	std::string chrname;
	std::ifstream fin(inpF);
	std::pair<std::string, Chrom> gm_pair;

	genomeMap.clear();
	general_assert(fin.is_open(), "Cannot open " + cstrtos(inpF) + "! It may not exist.");

	fin>> n;
	for (int i = 0; i < n; ++i) {
		fin>> chrname>> m;
		gm_pair.first = chrname;
		Chrom& chrom = genomeMap.insert(gm_pair).first->second;
		chrom.resize(m);
		for (int j = 0; j < m; ++j) {
			fin>> chrom.coords[j]>> chrom.types[j]>> s;
			chrom.exons[j].resize(s);
			for (int k = 0; k < s; ++k) fin>> chrom.exons[j][k].tid>> chrom.exons[j][k].eid;
		}
	}
}

bool cmp(const genome_map_iter& a, const genome_map_iter& b) {
	return a->first < b->first;
}

void GenomeMap::writeTo(const char* outF) {
	std::vector<genome_map_iter> iters;
	std::ofstream fout(outF);

	iters.clear();
	for (genome_map_iter it = genomeMap.begin(); it != genomeMap.end(); ++it) {
		iters.push_back(it);
	}
	std::sort(iters.begin(), iters.end(), cmp);

	fout<< genomeMap.size()<< std::endl;
	for (std::size_t i = 0; i < iters.size(); ++i) {
		Chrom& chrom = iters[i]->second;

		fout<< std::endl;
		fout<< iters[i]->first<< '\t'<< chrom.coords.size()<< std::endl;
		for (std::size_t j = 0; j < chrom.coords.size(); ++j) {
			fout<< chrom.coords[j]<< '\t'<< chrom.types[j]<< '\t'<< chrom.exons[j].size();
			for (std::size_t k = 0; k < chrom.exons[j].size(); ++k) fout<< ' '<< chrom.exons[j][k].tid<< ' '<< chrom.exons[j][k].eid;			
			fout<< std::endl;
		}
	}
}
