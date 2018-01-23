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

#ifndef SAMPARSER_H_
#define SAMPARSER_H_

#include <cctype>
#include <cassert>
#include <string>
#include <vector>
#include <map>

#include "htslib/sam.h"

#include "utils.h"

class SamParser {
public:
	static void buildRefHeader(const char* transListF);
	static void releaseRefHeader();

	// omitF contains RSEM tids that never appeared in any of the alignment files
	// FORMAT: one tid per line
	static void writeOmitFile(const char* omitF); 


	SamParser(const char* inpF, bool build_mapping = false); 
	~SamParser();

	const bam_hdr_t* getHeader() const { 
		return header;
	}

	const char* getProgramID(); // scan header to look up program ID, and convert it to lower case, slow

	bool read(bam1_t* b) {
		if (sam_read1(in, header, b) < 0) return false;
		if (remap) {
			if (!(b->core.flag & BAM_FUNMAP)) b->core.tid = get_tid(b->core.tid);
			if ((b->core.flag & BAM_FPAIRED) && !(b->core.flag & BAM_FMUNMAP)) b->core.mtid = get_tid(b->core.mtid);
		}
		return true;
	}
	
private:
	static const bam_hdr_t *ref_header = NULL;
	static std::map<std::string, int> tname2tid; // mapping from transcript name to transcript id
	static std::map<std::string, int>::iterator iter;
	static std::vector<bool> appeared; // vector recording if a RSEM tid appeared in the alignment file

	samFile *sam_in;
	bam_hdr_t *header;

	bool delete_header;
	char program_id[STRLEN];

	bool remap; // if the input BAM header is not consistent with the RSEM BAM header
	std::vector<int> sid2tid; // mapping from input BAM sequence id to RSEM BAM transcript id

	void buildMapping();

	// 0-based
	int get_tid(int sid) const {
		return sid2tid[sid];
	}

	bool set_program_id(const char *fr, const char *p) {
		if (p - fr > 3 && !strncmp(fr, "ID:", 3)) {
			fr += 3;
			assert(p - fr <= 100);
			for (int i = 0; i < p - fr; ++i)
				program_id[i] = tolower(fr[i]);
			return true;
		}
		return false;
	}
};

#endif
