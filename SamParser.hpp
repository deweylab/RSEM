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
#include <cstdint>
#include <cassert>
#include <string>
#include <vector>
#include <map>

#include "htslib/bgzf.h"
#include "htslib/sam.h"

#include "SamHeaderText.hpp"

class SamParser {
public:
	static void buildMapping(const char* transListF, const bam_hdr_t* header, const char* remapF);
	static void loadMapping(const char* remapF);

	SamParser(const char* inpF, htsThreadPool* p = NULL); 
	~SamParser();

	const bam_hdr_t* getHeader() const { 
		return header;
	}

	const std::string& getProgramID() { 
		return ht->getProgramID(); 
	}

	bool read(bam1_t* b) {
		if (sam_read1(sam_in, header, b) < 0) return false;
		if (remap) {
			if (!(b->core.flag & BAM_FUNMAP)) b->core.tid = get_tid(b->core.tid);
			if ((b->core.flag & BAM_FPAIRED) && !(b->core.flag & BAM_FMUNMAP)) b->core.mtid = get_tid(b->core.mtid);
		}
		return true;
	}
	
	int64_t tell() { return bgzf_tell(sam_in->fp.bgzf); }

	void seek(int64_t pos) {
		assert(bgzf_seek(sam_in->fp.bgzf, pos, SEEK_SET) == 0);
	}

	const char* getSeqName(int tid) {
		return header->target_name[tid];
	}

private:
	static bool remap; // if the input BAM header is not consistent with the RSEM BAM header
	static std::vector<int> sid2tid; // mapping from input BAM sequence id to RSEM BAM transcript id

	// 0-based
	static int get_tid(int sid) {
		return sid2tid[sid];
	}


	samFile* sam_in;
	bam_hdr_t* header;
	SamHeaderText* ht;
};

#endif
