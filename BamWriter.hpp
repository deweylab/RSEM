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

#ifndef BAMWRITER_H_
#define BAMWRITER_H_

#include <string>
#include "htslib/sam.h"
#include "my_assert.h"

class BamWriter {
public:
	/*
		@param   outF         output BAM file name
		@param   header       BAM header for output
		@param   n_threads     number of threads used to generate the BAM file, the number of compressing threads is n_threads - 1
	 */
	BamWriter(const char* outF, const bam_hdr_t* header, const char* program_id = NULL, int n_threads = 1);	
	~BamWriter();

	void write(bam1_t* b) {
		general_assert(sam_write1(bam_out, header, b) >= 0, "Fail to write alignments to BAM file!");    
	}

private:
	samFile* bam_out;
	bam_hdr_t* header;

	bam_hdr_t* header_duplicate_without_text(const bam_hdr_t* ori_h);
	void header_append_new_text(bam_hdr_t* header, const std::string& new_text);
};

#endif
