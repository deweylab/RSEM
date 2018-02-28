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

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "my_assert.h"
#include "SamHeaderText.hpp"

class BamWriter {
public:
	/*
		@param   outF         output BAM file name
		@param   header       input header
		@param   p 			  pointer to htsThreadPool, NULL if no multi-threading.
	 */
	BamWriter(const char* outF, const bam_hdr_t* header, htsThreadPool* p = NULL);	
	~BamWriter();

	// switch between genome references and transcriptome references
	void replaceReferences(const char* faiF) {
		ht->replaceSQ(faiF);
	}

	// add program line
	void addProgram(const std::string& pid, const std::string& version = "", const std::string& command = "") {
		ht->addProgram(pid, version, command);
	}

	void writeHeader() {
		header = ht->create_header();
		general_assert(sam_hdr_write(bam_out, header) == 0, "Cannot write header!");
	}
	
	int name2id(const char* ref) { return bam_name2id(header, ref); }

	void write(bam1_t* b) {
		general_assert(sam_write1(bam_out, header, b) >= 0, "Fail to write alignments to BAM file!");    
	}

private:
	samFile* bam_out;
	bam_hdr_t* header;
	SamHeaderText* ht;
};

#endif
