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
#include <cassert>
#include <cstdint> 
#include <cstring>
#include <cstdlib>
#include <string>
#include <sstream>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "my_assert.h"
#include "BamWriter.hpp"

// If header == NULL, just create an empty header with one target (to avoid free non null pointer due to calloc) and paste the program_id line
// if program_id == NULL, use the header passed
BamWriter::BamWriter(const char* outF, const bam_hdr_t* header, const char* program_id, htsThreadPool* p) {
	assert(header != NULL);
	this->header = (program_id == NULL ? bam_hdr_dup(header) : header_duplicate_without_text(header));

	if (program_id != NULL) {
		std::ostringstream strout;
		strout<< "@HD\tVN:1.4\tSO:unknown\n@PG\tID:"<< program_id<< std::endl;
		header_append_new_text(this->header, strout.str());
	}
	
	bam_out = sam_open(outF, "wb");
	general_assert(bam_out != 0, "Cannot write to " + cstrtos(outF) + "!");
	general_assert(sam_hdr_write(bam_out, this->header) == 0, "Cannot write header!");
	if (p != NULL) hts_set_opt(bam_out, HTS_OPT_THREAD_POOL, p);
}

BamWriter::~BamWriter() {
	bam_hdr_destroy(header);
	sam_close(bam_out);
}

bam_hdr_t* BamWriter::header_duplicate_without_text(const bam_hdr_t *ori_h) {
	bam_hdr_t *h = bam_hdr_init();
	h->n_targets = ori_h->n_targets;
	h->target_len = new uint32_t[h->n_targets];
	h->target_name = new char*[h->n_targets];
	for (int i = 0; i < h->n_targets; ++i) {
		h->target_len[i] = ori_h->target_len[i];
		h->target_name[i] = new char[strlen(ori_h->target_name[i]) + 1];
		strcpy(h->target_name[i], ori_h->target_name[i]);
	}
	return h;
}

void BamWriter::header_append_new_text(bam_hdr_t* header, const std::string& new_text) {
	if (new_text == "") return;
	int len = new_text.length();
	int max_size = (header->text == NULL ? 0 : header->l_text + 1);
	kroundup32(max_size);
	if (max_size < int(header->l_text + len + 1)) {
		max_size = header->l_text + len + 1;
		kroundup32(max_size);
		header->text = (char*)realloc(header->text, max_size);
	}
	strcpy(header->text + header->l_text, new_text.c_str());
	header->l_text += len;
}
