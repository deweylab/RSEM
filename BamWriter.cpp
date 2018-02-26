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

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "my_assert.h"
#include "SamHeaderText.hpp"
#include "BamWriter.hpp"

BamWriter::BamWriter(const char* outF, const bam_hdr_t* header, htsThreadPool* p) {
	bam_out = sam_open(outF, "wb");
	general_assert(bam_out != 0, "Cannot write to " + cstrtos(outF) + "!");
	ht = new SamHeaderText(header);
	if (p != NULL) hts_set_opt(bam_out, HTS_OPT_THREAD_POOL, p);
}

BamWriter::~BamWriter() {
	delete ht;
	bam_hdr_destroy(header);
	sam_close(bam_out);
}
