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

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <string>
#include <algorithm>

#include <stdint.h>
#include "htslib/sam.h"
#include "my_assert.h"
#include "SamParser.hpp"
#include "BamWriter.hpp"
#include "BamBufferedWriters.hpp"
#include "BamAlignment.hpp"

const uint8_t BamAlignment::rnt_table[16] = {0, 8, 4, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 15};

bool BamAlignment::read(SamParser *in, BamAlignment *o) {
	is_aligned = -1;

	if (b == NULL) b = bam_init1();
	if (!in->read(b)) return false;


	is_paired = bam_is_paired(b);
	// read the second mate
	if (is_paired) { 
		if (b2 == NULL) b2 = bam_init1();

		general_assert(in->read(b2) && bam_is_paired(b2), "Fail to read the other mate for a paired-end alignment!");
		general_assert(((b->core.flag & 0x00C0) == 0x0040 && (b2->core.flag & 0x00C0) == 0x0080) || 
			 ((b->core.flag & 0x00C0) == 0x0080 && (b2->core.flag & 0x00C0) == 0x0040), 
			 "Cannot detect both mates of a paired-end alignment!");
		
		if (bam_is_read2(b)) { bam1_t* tmp = b; b = b2; b2 = tmp; }
	}

	// calculate is_aligned
	is_aligned = bam_is_mapped(b);
	if (is_paired) is_aligned |= ((char)bam_is_mapped(b2) << 1);
		
	// The following four statements are grouped together
	int tmp_len = 0, tmp_len2 = 0;
	tmp_len = b->core.l_qseq <= 0 ? o->getSeqLength(1) : b->core.l_qseq;
	if (is_paired) tmp_len2 = b2->core.l_qseq <= 0 ? o->getSeqLength(2) : b2->core.l_qseq;
	assert(!(is_aligned & 1) || tmp_len == bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)));
	assert(!(is_aligned & 2) || tmp_len2 == bam_cigar2qlen(b2->core.n_cigar, bam_get_cigar(b2)));

	return true;
}

// choice: 0, do nothing; 1, delete read sequence and qual score; 2, add read sequence and qual score. o, the alignment that contain read sequence and quality score information
bool BamAlignment::write(BamWriter *out, int choice, BamAlignment *o) {
	assert(is_aligned >= 0 && b != NULL && (!is_paired || b2 != NULL));

	if (b->core.l_qname == 1) b->core.l_qseq = 0;
	if (is_paired && (b2->core.l_qname == 1)) b2->core.l_qseq = 0;

	switch(choice) {
	case 0: 
		break;
	case 1:
		if (b->core.l_qname > 1) compress(b);
		if (is_paired && (b2->core.l_qname > 1)) compress(b2);
		break;
	case 2:
		if (b->core.l_qname == 1) decompress(b, o->b);
		if (is_paired && (b2->core.l_qname == 1)) decompress(b2, o->b2);
		break;
	default: assert(false);
	}

	out->write(b);
	if (is_paired) out->write(b2);

	return true;
}

bool BamAlignment::writeToBuffer(BamBufferedWriters* buffer, int id, int choice) {
	assert(is_aligned >= 0 && b != NULL && (!is_paired || b2 != NULL));
	
	if (b->core.l_qname == 1) b->core.l_qseq = 0;
	if (is_paired && (b2->core.l_qname == 1)) b2->core.l_qseq = 0;

	bam1_t *target = NULL;
	
	target = buffer->get_bam1_t(id);
	copy_bam1_t(target, b, choice);

	if (is_paired) {
		target = buffer->get_bam1_t(id);
		copy_bam1_t(target, b2, choice);
	}

	return true;
}

void BamAlignment::compress(bam1_t* b) {
	int l_aux = bam_get_l_aux(b);
	b->data[0] = 0;
	memmove(b->data + 1, bam_get_cigar(b), b->core.n_cigar * 4);
	memmove(b->data + 1 + b->core.n_cigar * 4, bam_get_aux(b), l_aux);
	b->l_data = 1 + b->core.n_cigar * 4 + l_aux;
	b->core.l_qname = 1;
	b->core.l_qseq = 0;
}

void BamAlignment::decompress(bam1_t* b, bam1_t* other) {
	int l_aux = bam_get_l_aux(b);
	b->core.l_qname = other->core.l_qname;
	b->core.l_qseq = other->core.l_qseq;
	b->l_data = b->core.l_qname + b->core.n_cigar * 4 + (b->core.l_qseq + 1) / 2 + b->core.l_qseq + l_aux;
	expand_data_size(b);
	memmove(bam_get_aux(b), b->data + 1 + b->core.n_cigar * 4, l_aux); // move aux options
	memmove(bam_get_cigar(b), b->data + 1, b->core.n_cigar * 4); // move cigar string
	memcpy(bam_get_qname(b), bam_get_qname(other), b->core.l_qname); // copy qname

	if (bam_is_rev(b) == bam_is_rev(other)) {
		memcpy(bam_get_seq(b), bam_get_seq(other), (b->core.l_qseq + 1) / 2); // copy seq
		memcpy(bam_get_qual(b), bam_get_qual(other), b->core.l_qseq); // copy qual
	}
	else {
		copy_rc_seq(bam_get_seq(b), bam_get_seq(other), b->core.l_qseq); // copy reverse complement seq
		copy_r_qual(bam_get_qual(b), bam_get_qual(other), b->core.l_qseq); // copy reverse qual
	}
}

// choice can only be 0 or 1
void BamAlignment::copy_bam1_t(bam1_t *target, bam1_t *b, int choice) {
	target->core = b->core;
	if (choice == 0 || b->core.l_qname <= 1) { // copy
		target->l_data = b->l_data;
		expand_data_size(target);
		memcpy(target->data, b->data, target->l_data);
	}
	else { // compress
		int l_aux = bam_get_l_aux(b);
		target->l_data = 1 + b->core.n_cigar * 4 + l_aux;
		expand_data_size(target);
		target->data[0] = 0;
		memcpy(target->data + 1, bam_get_cigar(b), target->core.n_cigar * 4);
		memcpy(target->data + 1 + target->core.n_cigar * 4, bam_get_aux(b), l_aux);
		target->core.l_qname = 1;
		target->core.l_qseq = 0;
	}  
}
