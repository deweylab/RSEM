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
#include <cctype>
#include <cstring>
#include <cstdint>
#include <cassert>
#include <string>
#include <algorithm>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "utils.h"
#include "my_assert.h"
#include "SamParser.hpp"
#include "BamWriter.hpp"
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

		if (!(((b->core.flag & 0x00C0) == 0x0040 && (b2->core.flag & 0x00C0) == 0x0080) || 
			 ((b->core.flag & 0x00C0) == 0x0080 && (b2->core.flag & 0x00C0) == 0x0040))) {
			printf("%s\n", bam_get_qname(b));
			printf("%s\n", bam_get_qname(b2));
		}

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
		if (b->core.l_qname - b->core.l_extranul > 1) compress(b);
		if (is_paired && (b2->core.l_qname - b2->core.l_extranul > 1)) compress(b2);
		break;
	case 2:
		if (b->core.l_qname - b->core.l_extranul == 1) decompress(b, o->b);
		if (is_paired && (b2->core.l_qname - b->core.l_extranul == 1)) decompress(b2, o->b2);
		break;
	default: assert(false);
	}

	out->write(b);
	if (is_paired) out->write(b2);

	return true;
}

void BamAlignment::setConvertedInfo(int mate, int32_t tid, int32_t pos, bool is_rev, uint32_t n_cigar, const uint32_t* cigar) {
	assert(mate == 1 || (is_paired && mate == 2));
	bam1_t *bt = (mate == 1 ? b : b2);
	bt->core.tid = tid;
	bt->core.pos = pos;
	bt->core.flag = (bt->core.flag ^ (bt->core.flag & BAM_FREVERSE)) | (is_rev ? BAM_FREVERSE : 0);
	bt->l_data = bt->l_data - (bt->core.n_cigar<< 2) + (n_cigar << 2);
	bt->core.n_cigar = n_cigar;
	expand_data_size(bt);
	if (is_rev) copy_r_cigar(bam_get_cigar(bt), cigar, n_cigar);
	else memcpy(bam_get_cigar(bt), cigar, sizeof(uint32_t) * n_cigar);
}

void BamAlignment::completeAlignment(const BamAlignment* o, bool is_secondary) {
	transfer(b, o->b, is_secondary);
	if (is_paired) transfer(b2, o->b2, is_secondary);
	
	if (is_aligned == 3) {
		b->core.mtid = b2->core.tid; b2->core.mtid = b->core.tid;
		b->core.mpos = b2->core.pos; b2->core.mpos = b->core.pos;
		if (b->core.pos < b2->core.pos) { b->core.isize = bam_endpos(b2) - b->core.pos; b2->core.isize = -b->core.isize; }
		else { b2->core.isize = bam_endpos(b) - b2->core.pos; b->core.isize = -b2->core.isize; }
		b->core.flag = (b->core.flag ^ (b->core.flag & BAM_FMREVERSE)) | (bam_is_rev(b2) ? BAM_FMREVERSE : 0);
		b2->core.flag = (b2->core.flag ^ (b2->core.flag & BAM_FMREVERSE)) | (bam_is_rev(b) ? BAM_FMREVERSE : 0);
	}
}

void BamAlignment::asUnmap(const BamAlignment* o, const std::string& str1, const std::string& str2) {
	is_paired = o->is_paired;
	is_aligned = 0;
	if (b == NULL) b = bam_init1();
	if (is_paired && b2 == NULL) b2 = bam_init1();

	b->core.tid = b->core.pos = b->core.mtid = b->core.mpos = -1; b->core.isize = 0;
	b->core.bin = 0; b->core.qual = b->core.unused1 = 0; b->core.n_cigar = 0;
	b->core.l_qname = o->b->core.l_qname; b->core.l_extranul = o->b->core.l_extranul; b->core.l_qseq = o->b->core.l_qseq;
	b->core.flag = (o->b->core.flag & 0x06C1) | BAM_FUNMAP | (is_paired ? BAM_FMUNMAP : 0); // 0x06C1 = 0x0001, 0x0040, 0x0080, 0x0200, 0x0400
	b->l_data = b->core.l_qname + ((b->core.l_qseq + 1)>> 1) + b->core.l_qseq;
	expand_data_size(b);
	memcpy(bam_get_qname(b), bam_get_qname(o->b), b->core.l_qname); // copy qname
	if (!bam_is_rev(o->b)) {
		memcpy(bam_get_seq(b), bam_get_seq(o->b), ((b->core.l_qseq + 1)>> 1) + b->core.l_qseq); // copy seq + qual
	}
	else {
		copy_rc_seq(bam_get_seq(b), bam_get_seq(o->b), b->core.l_qseq); // copy reverse complement seq
		copy_r_qual(bam_get_qual(b), bam_get_qual(o->b), b->core.l_qseq); // copy reverse qual
	}
	bam_aux_append(b, "ZG", 'Z', str1.length() + 1, (uint8_t*)str1.c_str()); // append read type tag

	if (is_paired) {
		b2->core.tid = b2->core.pos = b2->core.mtid = b2->core.mpos = -1; b2->core.isize = 0;
		b2->core.bin = 0; b2->core.qual = b2->core.unused1 = 0; b2->core.n_cigar = 0;
		b2->core.l_qname = o->b2->core.l_qname; b2->core.l_extranul = o->b2->core.l_extranul; b2->core.l_qseq = o->b2->core.l_qseq;
		b2->core.flag = (o->b2->core.flag & 0x06C1) | BAM_FUNMAP | (is_paired ? BAM_FMUNMAP : 0); // 0x06C1 = 0x0001, 0x0040, 0x0080, 0x0200, 0x0400
		b2->l_data = b2->core.l_qname + ((b2->core.l_qseq + 1)>> 1) + b2->core.l_qseq;
		expand_data_size(b2);
		memcpy(bam_get_qname(b2), bam_get_qname(o->b2), b2->core.l_qname); // copy qname
		if (!bam_is_rev(o->b2)) {
			memcpy(bam_get_seq(b2), bam_get_seq(o->b2), ((b2->core.l_qseq + 1)>> 1) + b2->core.l_qseq); // copy seq + qual
		}
		else {
			copy_rc_seq(bam_get_seq(b2), bam_get_seq(o->b2), b2->core.l_qseq); // copy reverse complement seq
			copy_r_qual(bam_get_qual(b2), bam_get_qual(o->b2), b2->core.l_qseq); // copy reverse qual			
		}
		bam_aux_append(b2, "ZG", 'Z', str2.length() + 1, (uint8_t*)str2.c_str()); // append read type tag
	}
}

void BamAlignment::compress(bam1_t* b) {
	int l_aux = bam_get_l_aux(b);
	memset(b->data, 0, 4);
	memmove(b->data + 4, bam_get_cigar(b), b->core.n_cigar * 4);
	memmove(b->data + 4 + b->core.n_cigar * 4, bam_get_aux(b), l_aux);
	b->l_data = 4 + b->core.n_cigar * 4 + l_aux;
	b->core.l_qname = 4;
	b->core.l_extranul = 3;
	b->core.l_qseq = 0;
}

void BamAlignment::decompress(bam1_t* b, const bam1_t* other) {
	int l_aux = bam_get_l_aux(b);
	b->core.l_qname = other->core.l_qname;
	b->core.l_extranul = other->core.l_extranul;
	b->core.l_qseq = other->core.l_qseq;
	b->l_data = b->core.l_qname + b->core.n_cigar * 4 + (b->core.l_qseq + 1) / 2 + b->core.l_qseq + l_aux;
	expand_data_size(b);
	memmove(bam_get_aux(b), b->data + 4 + b->core.n_cigar * 4, l_aux); // move aux options
	memmove(bam_get_cigar(b), b->data + 4, b->core.n_cigar * 4); // move cigar string
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

void BamAlignment::transfer(bam1_t* b, const bam1_t* other, bool is_secondary) {
	memcpy(bam_get_qname(b), bam_get_qname(other), b->core.l_qname); // copy qname
	if (bam_is_unmapped(b) || bam_is_rev(b) == bam_is_rev(other)) {
		memcpy(bam_get_seq(b), bam_get_seq(other), other->l_data - other->core.l_qname - (other->core.n_cigar << 2)); // copy seq-qual-aux
	}
	else {
		copy_rc_seq(bam_get_seq(b), bam_get_seq(other), other->core.l_qseq); // copy reverse complement seq
		copy_r_qual(bam_get_qual(b), bam_get_qual(other), other->core.l_qseq); // copy reverse qual
		memcpy(bam_get_aux(b), bam_get_aux(other), bam_get_l_aux(other)); // copy aux
		reverse_MD(b);
	}
	if (bam_is_mapped(b)) {
		if (is_secondary) b->core.flag |= BAM_FSECONDARY;
		b->core.bin = hts_reg2bin(b->core.pos, b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)), 14, 5);
	}
}

void BamAlignment::reverse_MD(bam1_t* b) {
	uint8_t* p = bam_aux_get(b, "MD");
	if (p == NULL) return;
	char* mdstr = bam_aux2Z(p);
	assert(mdstr != 0);

	const uint32_t* cigar = bam_get_cigar(b);

	int len = strlen(mdstr), fr, to;
	int rpos, rsize, rc_pos;
	char* buffer = strdup(mdstr);

	fr = to = 0; rpos = len - 1; rc_pos = b->core.n_cigar - 1;
	while (to < len) {
		if (isdigit(buffer[to])) {
			fr = to++;
			while (isdigit(buffer[to])) ++to;
			rsize = to - fr;
			for (int i = 1; i <= rsize; ++i) mdstr[rpos - rsize + i] = buffer[fr + i - 1];
			rpos -= rsize;	
		} 
		else if (buffer[to] == '^') {
			while (rc_pos >= 0 && bam_cigar_op(cigar[rc_pos]) != BAM_CDEL) --rc_pos;
			assert(rc_pos >= 0);

			fr = ++to; 
			rsize = bam_cigar_oplen(cigar[rc_pos]);
			to += rsize;
			for (int i = 0; i < rsize; ++i) {
				assert(isalpha(buffer[fr + i]));
				mdstr[rpos--] = base2rbase[buffer[fr + i]];
			}
			mdstr[rpos--] = '^';
		}
		else {
			assert(isalpha(buffer[to]));
			mdstr[rpos--] = base2rbase[buffer[to++]];
		}
	}
	free(buffer);
}
