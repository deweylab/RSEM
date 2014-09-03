/* ReadType here means if the read is unalignable, alignable or aligned too much. It is NOT siheaderngle read or paired-end read */
#ifndef SAMPARSER_H_
#define SAMPARSER_H_

#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>

#include "sam/bam.h"
#include "sam/sam.h"

#include "utils.h"
#include "my_assert.h"

#include "SingleRead.h"
#include "SingleReadQ.h"
#include "PairedEndRead.h"
#include "PairedEndReadQ.h"
#include "SingleHit.h"
#include "PairedEndHit.h"

#include "Transcripts.h"

class SamParser {
public:
	SamParser(char, const char*, const char*, Transcripts&, const char*);
	~SamParser();

	/**
	 * return value
	 * -1 : no more alignment
	 * 0 : new read , type 0
	 * 1 : new read , type 1 with alignment
	 * 2 : new read , type 2
	 * 4 : discard this read
	 * 5 : new alignment but same read
	 */
	int parseNext(SingleRead&, SingleHit&);
	int parseNext(SingleReadQ&, SingleHit&);
	int parseNext(PairedEndRead&, PairedEndHit&);
	int parseNext(PairedEndReadQ&, PairedEndHit&);

	static void setReadTypeTag(const char* tag) {
		strcpy(rtTag, tag);
	}

private:
	samfile_t *sam_in;
	bam_header_t *header;
	bam1_t *b, *b2;

	Transcripts& transcripts;

	//tag used by aligner
	static char rtTag[STRLEN];

	//1 +, -1 - here
	int getDir(const bam1_t* b) {
		return ((b->core.flag & 0x0010) ? -1 : 1);
	}

	std::string getName(const bam1_t* b) {
		return std::string(bam1_qname(b));
	}

	std::string getReadSeq(const bam1_t*);
	std::string getQScore(const bam1_t*);

	//0 ~ N0 1 ~ N1 2 ~ N2
	int getReadType(const bam1_t*);
	int getReadType(const bam1_t*, const bam1_t*); // for paired-end reads

	bool check(bam1_t *b) {
		return (b->core.n_cigar == 1) && ((*bam1_cigar(b) & BAM_CIGAR_MASK) == BAM_CMATCH) && (b->core.l_qseq == (int32_t)(*bam1_cigar(b) >> BAM_CIGAR_SHIFT));
	}
};

char SamParser::rtTag[STRLEN] = ""; // default : no tag, thus no Type 2 reads

// aux, if not 0, points to the file name of fn_list
SamParser::SamParser(char inpType, const char* inpF, const char* aux, Transcripts& transcripts, const char* imdName)
	: transcripts(transcripts)
{
	switch(inpType) {
	case 'b': sam_in = samopen(inpF, "rb", aux); break;
	case 's': sam_in = samopen(inpF, "r", aux); break;
	default: assert(false);
	}

	general_assert(sam_in != 0, "Cannot open " + cstrtos(inpF) + "! It may not exist.");
	header = sam_in->header;
	general_assert(header != 0, "Fail to parse sam header!");

	transcripts.buildMappings(header->n_targets, header->target_name, imdName);

	b = bam_init1();
	b2 = bam_init1();
}

SamParser::~SamParser() {
	samclose(sam_in);
	bam_destroy1(b);
	bam_destroy1(b2);
}

// If sam_read1 returns 0 , what does it mean?
//Assume b.core.tid is 0-based
int SamParser::parseNext(SingleRead& read, SingleHit& hit) {
	int val; // return value
	bool canR = (samread(sam_in, b) >= 0);
	if (!canR) return -1;

	general_assert(!(b->core.flag & 0x0001), "Find a paired end read in the file!");

	int readType = getReadType(b);
	std::string name = getName(b);

	if (readType != 1 || (readType == 1 && read.getName().compare(name) != 0)) {
		val = readType;
		read = SingleRead(name, getReadSeq(b));
	}
	else val = 5;

	if (readType == 1) {
		if (!check(b)) { fprintf(stderr, "RSEM does not support gapped alignments, sorry!\n"); exit(-1); }

		if (getDir(b) > 0) {
			hit = SingleHit(transcripts.getInternalSid(b->core.tid + 1), b->core.pos);
		}
		else {
			hit = SingleHit(-transcripts.getInternalSid(b->core.tid + 1), header->target_len[b->core.tid] - b->core.pos - b->core.l_qseq);
		}
	}

	return val;
}

int SamParser::parseNext(SingleReadQ& read, SingleHit& hit) {
	int val;
	bool canR = (samread(sam_in, b) >= 0);
	if (!canR) return -1;

	general_assert(!(b->core.flag & 0x0001), "Find a paired end read in the file!");

	int readType = getReadType(b);
	std::string name = getName(b);

	if (readType != 1 || (readType == 1 && read.getName().compare(name) != 0)) {
		val = readType;
		read = SingleReadQ(name, getReadSeq(b), getQScore(b));
	}
	else val = 5;

	if (readType == 1) {
		if (!check(b)) { fprintf(stderr, "RSEM does not support gapped alignments, sorry!\n"); exit(-1); }

		if (getDir(b) > 0) {
			hit = SingleHit(transcripts.getInternalSid(b->core.tid + 1), b->core.pos);
		}
		else {
			hit = SingleHit(-transcripts.getInternalSid(b->core.tid + 1), header->target_len[b->core.tid] - b->core.pos - b->core.l_qseq);
		}
	}

	return val;
}

//Assume whether aligned or not , two mates of paired-end reads are always get together
int SamParser::parseNext(PairedEndRead& read, PairedEndHit& hit) {
	int val;
	bool canR = ((samread(sam_in, b) >= 0) && (samread(sam_in, b2) >= 0));
	if (!canR) return -1;

	general_assert((b->core.flag & 0x0001) && (b2->core.flag & 0x0001), \
			"One of the mate is not paired-end! (RSEM assumes the two mates of a paired-end read should be adjacent)");

	bam1_t *mp1 = NULL, *mp2 = NULL;

	// If lose mate info, discard. is it necessary?
	if (!((b->core.flag & 0x0040) && (b2->core.flag & 0x0080)) && !((b->core.flag & 0x0080) && (b2->core.flag & 0x0040))) return 4;
	// If only one mate is mapped, discard
	if (((b->core.flag & 0x0004) && !(b2->core.flag & 0x0004)) || (!(b->core.flag & 0x0004) && (b2->core.flag & 0x0004))) return 4;

	if (b->core.flag & 0x0040) {
		mp1 = b; mp2 = b2;
	}
	else  {
		mp1 = b2; mp2 = b;
	}

	general_assert(!strcmp(bam1_qname(mp1), bam1_qname(mp2)), "Detected a read pair whose two mates have different names: " + getName(mp1) + " , " + getName(mp2) + " !");

	int readType = getReadType(mp1, mp2);
	std::string name = getName(mp1);

	if (readType != 1 || (readType == 1 && read.getName().compare(name) != 0)) {
		val = readType;
		SingleRead mate1(getName(mp1), getReadSeq(mp1));
		SingleRead mate2(getName(mp2), getReadSeq(mp2));
		read = PairedEndRead(mate1, mate2);
	}
	else val = 5;

	if (readType == 1) {
		if (!check(mp1) || !check(mp2)) { fprintf(stderr, "RSEM does not support gapped alignments, sorry!\n"); exit(-1); }

		if (mp1->core.tid != mp2->core.tid) {
			fprintf(stderr, "The two reads do not come from the same pair!");
			exit(-1);
		}
		//assert(mp1->core.tid == mp2->core.tid);
		if (getDir(mp1) > 0) {
			hit = PairedEndHit(transcripts.getInternalSid(mp1->core.tid + 1), mp1->core.pos, mp2->core.pos + mp2->core.l_qseq - mp1->core.pos);
		}
		else {
			hit = PairedEndHit(-transcripts.getInternalSid(mp1->core.tid + 1), header->target_len[mp1->core.tid] - mp1->core.pos - mp1->core.l_qseq, mp1->core.pos + mp1->core.l_qseq - mp2->core.pos);
		}
	}

	return val;
}

int SamParser::parseNext(PairedEndReadQ& read, PairedEndHit& hit) {
	int val;
	bool canR = ((samread(sam_in, b) >= 0) && (samread(sam_in, b2) >= 0));
	if (!canR) return -1;

	general_assert((b->core.flag & 0x0001) && (b2->core.flag & 0x0001), \
			"One of the mate is not paired-end! (RSEM assumes the two mates of a paired-end read should be adjacent)");

	bam1_t *mp1 = NULL, *mp2 = NULL;

	// If lose mate info, discard. is it necessary?
	if (!((b->core.flag & 0x0040) && (b2->core.flag & 0x0080)) && !((b->core.flag & 0x0080) && (b2->core.flag & 0x0040))) return 4;
	// If only one mate is mapped, discard
	if (((b->core.flag & 0x0004) && !(b2->core.flag & 0x0004)) || (!(b->core.flag & 0x0004) && (b2->core.flag & 0x0004))) return 4;

	if (b->core.flag & 0x0040) {
		mp1 = b; mp2 = b2;
	}
	else  {
		mp1 = b2; mp2 = b;
	}

	general_assert(!strcmp(bam1_qname(mp1), bam1_qname(mp2)), "Detected a read pair whose two mates have different names: " + getName(mp1) + " , " + getName(mp2) + " !");

	int readType = getReadType(mp1, mp2);
	std::string name = getName(mp1);

	if (readType != 1 || (readType == 1 && read.getName().compare(name) != 0)) {
		val = readType;
		SingleReadQ mate1(getName(mp1), getReadSeq(mp1), getQScore(mp1));
		SingleReadQ mate2(getName(mp2), getReadSeq(mp2), getQScore(mp2));
		read = PairedEndReadQ(mate1, mate2);
	}
	else val = 5;

	if (readType == 1) {
		if (!check(mp1) || !check(mp2)) { fprintf(stderr, "RSEM does not support gapped alignments, sorry!\n"); exit(-1); }

		if (mp1->core.tid != mp2->core.tid) {
			fprintf(stderr, "The two reads do not come from the same pair!");
			exit(-1);
		}
		//assert(mp1->core.tid == mp2->core.tid);
		if (getDir(mp1) > 0) {
			hit = PairedEndHit(transcripts.getInternalSid(mp1->core.tid + 1), mp1->core.pos, mp2->core.pos + mp2->core.l_qseq - mp1->core.pos);
		}
		else {
			hit = PairedEndHit(-transcripts.getInternalSid(mp1->core.tid + 1), header->target_len[mp1->core.tid] - mp1->core.pos - mp1->core.l_qseq, mp1->core.pos + mp1->core.l_qseq - mp2->core.pos);
		}
	}

	return val;
}

inline std::string SamParser::getReadSeq(const bam1_t* b) {
	uint8_t *p = bam1_seq(b);
	std::string readseq = "";
	char base = 0;

	if (getDir(b) < 0) {
		for (int i = b->core.l_qseq - 1; i >= 0; i--) {
			switch(bam1_seqi(p, i)) {
			//case 0 : base = '='; break;
			case 1 : base = 'T'; break;
			case 2 : base = 'G'; break;
			case 4 : base = 'C'; break;
			case 8 : base = 'A'; break;
			case 15 : base = 'N'; break;
			default : assert(false);
			}
			readseq.append(1, base);
		}
	}
	else {
		for (int i = 0; i < b->core.l_qseq; i++) {
			switch(bam1_seqi(p, i)) {
			//case 0 : base = '='; break;
			case 1 : base = 'A'; break;
			case 2 : base = 'C'; break;
			case 4 : base = 'G'; break;
			case 8 : base = 'T'; break;
			case 15 : base = 'N'; break;
			default : assert(false);
			}
			readseq.append(1, base);
		}
	}

	return readseq;
}

inline std::string SamParser::getQScore(const bam1_t* b) {
	uint8_t *p = bam1_qual(b);
	std::string qscore = "";

	if (getDir(b) > 0) {
		for (int i = 0; i < b->core.l_qseq; i++) {
			qscore.append(1, (char)(*p + 33));
			++p;
		}
	}
	else {
		p = p + b->core.l_qseq - 1;
		for (int i = 0; i < b->core.l_qseq; i++) {
			qscore.append(1, (char)(*p + 33));
			--p;
		}
	}

	return qscore;
}

//0 ~ N0 , 1 ~ N1, 2 ~ N2
inline int SamParser::getReadType(const bam1_t* b) {
	if (!(b->core.flag & 0x0004)) return 1;

	if (!strcmp(rtTag, "")) return 0;

	uint8_t *p = bam_aux_get(b, rtTag);
	if (p == NULL) return 0;
	return (bam_aux2i(p) > 0 ? 2 : 0);
}

//For paired-end reads, do not print out type 2 reads
inline int SamParser::getReadType(const bam1_t* b, const bam1_t* b2) {
	if (!(b->core.flag & 0x0004) && !(b2->core.flag & 0x0004)) return 1;

	if (!strcmp(rtTag, "")) return 0;

	uint8_t *p = bam_aux_get(b, rtTag);
	if (p != NULL && bam_aux2i(p) > 0) return 2;

	p = bam_aux_get(b2, rtTag);
	if (p != NULL && bam_aux2i(p) > 0) return 2;

	return 0;
}

#endif /* SAMPARSER_H_ */
