#ifndef MAQMAP_H_
#define MAQMAP_H_

#ifdef MAQ_LONGREADS
#  define MAX_READLEN 128
#else
#  define MAX_READLEN 64
#endif

#define MAX_NAMELEN 36
#define MAQMAP_FORMAT_OLD 0
#define MAQMAP_FORMAT_NEW -1

#define PAIRFLAG_FF      0x01
#define PAIRFLAG_FR      0x02
#define PAIRFLAG_RF      0x04
#define PAIRFLAG_RR      0x08
#define PAIRFLAG_PAIRED  0x10
#define PAIRFLAG_DIFFCHR 0x20
#define PAIRFLAG_NOMATCH 0x40
#define PAIRFLAG_SW      0x80

#include <string.h>
#include <zlib.h>
#include "const.h"

/*
  name: read name
  size: the length of the read
  seq: read sequence (see also below)
  seq[MAX_READLEN-1]: single end mapping quality (equals to map_qual if not paired)
  map_qual: the final mapping quality
  alt_qual: the lower quality of the two ends (equals to map_qual if not paired)
  flag: status of the pair
  dist: offset of the mate (zero if not paired)
  info1: mismatches in the 24bp (higher 4 bits) and mismatches (lower 4 bits)
  info2: sum of errors of the best hit
  c[2]: count of all 0- and 1-mismatch hits on the reference
 */
typedef struct
{
	bit8_t seq[MAX_READLEN]; /* the last base is the single-end mapping quality. */
	bit8_t size, map_qual, info1, info2, c[2], flag, alt_qual;
	bit32_t seqid, pos;
	int dist;
	char name[MAX_NAMELEN];
} maqmap1_t;

typedef struct
{
	int format, n_ref;
	char **ref_name;
	bit64_t n_mapped_reads;
	maqmap1_t *mapped_reads;
} maqmap_t;

#define maqmap_read1(fp, m1) gzread((fp), (m1), sizeof(maqmap1_t))

#ifdef __cplusplus
extern "C" {
#endif
	maqmap_t *maq_new_maqmap();
	void maq_delete_maqmap(maqmap_t *mm);
	void maqmap_write_header(gzFile fp, const maqmap_t *mm);
	maqmap_t *maqmap_read_header(gzFile fp);
#ifdef __cplusplus
}
#endif

#endif
