#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cassert>
#include<string>
#include<set>

#include <stdint.h>
#include "htslib/sam.h"
#include "sam_utils.h"

#include "utils.h"
#include "my_assert.h"

using namespace std;

samFile *in;
bam_hdr_t *header;
bam1_t *b, *b2;

set<string> used;

bool isValid;

int main(int argc, char* argv[]) {
  if (argc != 2) {
    printf("Usage: rsem-sam-validator <input.sam/input.bam/input.cram>\n");
    exit(-1);
  }
  
  in = sam_open(argv[1], "r");
  general_assert(in != 0, "Cannot open input file!");
  header = sam_hdr_read(in);
  general_assert(header != 0, "Cannot load SAM header!");
  used.clear();
  b = bam_init1(); b2 = bam_init1();
  
  isValid = true;
  
  HIT_INT_TYPE cnt = 0;
  string cqname(""), qname;
  uint64_t creadlen = 0, readlen;
  bool cispaired = false, ispaired;
  
  printf("."); fflush(stdout);
  do {
    if (sam_read1(in, header, b) < 0) break;
    assert(b->core.l_qseq > 0);
    
    qname = bam_get_canonical_name(b);
    
    // if this is a paired-end read
    ispaired = bam_is_paired(b);
    if (ispaired) {
      
      isValid = (sam_read1(in, header, b2) >= 0) && (qname == bam_get_canonical_name(b2)) && bam_is_paired(b2);
      if (!isValid) { printf("\nOnly find one mate for paired-end read %s!\n", qname.c_str()); continue; }
      assert(b2->core.l_qseq > 0);
      isValid = !(bam_is_read1(b) && bam_is_read2(b)) && !(bam_is_read1(b2) && bam_is_read2(b2));
      if (!isValid) { printf("\nRead %s has more than 2 segments (e.g. tripled or more ended reads)!\n", qname.c_str()); continue; }
      
      if (!bam_is_read1(b)) { bam1_t *tmp = b; b = b2; b2 = tmp; }
      isValid = bam_is_read1(b) && bam_is_read2(b2);
      if (!isValid) { printf("\nThe two mates of paired-end read %s are not adjacent!\n", qname.c_str()); continue; }
      
      // both mates are mapped
      if (bam_is_mapped(b) && bam_is_mapped(b2)) {
	isValid = (b->core.tid == b2->core.tid) && (b->core.pos == b2->core.mpos) && (b2->core.pos == b->core.mpos);
	if (!isValid) { printf("\nOne of paired-end read %s's alignment does not have two mates adjacent to each other! If you're running covert-sam-for-rsem now, this might mean the read contains duplicate alignments.\n", qname.c_str()); continue; }
      }
      
      readlen = (uint64_t(b->core.l_qseq) << 32) + b2->core.l_qseq;
    }
    else readlen = b->core.l_qseq;
    
    if (cqname != qname) {
      isValid = used.find(qname) == used.end();
      if (!isValid) { printf("\nThe alignments of read %s are not grouped together!", qname.c_str()); continue; }
      used.insert(cqname);
      cqname = qname;
      creadlen = readlen;
      cispaired = ispaired;
    }
    else {
      assert(cqname != "");
      isValid = (creadlen == readlen);
      if (!isValid) { printf("\nRead %s have different read/mate lengths!\n", qname.c_str()); continue; }
      isValid = (cispaired == ispaired);
      if (!isValid) { printf("\nRead %s is detected as both single-end read and paired-end read!\n", qname.c_str()); continue; }
    }
    
    ++cnt;
    if (cnt % 1000000 == 0) { printf("."); fflush(stdout); }
    
  } while(isValid);
  
  bam_destroy1(b); bam_destroy1(b2);
  bam_hdr_destroy(header);
  sam_close(in);
  
  if (isValid) printf("\nThe input file is valid!\n");
  else printf("The input file is not valid!\n");
  
  return 0;
}
