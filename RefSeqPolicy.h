#ifndef REFSEQPOLICY
#define REFSEQPOLICY

#include<string>

/**
Convert reference sequences to RSEM format
 */
class RefSeqPolicy {
 public:
  std::string convert(const std::string& rawseq) {
    int size = (int)rawseq.size();
    std::string seq = rawseq;
    for (int i = 0; i < size; i++) {
      seq[i] = toupper(rawseq[i]);
      if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T') seq[i] = 'N';
    }
    return seq;
  }
};

#endif
