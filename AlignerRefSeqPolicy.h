#ifndef ALIGNERREFSEQPOLICY
#define ALIGNERREFSEQPOLICY

#include<string>

#include "RefSeqPolicy.h"

class AlignerRefSeqPolicy : public RefSeqPolicy {
 public :
  std::string convert(const std::string& rawseq) {
    int size = (int)rawseq.size();
    std::string seq = rawseq;
    for (int i = 0; i < size; i++)
      if (seq[i] == 'N') seq[i] = 'G';
    return seq;
  }
};

#endif
