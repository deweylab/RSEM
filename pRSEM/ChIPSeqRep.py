__doc__="""

  pliu 20150510

  python module for one ChIP-seq replicate
"""

import File

class ChIPSeqRep:
  def __init__(self):
    self.param     = None  ## reference to parameters
    self.rep_index = None
    self.is_target = None
    self.ffq       = None  ## File object for fastq
    self.fta       = None  ## File object for tagAlign
    self.chipseqexp = None ## reference to ChIPSeqExp object
    self.encoding   = None ## fastq encoding, not sure if needed

  @classmethod
  def initFromFastqFile(cls, ffq):
    csr = cls()
    csr.ffq = File.initFromFullFileName(ffq)

    return csr


def initFromFastqFile(ffq):
  return ChIPSeqRep.initFromFastqFile(ffq)
