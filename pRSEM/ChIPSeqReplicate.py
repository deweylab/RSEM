__doc__="""

  pliu 20150510

  python module for one ChIP-seq replicate
"""

import File

class ChIPSeqReplicate:
  def __init__(self):
    self.fastq      = None  ## File object for fastq
    self.name       = None  ## default: fastq file's basename
    self.index      = None  ## replicate's index number
    self.tagalign   = None  ## File object for tagAlign
    self.encoding   = None  ## fastq encoding, not sure if needed

    self.param      = None  ## reference to parameters
    self.chipseqexp = None  ## reference to ChIPSeqExperiment object

 #def __str__(self):
 #  return "%s %s %d %s" % (self.fastq.fullname, self.name, self.index,
 #                          self.encoding)

  @classmethod
  def initFromFastqFile(cls, ffq):
    csr = cls()
    csr.fastq = File.initFromFullFileName(ffq)
    csr.name = csr.fastq.basename
    return csr

  @classmethod
  def initFromBedFile(cls, fbed):
    csr = cls()
    csr.tagalign = File.initFromFullFileName(fbed)
    csr.name = csr.tagalign.basename
    return csr

def initFromFastqFile(ffq):
  return ChIPSeqReplicate.initFromFastqFile(ffq)

def initFromBedFile(fbed):
  return ChIPSeqReplicate.initFromBedFile(fbed)
