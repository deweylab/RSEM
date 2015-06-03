__doc__="""

  pliu 20150511

  python module for a ChIP-seq experiment that contains
  replicates of ChIP-seq data for target and/or control
"""

import ChIPSeqReplicate


class ChIPSeqExperiment:
  def __init__(self):
    self.param = None ## reference to input parameters
    self.target_reps  = []
    self.control_reps = []
    self.fta_target_pool  = None  ## pooled tagAlign File object for target
    self.fta_control_pool = None  ## pooled tagAlign File object for control

  @classmethod
  def initFromParam(cls, param):
    cse = cls()
    cse.param = param
    ftgts = param.chipseq_target_read_files.split(',')
    cse.target_reps = [ ChIPSeqRep.initFromFastqFile(ffq) for ffq in ftgts ]
    for rep in cse.target_reps:
      rep.chipseqexp = cse

    if param.chipseq_control_read_files is not None:
      fctrls = param.chipseq_control_read_files.split(',')
      cse.control_reps = [ ChIPSeqRep.initFromFastqFile(ffq) for ffq in fctrls ]
      for rep in cse.control_reps:
        rep.chipseqexp = cse

    return cse


def initFromParam(param):
  return ChIPSeqExperiment.initFromParam(param)
