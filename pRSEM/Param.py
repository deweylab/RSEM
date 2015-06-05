__doc__="""

  pliu 20150511

  python module for all parameters, input arguments
"""

class Param:
  IDR_THRESHOLD = 0.05

  def __init__(self):
    self.argdict = None

    ## has to be in the same naming convention as prsem-calculate-expression
    self.num_threads                = None
    self.chipseq_target_read_files  = None
    self.chipseq_control_read_files = None
    self.bowtie_bin_for_chipseq     = None
    self.samtools_bin               = None
    self.bedtools_bin_for_chipseq   = None
    self.quiet                      = None
    ##

    self.ref_name    = None
    self.sample_name = None
    self.stat_name   = None
    self.imd_name    = None

    self.temp_dir    = None     ## dir to save RSEM/pRSEM intermediate files
    self.pRSEM_scr_dir   = None ## pRSEM scripts dir
    self.chipseq_rscript = None ## fullname of prsem-process-chipseq.R


  def __str__(self):
    ss = [ "%-28s %s\n" % (key, val) for (key, val) in self.argdict.items()] + \
         [ "%-28s %s\n" % ('RSEM_temp_dir', self.temp_dir ) ] + \
         [ "%-28s %s\n" % ('pRSEM_scr_dir', self.pRSEM_scr_dir) ]
    return ''.join(ss)


  @classmethod
  def initFromCommandLineArguments(cls, argdict):
    import os
    prm = cls()
    prm.argdict = argdict
    for (key, val) in argdict.items():
      setattr(prm, key, val)
    prm.temp_dir = os.path.split(prm.imd_name)[0] + '/'
    prm.pRSEM_scr_dir = os.path.dirname(os.path.realpath(__file__)) + '/'
    prm.chipseq_rscript = prm.pRSEM_scr_dir + 'prsem-process-chipseq.R'
    return prm


def initFromCommandLineArguments(argdict):
  return Param.initFromCommandLineArguments(argdict)
