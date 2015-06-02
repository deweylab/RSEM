__doc__="""

  pliu 20150511

  python module for all parameters, input arguments
"""

class Param:
  IDR_THRESHOLD = 0.05

  def __init__(self):
    self.argdict = None
    self.num_threads                = None
    self.chipseq_target_read_files  = None
    self.chipseq_control_read_files = None
    self.bowtie_bin_for_chipseq     = None
    self.samtools_bin               = None
    self.bedtools_bin_for_chipseq   = None
    self.quiet                      = None
    self.ref_name    = None
    self.sample_name = None
    self.stat_name   = None
    self.imd_name    = None

  def __str__(self):
    ss = [ "%-28s %s\n" % (key, val) for (key, val) in self.argdict.items()]
    return ''.join(ss)

  @classmethod
  def initFromCommandLineArguments(cls, argdict):
    prm = cls()
    prm.argdict = argdict
    for (key, val) in argdict.items():
      setattr(prm, key, val)
    return prm


def initFromCommandLineArguments(argdict):
  return Param.initFromCommandLineArguments(argdict)
