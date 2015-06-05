__doc__="""

  pliu 20150511

  python module for all parameters, input arguments
"""

class Param:
  IDR_THRESHOLD  = 0.05
  N_PEAK         = 300000
  PEAK_TYPE      = '-savr'
  EXCLUSION_ZONE = '-500:80'  ## Anshul recommend -500:85, but no results for
                              ## GM12878, have to loose the criterion

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
    self.prsem_scr_dir   = None ## pRSEM scripts dir
    self.chipseq_rscript = None ## fullname of prsem-process-chipseq.R

    ## for peak calling and IDR calculation
    self.spp_script    = None
    self.idr_script    = None
    self.fgenome_table = None


  def __str__(self):
    ss = [ "%-28s %s\n" % (key, val) for (key, val) in self.argdict.items()] + \
         [ "%-28s %s\n" % ('RSEM_temp_dir', self.temp_dir ) ] + \
         [ "%-28s %s\n" % ('pRSEM_scr_dir', self.prsem_scr_dir) ]
    return ''.join(ss)


  @classmethod
  def initFromCommandLineArguments(cls, argdict):
    import os
    prm = cls()
    prm.argdict = argdict
    for (key, val) in argdict.items():
      setattr(prm, key, val)
    prm.temp_dir = os.path.split(prm.imd_name)[0] + '/'
    prm.prsem_scr_dir = os.path.dirname(os.path.realpath(__file__)) + '/'
    prm.chipseq_rscript = prm.prsem_scr_dir + 'prsem-process-chipseq.R'
    prm.spp_script = prm.prsem_scr_dir + 'phantompeakqualtools/run_spp.R'
    prm.idr_script = prm.prsem_scr_dir + 'idrCode/batch-consistency-analysis.r'
    prm.fgenome_table = prm.ref_name + '.chrlist'
    return prm


def initFromCommandLineArguments(argdict):
  return Param.initFromCommandLineArguments(argdict)
