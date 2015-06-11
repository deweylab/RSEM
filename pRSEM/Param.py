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
  TRAINING_GENE_MIN_LEN   = 1003
  TRAINING_MIN_MAPPABILITY = 0.8
  FLANKING_WIDTH = 500  ## in nt, flanking region around TSS and TES

  def __init__(self):
    self.argdict = None

    ## has to be in the same naming convention as prsem-calculate-expression
    self.num_threads                = None
    self.chipseq_target_read_files  = None
    self.chipseq_control_read_files = None
    self.bowtie_bin_for_chipseq     = None
    self.samtools_bin               = None
    self.bedtools_bin_for_chipseq   = None
    self.chipseq_peak_file          = None
    self.mappability_bigwig_file    = None
    self.partition_model            = None
    self.gibbs_burnin               = None
    self.gibbs_number_of_samples    = None
    self.gibbs_sampling_gap         = None
    self.quiet                      = None

    ## arguments
    self.ref_name    = None
    self.sample_name = None
    self.stat_name   = None
    self.imd_name    = None

    ## path and pRSEM scripts
    self.temp_dir       = None ## dir to save RSEM/pRSEM intermediate files
    self.prsem_scr_dir  = None ## pRSEM scripts dir
    self.prsem_rlib_dir = None ## place to install pRSEM required R libraries

    ## ChIP-seq
    self.chipseqexperiment_target  = None ## reference to ChIP-seq experiment
    self.chipseqexperiment_control = None ## reference to ChIP-seq experiment
    self.chipseq_rscript   = None ## fullname of process-chipseq.R
    self.spp_tgz           = None
    self.spp_script        = None
    self.idr_scr_dir       = None
    self.idr_script        = None
    self.fgenome_table     = None
    self.fidr_chipseq_peaks = None
    self.fall_chipseq_peaks = None
    self.fchipseq_target_signals  = None
    self.fchipseq_control_signals = None

    ## transcripts and RNA-seq
    self.transcripts = None ## reference to all transcripts to be quantified
    self.genes       = None ## reference to all genes to be quantified

    self.rnaseq_rscript     = None ## fullname of R script for dealing RNA-seq
    self.fti                = None ## RSEM's reference .ti file
    self.bigwigsummary_bin  = None ## bigWigSummary binary
    self.fall_tr_crd        = None ## tr info + mappability
    self.ftraining_tr_crd   = None ## training set tr
    self.fall_tr_features   = None ## file for all isoforms' features
    self.fall_tr_prior      = None ## file for all isoforms' priors
    self.fisoforms_results  = None ## file for RSEM .isoforms.results


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
    prm.prsem_rlib_dir = prm.prsem_scr_dir + 'Rlib/'
    if not os.path.exists(prm.prsem_rlib_dir):
      os.mkdir(prm.prsem_rlib_dir)

    ## ChIP-seq
    prm.chipseq_rscript = prm.prsem_scr_dir + 'process-chipseq.R'
    prm.spp_tgz = prm.prsem_scr_dir + 'phantompeakqualtools/spp_1.10.1.tar.gz'
    prm.spp_script    = prm.prsem_scr_dir + 'phantompeakqualtools/run_spp.R'
    prm.idr_scr_dir   = prm.prsem_scr_dir + 'idrCode/'
    prm.idr_script    = prm.idr_scr_dir + 'batch-consistency-analysis.r'
    prm.fgenome_table = prm.ref_name + '.chrlist'
    prm.fidr_chipseq_peaks = "%s/%s" % (prm.temp_dir,
                                        'idr_target_vs_control.regionPeak.gz')
    ## have to name it this way due to run_spp.R's wired naming convention
    ## this names depens on the next two names
    prm.fall_chipseq_peaks = "%s/%s" % (prm.temp_dir,
                            'target.tagAlign_VS_control.tagAlign.regionPeak.gz')
    prm.fchipseq_target_signals  = prm.temp_dir + 'target.tagAlign.gz'
    prm.fchipseq_control_signals = prm.temp_dir + 'control.tagAlign.gz'

    ## transcripts and RNA-seq
    prm.rnaseq_rscript = prm.prsem_scr_dir + 'process-rnaseq.R'
    prm.fti            = prm.ref_name + '.ti'
    prm.ffasta         = prm.ref_name + '.transcripts.fa'
    prm.bigwigsummary_bin  = prm.prsem_scr_dir + 'bigWigSummary'
    prm.fall_tr_crd        = prm.imd_name  + '_prsem.all_tr_crd'
    prm.fall_tr_gc         = prm.imd_name  + '_prsem.all_tr_gc'
    prm.ftraining_tr_crd   = prm.imd_name  + '_prsem.training_tr_crd'
    prm.fall_tr_features   = prm.stat_name + '_prsem.all_tr_features'
    prm.fall_tr_prior      = prm.stat_name + '_prsem.all_tr_prior'
    prm.fisoforms_results  = prm.sample_name + '.isoforms.results'

    return prm


def initFromCommandLineArguments(argdict):
  return Param.initFromCommandLineArguments(argdict)