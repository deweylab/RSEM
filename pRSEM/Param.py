__doc__="""

  pliu 20150511

  python module for all parameters, input arguments
"""

class Param:
  IDR_THRESHOLD  = 0.05
  N_PEAK         = 300000
  PEAK_TYPE      = '-savr'
  EXCLUSION_ZONE = '-500:85'  ## Anshul recommend -500:85
  TRAINING_GENE_MIN_LEN   = 1003
  TRAINING_MIN_MAPPABILITY = 0.8
  FLANKING_WIDTH = 500  ## in nt, flanking region around TSS and TES
  INFORMATIVE_DATA_MAX_P_VALUE = 0.01 ## external data set is informative if
                                      ## p-value is not more than this value

  def __init__(self):
    self.argdict = None

    ## has to be in the same naming convention as prsem-calculate-expression
    self.num_threads                = None
    self.chipseq_target_read_files  = None
    self.chipseq_control_read_files = None
    self.chipseq_read_files_multi_targets = None
    self.chipseq_bed_files_multi_targets  = None
    self.cap_stacked_chipseq_reads   = None
    self.n_max_stacked_chipseq_reads = None
    self.bowtie_path                = None
    self.chipseq_peak_file          = None
    self.mappability_bigwig_file    = None
    self.partition_model            = None
    self.gibbs_burnin               = None
    self.gibbs_number_of_samples    = None
    self.gibbs_sampling_gap         = None
    self.quiet                      = False

    ## arguments
    self.ref_fasta   = None
    self.ref_name    = None
    self.sample_name = None
    self.stat_name   = None
    self.imd_name    = None

    ## path and pRSEM scripts
    self.temp_dir       = None ## dir to save RSEM/pRSEM intermediate files
    self.prsem_scr_dir  = None ## pRSEM scripts dir
    self.prsem_rlib_dir = None ## place to install pRSEM required R libraries

    ## genome reference: training set isoforms
    self.fall_exon_crd     = None
    self.fall_tr_crd       = None ## tr info + mappability
    self.ftraining_tr_crd  = None ## training set tr

    ## ChIP-seq
    self.chipseqexperiment_target  = None ## reference to ChIP-seq experiment
    self.chipseqexperiment_control = None ## reference to ChIP-seq experiment
    self.chipseq_rscript   = None ## full name of process-chipseq.R
    self.filterSam2Bed     = None ## full name of filterSam2Bed binary
    self.spp_tgz           = None
    self.spp_script        = None
    self.idr_scr_dir       = None
    self.idr_script        = None
    self.fgenome_table     = None
    self.fidr_chipseq_peaks = None
    self.fall_chipseq_peaks = None
    self.fchipseq_peaks    = None ## full name of user supplied ChIP-seq peak
                                  ## file, otherwise is fidr_chipseq_peaks
    self.chipseq_target_fraglen = None ## spp-estimated fragment length
    self.fsppout_target = None ## full name of SPP output
                               ## this implementation needs to be refined since
                               ## the var is define in both Param and ChIPSeqExp
    self.fchipseq_target_signals  = None
    self.fchipseq_control_signals = None

    ## transcripts and RNA-seq
    self.transcripts = None ## reference to all transcripts to be quantified
    self.genes       = None ## reference to all genes to be quantified

    self.rnaseq_rscript    = None ## fullname of R script for dealing RNA-seq
    self.fti               = None ## RSEM's reference .ti file
    self.bigwigsummary_bin = None ## bigWigSummary binary
    self.fall_tr_features  = None ## file for all isoforms' features
    self.fall_tr_prior     = None ## file for all isoforms' priors
    self.fisoforms_results = None ## file for RSEM .isoforms.results
    self.fpvalLL           = None ## file for p-value on if informative
                                  ## and for log-likelihood
    self.fall_pvalLL = None ## file to store all the p-val and log-likelihood

    ## for multiple external data sets
    self.targetid2fchipseq_alignment = {}
    self.finfo_multi_targets = None
    self.flgt_model_multi_targets = None

    ## for testing procedure
    self.targetids = []


  def __str__(self):
    ss = [ "%-33s %s\n" % (key, val) for (key, val) in self.argdict.items()] + \
         [ "%-33s %s\n" % ('RSEM_temp_dir', self.temp_dir ) ] + \
         [ "%-33s %s\n" % ('pRSEM_scr_dir', self.prsem_scr_dir) ]
    return ''.join(ss)


  @classmethod
  def initFromCommandLineArguments(cls, argdict):
    import os
    prm = cls()
    prm.argdict = argdict
    for (key, val) in argdict.items():
      setattr(prm, key, val)

    if prm.imd_name is not None:
      prm.temp_dir = os.path.split(prm.imd_name)[0] + '/'
    prm.prsem_scr_dir = os.path.dirname(os.path.realpath(__file__)) + '/'
    prm.prsem_rlib_dir = prm.prsem_scr_dir + 'RLib/'
    if not os.path.exists(prm.prsem_rlib_dir):
      os.mkdir(prm.prsem_rlib_dir)

    ## genome reference: pRSEM training set isoforms
    prm.fall_exon_crd    = prm.ref_name + '_prsem.all_exon_crd'
    prm.fall_tr_crd      = prm.ref_name + '_prsem.all_tr_crd'
    prm.ftraining_tr_crd = prm.ref_name + '_prsem.training_tr_crd'

    ## ChIP-seq
    prm.chipseq_rscript = prm.prsem_scr_dir + 'process-chipseq.R'
    prm.filterSam2Bed   = prm.prsem_scr_dir + 'filterSam2Bed'
    prm.spp_tgz = prm.prsem_scr_dir + 'phantompeakqualtools/spp_1.10.1.tar.gz'
    prm.spp_script    = prm.prsem_scr_dir + 'phantompeakqualtools/run_spp.R'
    prm.idr_scr_dir   = prm.prsem_scr_dir + 'idrCode/'
    prm.idr_script    = prm.idr_scr_dir + 'batch-consistency-analysis.r'
    prm.fgenome_table = prm.ref_name + '.chrlist'

    if prm.temp_dir is not None:
      prm.fsppout_target           = prm.temp_dir + 'target_phantom.tab'
      prm.fchipseq_target_signals  = prm.temp_dir + 'target.tagAlign.gz'
      prm.fchipseq_control_signals = prm.temp_dir + 'control.tagAlign.gz'
      prm.fidr_chipseq_peaks = "%s/%s" % (prm.temp_dir,
                                          'idr_target_vs_control.regionPeak.gz')
      ## have to name it this way due to run_spp.R's wired naming convention
      ## this names depens on the next two names
      prm.fall_chipseq_peaks = "%s/%s" % (prm.temp_dir,
                            'target.tagAlign_VS_control.tagAlign.regionPeak.gz')

    if prm.chipseq_peak_file is not None:
      prm.fchipseq_peaks = prm.chipseq_peak_file
    else:
      prm.fchipseq_peaks = prm.fidr_chipseq_peaks


    ## transcripts and RNA-seq
    prm.rnaseq_rscript = prm.prsem_scr_dir + 'process-rnaseq.R'
    prm.fti            = prm.ref_name + '.ti'
    prm.ffasta         = prm.ref_name + '.transcripts.fa'
    prm.bigwigsummary_bin = prm.prsem_scr_dir + 'bigWigSummary'
   #prm.fall_exon_crd     = prm.imd_name  + '_prsem.all_exon_crd'
   #prm.fall_tr_crd       = prm.imd_name  + '_prsem.all_tr_crd'
   #prm.ftraining_tr_crd  = prm.imd_name  + '_prsem.training_tr_crd'
    if prm.sample_name is not None: ## for calc-expr
      prm.fall_tr_gc        = prm.imd_name  + '_prsem.all_tr_gc'
      prm.fall_tr_features  = prm.stat_name + '_prsem.all_tr_features'
      prm.fall_tr_prior     = prm.stat_name + '_prsem.all_tr_prior'
      prm.fpvalLL           = prm.stat_name + '_prsem.pval_LL'

      prm.fisoforms_results = prm.sample_name + '.isoforms.results'
      prm.fall_pvalLL       = prm.sample_name + '.all.pval_LL'

      ## for multiple external data sets
      prm.finfo_multi_targets      = prm.temp_dir + 'multi_targets.info'
      prm.flgt_model_multi_targets = prm.stat_name + '_prsem.lgt_mdl.RData'

    return prm


def initFromCommandLineArguments(argdict):
  return Param.initFromCommandLineArguments(argdict)
