#!/bin/env python

__doc__="""

  pliu 20150304

  python function for pRSEM
"""

import os
import sys
import Util


def genChIPSeqSignalFilesFromBed(param):
  import ChIPSeqReplicate
  fbeds = param.chipseq_bed_files_multi_targets.split(',')
  for fbed in fbeds:
    csr = ChIPSeqReplicate.initFromBedFile(fbed)
    ta = csr.tagalign
    param.targetid2fchipseq_alignment[ta.basename] = ta.fullname


def genChIPSeqSignalFilesFromReads(param):
  import ChIPSeqExperiment
  cse_target = ChIPSeqExperiment.initFromParam(param, 'multi-targets')
  cse_target.getFastqEncoding()
  cse_target.alignReadByBowtie()

  param.chipseqexperiment_target = cse_target
  for rep in cse_target.reps:
    ta = rep.tagalign
    param.targetid2fchipseq_alignment[ta.basename] = ta.fullname



def genChIPSeqPeakFileBySPPIDR(param):
  import ChIPSeqExperiment

  cse_target = ChIPSeqExperiment.initFromParam(param, 'target')
  cse_target.getFastqEncoding()
  cse_target.alignReadByBowtie()
  cse_target.poolTagAlign()

  param.chipseqexperiment_target = cse_target

  if param.chipseq_control_read_files is not None:
    cse_control = ChIPSeqExperiment.initFromParam(param, 'control')
    cse_control.getFastqEncoding()
    cse_control.alignReadByBowtie()
    cse_control.poolTagAlign()
    cse_target.callPeaksBySPP(cse_control.pooled_tagalign)
    cse_target.getPeaksByIDR(cse_control.pooled_tagalign)

    param.chipseq_peak_file = cse_target.final_peaks.fullname
    param.chipseqexperiment_control = cse_control
  else:
    pass  ## to-be-implemented, call peaks by MOSAiCS without ChIP-seq control


def buildTrainingSet(prm):
  """
  write training set in file Param.ftraining_tr_crd
  transcript as listed in the same order as RSEM's .ti file
  The order is required by rsem-run-gibbs so that prior can be assigned to
  transcript correctly
  """
  ogot_genes = filter(lambda g: len(g.transcripts) == 1 and
                                (g.end - g.start + 1) >=
                                prm.TRAINING_GENE_MIN_LEN, prm.genes)

  trs = [tr for g in ogot_genes for tr in g.transcripts]

  trid2mpps = Util.runMPOverAList(prm.num_threads, calTSSBodyTESMappability,
                                  [trs, prm])

  with open(prm.fall_tr_crd, 'w') as f_fout:
    f_fout.write("geneid\ttrid\tchrom\tstrand\tstart\tend\t")
    f_fout.write("tss_mpp\tbody_mpp\ttes_mpp\n")
    for tr in prm.transcripts: ## in the same order as RSEM's .ti file
      f_fout.write("%s\t%s\t%s\t%s\t%d\t%d\t" % ( tr.gene_id,
                   tr.transcript_id, tr.chrom, tr.strand, tr.start, tr.end))
      if tr.transcript_id in trid2mpps:
        mpps = trid2mpps[tr.transcript_id]
        f_fout.write("%5.3f\t%5.3f\t%5.3f\n" % mpps)
      else:
        f_fout.write("NA\tNA\tNA\n")

  with open(prm.fall_exon_crd, 'w') as f_fexon:
    f_fexon.write("trid\texon_index\tchrom\tstrand\tstart\tend\n")
    for tr in prm.transcripts:
      for (i, (exon_start, exon_end)) in enumerate(tr.exon_ranges):
        f_fexon.write("%s\t%d\t%s\t%s\t%d\t%d\n" % (tr.transcript_id, i+1,
                      tr.chrom, tr.strand, exon_start, exon_end))

  Util.runCommand('/bin/env', 'Rscript', prm.rnaseq_rscript, 'selTrainingTr',
                  prm.prsem_rlib_dir, prm.fall_tr_crd, prm.fall_exon_crd,
                  prm.TRAINING_MIN_MAPPABILITY, prm.FLANKING_WIDTH,
                  prm.ftraining_tr_crd, quiet=prm.quiet)

  if not os.path.exists(prm.ftraining_tr_crd):
    sys.exit("Failed to generate file: %s\n" % prm.ftraining_tr_crd)


def calTSSBodyTESMappability(trs, prm, out_q):
  """
  calculate average mappability around TSS, body, and TES for all transcripts of
  given list of genes

  save results in transcript's attribute
  """
  outdict = {}
  for tr in trs:
    tr.calculateMappability(prm.bigwigsummary_bin, prm.mappability_bigwig_file,
                            prm.FLANKING_WIDTH, prm.quiet)
    outdict[tr.transcript_id] = (tr.ave_mpp_around_TSS, tr.ave_mpp_around_body,
                                 tr.ave_mpp_around_TES)
  out_q.put(outdict)


def genPriorByCombinedTSSSignals(prm):
  """
  calculate TSS signals for all external data sets
  compute informative p-value, LL for individual data set and combined one
  learn prior from training set partitioned by combined TSS signals
  derive priors for all isoforms
  """
  f_fout = open(prm.finfo_multi_targets, 'w')
  f_fout.write("targetid\tfaln\tfftrs\n")
  for (tgtid, faln) in prm.targetid2fchipseq_alignment.items():
    fftrs = prm.imd_name + '_prsem.' + tgtid + '.all_tr_features'
    f_fout.write("%s\t%s\t%s\n" % (tgtid, faln, fftrs))
  f_fout.close()

  Util.runCommand('/bin/env', 'Rscript', prm.rnaseq_rscript,
                  'prepMultiTargetsFeatures', prm.prsem_rlib_dir,
                  prm.fall_tr_crd, prm.ftraining_tr_crd,
                  prm.fisoforms_results, prm.FLANKING_WIDTH,
                  prm.cap_stacked_chipseq_reads,
                  prm.n_max_stacked_chipseq_reads,
                  prm.finfo_multi_targets, prm.num_threads, quiet=prm.quiet)

  ## learn prior from partitioning by combined external data set
  Util.runCommand('/bin/env', 'Rscript', prm.rnaseq_rscript,
                  'genPriorByCombinedTSSSignals', prm.prsem_rlib_dir,
                  prm.finfo_multi_targets, prm.flgt_model_multi_targets,
                  prm.fall_tr_features, prm.fpvalLL, prm.fall_tr_prior,
                  quiet=prm.quiet)

  pval = float(Util.readFile(prm.fpvalLL)[1].split("\t")[0])

  if pval > prm.INFORMATIVE_DATA_MAX_P_VALUE:
    err_msg = "\nError: current external data is NOT informative for RNA-seq quantification\n" + \
      "\tp-value %.10e > %.3f\n" % (pval, prm.INFORMATIVE_DATA_MAX_P_VALUE) + \
      "pRSEM STOPs here. Please use other external data set(s)\n\n"
    sys.stderr.write(err_msg)
    sys.exit(0)

  if not os.path.exists(prm.fall_tr_prior):
    sys.exit("Failed to generate file: %s\n" % prm.fall_tr_prior)



def genPriorByPeakSignalGCLen(prm):
  """
  calculate peaks/signals for the TSS, body, and TES regions
  calculate GC contenct and effective length
  learn prior from training set and derived priors for all isoforms
  """
  ## calculate GC contect for isoforms
  trid2seq = Util.getFastaID2Seq(prm.ffasta)
  with open(prm.fall_tr_gc, 'w') as f_fall_tr_gc:
    f_fall_tr_gc.write("trid\tGC_fraction\n")
    for tr in prm.transcripts:
      gc_frac = Util.getGCFraction(trid2seq[tr.transcript_id])
      f_fall_tr_gc.write("%s\t%.2f\n" % (tr.transcript_id, gc_frac) )

  with open(prm.fsppout_target, 'r') as f_fsppout_target:
    words = f_fsppout_target.read().split("\t")
    prm.chipseq_target_fraglen = int(words[2])

  ## prepare a feature file of peaks and signals for all isoforms,
  ## isoforms in training set will be labeled
  if not os.path.exists(prm.fchipseq_peaks):
    sys.exit("File not exists: %s\n" % prm.fchipseq_peaks)
  Util.runCommand('/bin/env', 'Rscript', prm.rnaseq_rscript,
                  'prepPeakSignalGCLenFeatures', prm.prsem_rlib_dir,
                  prm.fall_tr_crd, prm.ftraining_tr_crd, prm.fall_tr_features,
                  prm.fisoforms_results, prm.FLANKING_WIDTH,
                  prm.partition_model, prm.fchipseq_peaks,
                  prm.fchipseq_target_signals, prm.fall_tr_gc, prm.num_threads,
                  prm.chipseq_target_fraglen, quiet=prm.quiet)

  if not os.path.exists(prm.fall_tr_gc):
    sys.exit("Failed to generate file: %s\n" % prm.fall_tr_gc)

  ## learn and generate prior for all transcripts
  Util.runCommand('/bin/env', 'Rscript', prm.rnaseq_rscript,
                  'genPriorByPeakSignalGCLen', prm.prsem_rlib_dir,
                  prm.fall_tr_features, prm.partition_model, prm.fall_tr_prior,
                  quiet=prm.quiet)

  if not os.path.exists(prm.fall_tr_prior):
    sys.exit("Failed to generate file: %s\n" % prm.fall_tr_prior)


def genPriorByTSSPeak(prm):
  """
  determine if isoform have TSS peak or not
  learn priors from training set and derived priors for all isoforms
  """
  ## prepare a feature file of TSS peaks for all isoforms,
  ## isoforms in training set will be labeled
  if not os.path.exists(prm.fchipseq_peaks):
    sys.exit("File not exists: %s\n" % prm.fchipseq_peaks)
  Util.runCommand('/bin/env', 'Rscript', prm.rnaseq_rscript,
                  'prepTSSPeakFeatures', prm.prsem_rlib_dir,
                  prm.fall_tr_crd, prm.ftraining_tr_crd, prm.fall_tr_features,
                  prm.fisoforms_results, prm.FLANKING_WIDTH,
                  prm.fchipseq_peaks, quiet=prm.quiet)

  if not os.path.exists(prm.fall_tr_features):
    sys.exit("Failed to generate file: %s\n" % prm.fall_tr_features)

  Util.runCommand('/bin/env', 'Rscript', prm.rnaseq_rscript,
                  'genPriorByTSSPeak', prm.prsem_rlib_dir,
                  prm.fall_tr_features,  prm.fpvalLL, prm.fall_tr_prior,
                  quiet=prm.quiet)

  pval = float(Util.readFile(prm.fpvalLL)[1].split("\t")[0])

  if pval > prm.INFORMATIVE_DATA_MAX_P_VALUE:
    err_msg = "\nError: current external data is NOT informative for RNA-seq quantification\n" + \
      "\tp-value %.10e > %.3f\n" % (pval, prm.INFORMATIVE_DATA_MAX_P_VALUE) + \
      "pRSEM STOPs here. Please use other external data set(s)\n\n"
    sys.stderr.write(err_msg)
    sys.exit(0)

  if not os.path.exists(prm.fall_tr_prior):
    sys.exit("Failed to generate file: %s\n" % prm.fall_tr_prior)


def runGibbsSampling(prm):
  if prm.quiet:
    run_gibbs_quiet = '-q'
  else:
    run_gibbs_quiet = ''
  Util.runCommand("%s/../rsem-run-gibbs" % prm.prsem_scr_dir,
                  prm.ref_name, prm.imd_name, prm.stat_name, prm.gibbs_burnin,
                  prm.gibbs_number_of_samples, prm.gibbs_sampling_gap,
                  '-p', prm.num_threads, run_gibbs_quiet,
                  '--prior', prm.fall_tr_prior,
                  quiet=prm.quiet)
