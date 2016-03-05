#
#  pliu 20150608
#
#  module for processing RNA-seq data
#

main <- function() {
  name2func <- list(
    'selTrainingTr'                = selTrainingTr,
    'prepTSSPeakFeatures'          = prepTSSPeakFeatures,
    'prepPeakSignalGCLenFeatures'  = prepPeakSignalGCLenFeatures,
    'prepMultiTargetsFeatures'     = prepMultiTargetsFeatures,
    'genPriorByTSSPeak'            = genPriorByTSSPeak,
    'genPriorByPeakSignalGCLen'    = genPriorByPeakSignalGCLen,
    'genPriorByCombinedTSSSignals' = genPriorByCombinedTSSSignals
  )

  argv <- commandArgs(trailingOnly=T)
  name2func[[argv[1]]](argv[2:length(argv)])
}


genPriorByCombinedTSSSignals <- function(argv=NA){
  libloc <- argv[1]
  finfo  <- argv[2]
  fout_glmmdl <- argv[3]
  fout_ftrs   <- argv[4]
  fout_pvalLL <- argv[5]
  fout_prior  <- argv[6]

# libloc <- '/ua/pliu/dev/RSEM/pRSEM/RLib/'
# finfo  <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.temp/multi_targets.info'
# fout_glmmdl <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.stat/example_prsem.lgt_mdl.RData'
# fout_ftrs   <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.stat/example_prsem.all_tr_features'
# fout_pvalLL <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.stat/example_prsem.pval_LL'
# fout_prior  <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.stat/example_prsem.all_tr_prior'

  .libPaths(c(libloc, .libPaths()))
  suppressMessages(library(data.table))

  infodt <- fread(finfo, header=T, sep="\t")
  tgtids <- infodt[, targetid]
  ftrsdt <- rbindlist(lapply(tgtids,
              function(tgtid) {
                fin <- subset(infodt, targetid==tgtid)[, fftrs]
                dt <- fread(fin, header=T, sep="\t")
                dt[, `:=`( log10_tss_sig = ifelse(tss_sig > 0, log10(tss_sig),
                                                  -4),
                           targetid      = tgtid )]
                dt[, `:=`(tss_sig = NULL, nrd = NULL) ]
                return(dt)
              }))

  alldt <- dcast(ftrsdt, ... ~ targetid, value.var = 'log10_tss_sig')
  trndt <- subset(alldt, is_training == 1)
  frm <- paste0('is_expr ~ ', paste(sort(tgtids), collapse=' + '))

  glmmdl <- glm(frm, family='binomial', data=trndt)
  save(glmmdl, file=fout_glmmdl)

  alldt[, prd_expr_prob := predict(glmmdl, alldt, type='response')]
  alldt[, partition := factor(ifelse(prd_expr_prob > 0.5, 1, 0))]

  trn_prtdt <- subset(alldt, is_training == 1)
  fit <- getFitByMLDM(trn_prtdt[, pme_count], trn_prtdt[, partition])
  alldt[, prior:= fit$par[partition]]

  orig_ordered_trids <- subset(ftrsdt, targetid == tgtids[1])[, trid]
  setkey(alldt, trid)

  write.table(alldt[orig_ordered_trids], fout_ftrs, quote=F, sep='\t',
              col.names=T, row.names=F)

  not_expr_cnt <- subset(trn_prtdt, partition == 0)[, pme_count]
  expr_cnt     <- subset(trn_prtdt, partition == 1)[, pme_count]
  wrs <- suppressWarnings(wilcox.test(expr_cnt, not_expr_cnt,
                                      alternative='greater', paired=F, exact=T))
  pval <- wrs$p.value
  loglikelihood <- fit$value
  pvalLLdt <- data.table(pvalue = pval, loglikelihood = loglikelihood)
  write.table(pvalLLdt, fout_pvalLL, quote=F, sep="\t", col.names=T,
              row.names=F)

  out_priordt <- alldt[orig_ordered_trids, list(prior, trid)]
  write.table(out_priordt, fout_prior, quote=F, sep=' # ', col.names=F,
              row.names=F)
}


genPriorByPeakSignalGCLen <- function(argv=NA) {
  libloc           <- argv[1]
  fall_tr_features <- argv[2]
  partition_model  <- argv[3]
  fout             <- argv[4]

# libloc <- '/ua/pliu/dev/RSEM/pRSEM/RLib/'
# fall_tr_features <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.all_tr_features'
# partition_model <- 'pk_lgtnopk'
# fout <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.all_tr_prior'

  partition_model2func <- list(
    ## not needed in this function
   #'pk'   = getSampleAndPriorByPeak
   #'1prt' = getSampleAndPriorByOnePartition,

    ## w/ peak + logit on no peak
    'pk_lgtnopk' = getSampleAndPriorByPeakLogitNoPeak,

    ## linear regression on all with different number of bins
    'lm3'        = getSampleAndPriorByLM3,  # 3 bins
    'lm4'        = getSampleAndPriorByLM4,  # 4 bins
    'lm5'        = getSampleAndPriorByLM5,  # 5 bins
    'lm6'        = getSampleAndPriorByLM6,  # 6 bins

	  ## no peak + lm on w/ peak with different number of bins
    'nopk_lm2pk' = getSampleAndPriorByNoPeakLM2Peak,
    'nopk_lm3pk' = getSampleAndPriorByNoPeakLM3Peak,
    'nopk_lm4pk' = getSampleAndPriorByNoPeakLM4Peak,
    'nopk_lm5pk' = getSampleAndPriorByNoPeakLM5Peak,

    ## w/ peak + lm on no peak with different number of bins
    'pk_lm2nopk' = getSampleAndPriorByPeakLM2NoPeak,
    'pk_lm3nopk' = getSampleAndPriorByPeakLM3NoPeak,
    'pk_lm4nopk' = getSampleAndPriorByPeakLM4NoPeak,
    'pk_lm5nopk' = getSampleAndPriorByPeakLM5NoPeak
  )

  .libPaths(c(libloc, .libPaths()))
  suppressMessages(library(data.table))

  all_trdt <- fread(fall_tr_features, header=T, sep="\t")
  GC_mean <- mean(all_trdt[, GC_fraction])
  all_trdt[, `:=`( log10_count    = log10(pme_count + 1),
                   log10_tss_sig  = ifelse(tss_sig > 0,  log10(tss_sig),  -4.0),
                   log10_body_sig = ifelse(body_sig > 0, log10(body_sig), -4.0),
                   log10_tes_sig  = ifelse(tes_sig > 0,  log10(tes_sig),  -4.0),
                   log10_eff_len  = ifelse(efflen > 0,   log10(efflen),   -4.0),
                   log10_GC_ov_mean = ifelse(GC_fraction > 0,
                                             log10(GC_fraction/GC_mean),  -4.0),
                   no_tss_pk      = 1 - tss_pk,
                   no_body_pk     = 1 - body_pk,
                   no_tes_pk      = 1 - tes_pk
                 )]

  training_trdt <- subset(all_trdt, is_training==1)

  func <- partition_model2func[[partition_model]]

  outdt <- func(training_trdt, all_trdt)
  write.table(outdt[, list(prior, trid)], fout, quote=F, sep=' # ', col.names=F,
              row.names=F)
}


prepMultiTargetsFeatures <- function(argv=NA){
  libloc            <- argv[1]
  fall_tr_crd       <- argv[2]
  ftraining_tr_crd  <- argv[3]
  fisoforms_results <- argv[4]
  flanking_width    <- as.numeric(argv[5])
  cap_stacked_chipseq_reads   <- argv[6]
  n_max_stacked_chipseq_reads <- argv[7]
  finfo             <- argv[8]
  nthr              <- argv[9]

# libloc            <- '/ua/pliu/dev/RSEM/pRSEM/RLib/'
# fall_tr_crd       <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.temp/example_prsem.all_tr_crd'
# ftraining_tr_crd  <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.temp/example_prsem.training_tr_crd'
# fisoforms_results <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.isoforms.results'
# flanking_width    <- 500
# cap_stacked_chipseq_reads   <- 'True'
# n_max_stacked_chipseq_reads <- 5
# finfo             <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.temp/multi_targets.info'
# nthr              <- 16

# tmpdir <- '/tier2/deweylab/scratch/pliu/test/pRSEM/histone/03_rsem_expr/LTHSCRep1/LTHSCRep1.temp/'
# libloc            <- '/ua/pliu/dev/RSEM/pRSEM/RLib/'
# fall_tr_crd       <- paste0(tmpdir, 'LTHSCRep1_prsem.all_tr_crd')
# ftraining_tr_crd  <- paste0(tmpdir, 'LTHSCRep1_prsem.training_tr_crd')
# fisoforms_results <- paste0(tmpdir, '../LTHSCRep1.isoforms.results')
# flanking_width    <- 500
# cap_stacked_chipseq_reads   <- 'False'
# n_max_stacked_chipseq_reads <- 5
# finfo             <- paste0(tmpdir, 'multi_targets.info')
# nthr              <- 16

  .libPaths(c(libloc, .libPaths()))
  suppressMessages(library(data.table))
  suppressMessages(library(GenomicRanges))

  all_trdt <- fread(fall_tr_crd, header=T, sep="\t")
  training_trdt <- fread(ftraining_tr_crd, header=T, sep="\t", select='trid')

  rsemdt <- fread(fisoforms_results, header=T, sep="\t", select=c(
                  'transcript_id', 'posterior_mean_count', 'pme_TPM'))
  setnames(rsemdt, 1:2, c('trid', 'pme_count'))

  trdt <- merge(all_trdt, rsemdt, by='trid', all.x=T)
  trdt[, `:=`( tss         = ifelse(strand == '+', start, end),
               is_training = ifelse(trid %in% training_trdt[, trid], 1, 0) )]

  tssdt  <- trdt[, list(chrom, tss, trid)]
  tssdt[,  `:=`( start = tss - flanking_width,
                 end   = tss + flanking_width ) ]

  infodt <- fread(finfo, header=T, sep="\t")
  dum <- mclapply(infodt[, targetid], prepTSSSignalsFeatures, tssdt,
                  infodt, trdt, all_trdt, flanking_width,
                  cap_stacked_chipseq_reads, n_max_stacked_chipseq_reads,
                  mc.cores=nthr)
 #dum <- lapply(infodt[, targetid], prepTSSSignalsFeatures, tssdt,
 #              infodt, trdt, all_trdt, flanking_width,
 #              cap_stacked_chipseq_reads, n_max_stacked_chipseq_reads)
}


prepTSSSignalsFeatures <- function(tgtid, tssdt, infodt, trdt, all_trdt,
                                   flanking_width, is_cap, n_max_cap) {
  faln <- subset(infodt, targetid == tgtid)[, faln]
  fout <- subset(infodt, targetid == tgtid)[, fftrs]
  allrddt <- fread(paste0('zcat ',faln), header=F, sep="\t",
                   select=c('V1', 'V2', 'V3', 'V6'))
  setnames(allrddt, 1:4, c('chrom', 'start', 'end', 'strand'))

  if ( is_cap == 'True' ) {
    ## keep at most 5 reads per strand-specific interval
    allrddt[, dupi := seq_len(.N), by=list(chrom, start, end, strand)]
    rddt <- subset(allrddt, dupi <= n_max_cap)
    rddt[, dupi := NULL]
  } else {
    rddt <- allrddt
  }

  ## since no peak is called here, just use the average read length as fraglen
  ## count # of reads as signals rather than # of overlapping nucleotide
  ## normalize by TSS interval length and read depth to RPKM
  tssgrs <- makeGRangesFromDataFrame(tssdt, keep.extra.columns=T)
  rdgrs  <- makeGRangesFromDataFrame(rddt,  keep.extra.columns=T)

  ol <- findOverlaps(rdgrs, tssgrs, type='within', ignore.strand=T)
  oldt <- data.table(query=queryHits(ol), subject=subjectHits(ol))
  oldt[, trid := tssdt[, trid][subject]]
  nrddt <- oldt[, list(nrd = .N), by=trid]

  trdt <- merge(trdt, nrddt, by='trid', all.x=T)

  n_tot_rds <- length(rdgrs)
  trdt[, `:=`( tss_sig = ifelse(is.na(nrd), 0,
                                nrd*1e+9/(flanking_width*2+1)/n_tot_rds),
               is_expr = ifelse(pme_count > 0 & pme_TPM >= 1, 1, 0) )]

  setkey(trdt, trid)
  trdt <- trdt[all_trdt[, trid]] ## keep the order of original trids
  write.table(trdt, fout, quote=F, sep="\t", col.names=T, row.names=F)
}


prepPeakSignalGCLenFeatures <- function(argv=NA){
  libloc            <- argv[1]
  fall_tr_crd       <- argv[2]
  ftraining_tr_crd  <- argv[3]
  fout              <- argv[4]
  fisoforms_results <- argv[5]
  flanking_width    <- as.numeric(argv[6])
  partition_model   <- argv[7]
  fchipseq_peaks    <- argv[8]
  fchipseq_target_signals <- argv[9]
  fall_tr_gc        <- argv[10]
  nthr              <- argv[11]
  fraglen           <- as.numeric(argv[12])

# libloc            <- '/ua/pliu/dev/RSEM/pRSEM/RLib/'
# fall_tr_crd       <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/test.temp/test_prsem.all_tr_crd'
# ftraining_tr_crd  <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/test.temp/test_prsem.training_tr_crd'
# fout              <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/test.temp/test_prsem.all_tr_features'
# fisoforms_results <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/test.isoforms.results'
# flanking_width    <- 500
# partition_model   <- 'lm4'
# fchipseq_peaks    <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/test.temp/idr_target_vs_control.regionPeak.gz'
# fchipseq_target_signals <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/test.temp/target.tagAlign.gz'
# fall_tr_gc        <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/test.temp/test_prsem.all_tr_gc'
# nthr              <- 20
# fraglen           <- 110

# runid  <- 'lm4_k562rep1'
# exprdir <- paste0('/tier2/deweylab/scratch/pliu/test/pRSEM/', runid,
#                   '/rsem_expr/')
# tempdir <- paste0(exprdir, runid, '.temp/')
# libloc            <- '/ua/pliu/dev/RSEM/pRSEM/RLib/'
# fall_tr_crd       <- paste0(tempdir, runid, '_prsem.all_tr_crd')
##fall_tr_crd       <- paste0(tempdir, 'tmp.all_tr_crd')
# ftraining_tr_crd  <- paste0(tempdir, runid, '_prsem.training_tr_crd')
# fout              <- paste0('./', runid, '_prsem.all_tr_features')
# fisoforms_results <- paste0(exprdir, runid, '.isoforms.results')
# flanking_width    <- 500
# partition_model   <- 'lm4'
# fchipseq_peaks    <- paste0(tempdir, 'idr_target_vs_control.regionPeak.gz')
# fchipseq_target_signals <- paste0(tempdir, 'target.tagAlign.gz')
##fchipseq_target_signals <- paste0(tempdir, 'tmp.chrX.tagAlign.gz')
# fall_tr_gc        <- paste0(tempdir, runid, '_prsem.all_tr_gc')
# nthr              <- 20
# fraglen           <- 110

  .libPaths(c(libloc, .libPaths()))
  suppressMessages(library(data.table))
  suppressMessages(library(GenomicRanges))

  all_trdt <- fread(fall_tr_crd, header=T, sep="\t")

  rsemdt <- fread(fisoforms_results, header=T, sep="\t", select=c(
                  'transcript_id', 'effective_length', 'posterior_mean_count'))
  setnames(rsemdt, 1:3, c('trid', 'efflen',  'pme_count'))

  gcdt <- fread(fall_tr_gc, header=T, sep="\t")

  trdt <- merge(all_trdt, rsemdt, by='trid', all.x=T)
  trdt <- merge(trdt, gcdt, by='trid', all.x=T)

  trdt[, `:=`( tss = ifelse(strand == '+', start, end  ),
               tes = ifelse(strand == '+', end,   start))]
  tssdt  <- trdt[, list(chrom, tss, trid)]
  bodydt <- trdt[, list(chrom, start, end, trid)]
  tesdt  <- trdt[, list(chrom, tes, trid)]

  tssdt[,  `:=`( start = tss - flanking_width,
                 end   = tss + flanking_width ) ]
  bodydt[, `:=`( body_start = start + flanking_width + 1,
                 body_end   = end   - flanking_width - 1)]
  bodydt[, `:=`( start = ifelse(body_start <= body_end, body_start, body_end),
                 end   = ifelse(body_start <= body_end, body_end, body_start))]
  tesdt[,  `:=`( start = tes - flanking_width,
                 end   = tes + flanking_width ) ]

  pkdt <- data.table(read.table(gzfile(fchipseq_peaks), header=F, sep="\t",
                                colClasses=c('character', 'numeric', 'numeric',
                                             rep('NULL', 7))))
  setnames(pkdt, 1:3, c('chrom', 'start', 'end'))

  has_tss_pk_trids  <- getRegionPeakOLTrID(tssdt,  pkdt)
  has_body_pk_trids <- getRegionPeakOLTrID(bodydt, pkdt)
  has_tes_pk_trids  <- getRegionPeakOLTrID(tesdt,  pkdt)

  rddt <- data.table(read.table(gzfile(fchipseq_target_signals), header=F,
                                sep="\t", colClasses=c('character', 'numeric',
                                'numeric', rep('NULL', 2), 'character')))
  setnames(rddt, 1:4, c('chrom', 'start', 'end', 'strand'))
  tss_sigdt  <- countRegionSignal(tssdt,  rddt, fraglen, nthr, 'tss')
  body_sigdt <- countRegionSignal(bodydt, rddt, fraglen, nthr, 'body')
  tes_sigdt  <- countRegionSignal(tesdt,  rddt, fraglen, nthr, 'tes')

  trdt <- merge(trdt, tss_sigdt,  by='trid', all.x=T)
  trdt <- merge(trdt, body_sigdt, by='trid', all.x=T)
  trdt <- merge(trdt, tes_sigdt,  by='trid', all.x=T)

  trdt[, `:=`(tss_sig  = ifelse(is.na(tss_sig), 0.0, tss_sig),
              body_sig = ifelse(is.na(body_sig), 0.0, body_sig),
              tes_sig  = ifelse(is.na(tes_sig),  0.0, tes_sig ))]

  training_trdt <- fread(ftraining_tr_crd, header=T, sep="\t", select='trid')
  trdt[, `:=`( tss_pk  = ifelse(trid %in% has_tss_pk_trids,  1, 0),
               body_pk = ifelse(trid %in% has_body_pk_trids, 1, 0),
               tes_pk  = ifelse(trid %in% has_tes_pk_trids,  1, 0),
               is_training = ifelse(trid %in% training_trdt[, trid], 1, 0))]

  setkey(trdt, trid)
  trdt <- trdt[all_trdt[, trid]] ## keep the order of original trids
  write.table(trdt, fout, quote=F, sep="\t", col.names=T, row.names=F)
}

#
#  need to modify the way to calculate signal as # of nuc from fragment rather
#  than the number of read overlapping with selected region
#
#  tagAlign list the read
#  '+' strand, [start, start + read_length]
#  '-' strand, [end - read_length, end]
#
#  1. extend read to fragment
#  2. find fragmens overlapping target region
#  3. remove fragments that have middle position outside target region (as how
#     dpeak works)
#  4. count number of fragment nucleotide and average it by target region's
#     width to get signal
#
countRegionSignal <- function(regiondt, readdt, fraglen, nthr, prefix=''){
  regiondtl <- split(regiondt[, list(chrom, start, end, trid)],
                     regiondt[, chrom])
  readdtl   <- split(readdt, readdt[, chrom])

 #outdt <- rbindlist( lapply(names(regiondtl), countRegionSignalByChrom,
 #                           regiondtl, readdtl, fraglen, prefix))

  outdt <- rbindlist( mclapply(names(regiondtl), countRegionSignalByChrom,
                               regiondtl, readdtl, fraglen, prefix,
                               mc.cores = nthr))
  return(outdt)
}


countRegionSignalByChrom <- function(chrom, regiondtl, readdtl, fraglen,
                                     prefix) {
  regiondt <- copy(regiondtl[[chrom]])
  readdt   <- copy(readdtl[[chrom]])

  readdt[, frag_start := ifelse(strand == '+', start, end-fraglen)]
  readdt[, frag_end := frag_start + fraglen - 1]

  fragdt <- readdt[, list(chrom, frag_start, frag_end)]
  setnames(fragdt, 2:3, c('start', 'end'))
  fraggrs <- makeGRangesFromDataFrame(fragdt, ignore.strand=T)

  regiongrs <- makeGRangesFromDataFrame(regiondt[, list(chrom, start, end,
                                                        trid)],
                                        keep.extra.columns=T, ignore.strand=T)

  ol <- findOverlaps(regiongrs, fraggrs, type='any', ignore.strand=T)
  oldt <- data.table(query=queryHits(ol), subject=subjectHits(ol))

  oldt[, `:=`( region_start = regiondt[, start][query],
               region_end   = regiondt[, end][query],
               trid         = regiondt[, trid][query],
               frag_start   = fragdt[, start][subject],
               frag_end     = fragdt[, end][subject] )]

  oldt[, frag_mid := (frag_start + frag_end)/2]

  ## as dpeak, only select fragment which has mid position falling into region
  seloldt <- subset(oldt, (frag_mid >= region_start) & (frag_mid <= region_end))

  seloldt[, `:=`(start = ifelse(frag_start < region_start, region_start,
                                frag_start),
                 end   = ifelse(frag_end   > region_end,   region_end,
                                frag_end))]

  sigdt <- seloldt[, list(nuc = sum(end - start + 1)), by=trid]
  sigdt <- merge(sigdt, regiondt, by='trid', all.x=T)
  colname <- paste0(prefix, '_sig')
  sigdt[, eval(colname) := nuc/(end - start + 1)]
  sigdt[, `:=`(nuc=NULL, chrom=NULL, start=NULL, end=NULL)]
  return(sigdt)
}


genPriorByTSSPeak <- function(argv=NA){
  libloc           <- argv[1]
  fall_tr_features <- argv[2]
  fpval_LL         <- argv[3]
  fall_tr_prior    <- argv[4]

# libloc           <- '/ua/pliu/dev/RSEM/pRSEM/RLib/'
# fall_tr_features <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.stat/example_prsem.all_tr_features'
# fpval_LL         <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.stat/example_prsem.pval_LL'
# fall_tr_prior    <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.stat/example_prsem.all_tr_prior'

  .libPaths(c(libloc, .libPaths()))
  suppressMessages(library(data.table))

  alldt <- fread(fall_tr_features, header=T, sep="\t")
  alldt[, partition := tss_pk]

  trndt <- subset(alldt, is_training == 1)
  outdt <- getSampleAndPriorByTSSPeak(trndt, alldt)

  wpk_cnt  <- subset(trndt, tss_pk == 1)[, pme_count]
  nopk_cnt <- subset(trndt, tss_pk == 0)[, pme_count]
  wrs <- suppressWarnings(wilcox.test(wpk_cnt, nopk_cnt, alternative='greater',
                                      paired=F, exact=T))
  pval <- wrs$p.value
  loglikelihood <- unique(outdt[, loglikelihood])
  pval_LLdt <- data.table(pvalue=pval, loglikelihood=loglikelihood)

  write.table(alldt, fall_tr_features, quote=F, sep="\t", col.names=T,
              row.names=F)
  write.table(pval_LLdt, fpval_LL, quote=F, sep="\t", col.names=T, row.names=F)
  write.table(outdt[, list(prior, trid)], fall_tr_prior, quote=F, sep="  # ",
              col.names=F, row.names=F)
}


prepTSSPeakFeatures <- function(argv=NA) {
  libloc            <- argv[1]
  fall_tr_crd       <- argv[2]
  ftraining_tr_crd  <- argv[3]
  fout              <- argv[4]
  fisoforms_results <- argv[5]
  flanking_width    <- as.numeric(argv[6])
  fchipseq_peaks    <- argv[7]

# libloc            <- '/ua/pliu/dev/RSEM/pRSEM/RLib/'
# fall_tr_crd       <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.temp/example_prsem.all_tr_crd'
# ftraining_tr_crd  <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.temp/example_prsem.training_tr_crd'
# fout              <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.temp/example_prsem.all_tr_features'
# fisoforms_results <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.isoforms.results'
# flanking_width <- 500
# fchipseq_peaks <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/example.temp/idr_target_vs_control.regionPeak.gz'

  .libPaths(c(libloc, .libPaths()))
  suppressMessages(library(data.table))
  suppressMessages(library(GenomicRanges))

  rsemdt <- fread(fisoforms_results, header=T, sep="\t",
                  select=c('transcript_id', 'posterior_mean_count'))
  setnames(rsemdt, 1:2, c('trid', 'pme_count'))

  intrdt <- fread(fall_tr_crd, header=T, sep="\t")
  trdt <- merge(intrdt, rsemdt, by='trid', all.x=T)
  trdt[, tss := ifelse(strand=='+', start, end)]
  tssdt <- trdt[, list(chrom, tss, trid)]
  tssdt[, `:=`( start = tss - flanking_width,
                end   = tss + flanking_width)]

  pkdt <- tryCatch({
    data.table(read.table(gzfile(fchipseq_peaks), header=F, sep="\t",
                          colClasses=c('character', 'numeric',
                                       'numeric', rep('NULL', 7))))
  }, error = function(err) {
    message(paste0("\nFail to read file: ", fchipseq_peaks, "\n"))
    message(err)
    return(NA)
  })

  setnames(pkdt, 1:3, c('chrom', 'start', 'end'))

  has_pk_trids <- getRegionPeakOLTrID(tssdt, pkdt)
  training_trids <- fread(ftraining_tr_crd, header=T, sep="\t", select='trid'
                         )[, trid]

  trdt[, `:=`( tss_pk      = ifelse(trid %in% has_pk_trids,   1, 0),
               is_training = ifelse(trid %in% training_trids, 1, 0) )]

  setkey(trdt, trid)
  trdt <- trdt[intrdt[, trid]]  ## keep the order of original trid
  write.table(trdt, fout, quote=F, sep="\t", col.names=T, row.names=F)
}


getRegionPeakOLTrID <- function(regiondt, peakdt) {
  regiongrs <- makeGRangesFromDataFrame(regiondt[, list(chrom, start, end,
                                                        trid)],
                                        keep.extra.columns=T, ignore.strand=T)

  peakgrs <- makeGRangesFromDataFrame(peakdt[, list(chrom, start, end)],
                                    ignore.strand=T)
  olgrs <- subsetByOverlaps(regiongrs, peakgrs, type='any', ignore.strand=T)
  has_peak_trids <- unique(mcols(olgrs)$trid)
  return(has_peak_trids)
}


selTrainingTr <- function(argv=NA) {
  libloc   <- argv[1]
  fin_tr   <- argv[2]
  fin_exon <- argv[3]
  min_mpp  <- as.numeric(argv[4])
  flanking_width <- as.numeric(argv[5])
  fout     <- argv[6]

# libloc   <- '/ua/pliu/dev/RSEM/pRSEM/RLib/'
# fin_tr   <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/test.temp/test_prsem.all_tr_crd'
# fin_exon <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/test.temp/test_prsem.all_exon_crd'
# min_mpp  <- 0.8
# flanking_width <- 500
# fout     <- '/tier2/deweylab/scratch/pliu/dev/pRSEM/rsem_expr/test.temp/test_prsem.training_tr_crd'

  .libPaths(c(libloc, .libPaths()))
  suppressMessages(library(data.table))
  suppressMessages(library(GenomicRanges))

  alltrdt <- fread(fin_tr, header=T, sep="\t")
  alltrdt[, tss := ifelse(strand == '+', start, end)]
  highmppdt <- subset(alltrdt, (! is.na(tss_mpp )) & ( tss_mpp  > min_mpp ) &
                               (! is.na(body_mpp)) & ( body_mpp > min_mpp ) &
                               (! is.na(tes_mpp )) & ( tes_mpp  > min_mpp ) )

  ## select tr that are not nested with other tr regardless of strand
  nested_trids <- getTrTrOLTrID(highmppdt, alltrdt, oltype='within',
                                ignore_strand=T)
  not_nested_trdt <- subset(highmppdt, ! trid %in% nested_trids)

  ## select tr that not have exon all nested with the union of other tr's exons
  ## regardless of strand
  allexondt <- fread(fin_exon, header=T, sep="\t")
  exon_all_ol_trids <- getExonsAllOLTrID(not_nested_trdt[, trid], allexondt)
  seltrdt <- subset(not_nested_trdt, ! trid %in% exon_all_ol_trids)

  seltr_tss_region_dt <- copy(seltrdt)
  seltr_tss_region_dt[, `:=`( start = tss - flanking_width,
                              end   = tss + flanking_width )]

  alltr_tss_dt <- copy(alltrdt)
  alltr_tss_dt[, `:=`(start = tss, end=tss)]

  tss_region_ol_trids <- getTrTrOLTrID(seltr_tss_region_dt, alltr_tss_dt,
                                       'any', ignore_strand=T)

  outdt <- subset(seltrdt, ! trid %in% tss_region_ol_trids)
  write.table(outdt, fout, quote=F, sep="\t", col.names=T, row.names=F)
}


getExonsAllOLTrID <- function(query_trids, allexondt) {
  allexongrs <- makeGRangesFromDataFrame(allexondt, keep.extra.columns=T,
                                         ignore.strand=F)
  queryexongrs <- subset(allexongrs, mcols(allexongrs)$trid %in% query_trids)
  ol <- findOverlaps(queryexongrs, allexongrs, type='within', ignore.strand=T)

  oldt <- data.table(query = queryHits(ol), subject=subjectHits(ol))
  oldt[, `:=`( query_trid = mcols(queryexongrs)$trid[query],
               query_exon_index = mcols(queryexongrs)$exon_index[query],
               subject_trid = mcols(allexongrs)$trid[subject]
             )]
  oldt <- subset(oldt, query_trid != subject_trid)
  nolexondt <- oldt[, list(nolexon = length(unique(query_exon_index))),
                      by=query_trid]
  subexondt  <- subset(allexondt, trid %in% nolexondt[, query_trid])
  subnexondt <- subexondt[, list(nexon = length(unique(exon_index))), by=trid]
  setnames(nolexondt, 1, 'trid')
  nolexondt <- merge(nolexondt, subnexondt, by='trid', all.x=T)
  exon_all_ol_trids <- subset(nolexondt, nolexon == nexon)[, trid]

  return(exon_all_ol_trids)
}


getTrTrOLTrID <- function(querydt, subjectdt, oltype, ignore_strand) {
  querygrs <- makeGRangesFromDataFrame(
                  querydt[, list(chrom, strand, start, end, trid)],
                  keep.extra.columns=T, ignore.strand=F)
  subjectgrs <- makeGRangesFromDataFrame(
                  subjectdt[, list(chrom, strand, start, end, trid)],
                  keep.extra.columns=T, ignore.strand=F)

  ol <- findOverlaps(querygrs, subjectgrs, type=oltype,
                     ignore.strand=ignore_strand)

  oldt <- data.table(query=queryHits(ol), subject=subjectHits(ol))
  oldt[, `:=`(query_trid   = mcols(querygrs)$trid[query],
              subject_trid = mcols(subjectgrs)$trid[subject] )]
  ol_trids <- subset(oldt, query_trid != subject_trid)[, query_trid]
  return(unique(ol_trids))
}


rdirichlet_multinomial <- function(alpha, n) {
  theta <- rdirichlet(alpha)
  return(as.vector(rmultinom(1, n, theta)))
}


rdirichlet <- function(alpha) {
  x <- rgamma(length(alpha), alpha)
  return(x / sum(x))
}


getFitByMLDM <- function(counts, partition) {
  initial_alpha <- rep(1, nlevels(partition))
  fit <- optim(initial_alpha,
               partitioned_log_likelihood,
               partitioned_log_likelihood_gradient,
               counts=counts,
               partition=partition,
               method="L-BFGS-B",
               lower=0.0001,
               upper=1e+4,
               control=list(fnscale=-1))

  return(fit)
}


getSampleByDM <- function(par, counts, partition) {
  ml_alpha <- par
  names(ml_alpha) <- levels(partition)

  alpha <- ml_alpha[partition]
  sample <- rdirichlet_multinomial(alpha, sum(counts))
  return(sample)
}


partitioned_log_likelihood <- function(alpha, counts, partition) {

  component_counts <- table(partition)
  N <- sum(counts)

  return(lgamma(N + 1) - sum(lgamma(counts + 1)) +
         lgamma(component_counts %*% alpha) -
         lgamma(N + component_counts %*% alpha) +
         sum(lgamma(counts + alpha[partition])) -
         component_counts %*% lgamma(alpha))
}


partitioned_log_likelihood_gradient <- function(alpha, counts, partition) {

  component_counts <- table(partition)

  N <- sum(counts)
  alpha_sum <- sum(component_counts * alpha)

  return(component_counts * (digamma(alpha_sum) -
                             digamma(N + alpha_sum) -
                             digamma(alpha)) +
         tapply(digamma(counts + alpha[partition]), partition, sum))
}


#
# create a new partition with give parition and data
# the breaks from old partition will be kept
# lower and upper bounds will be from give data
#
createPartitionForNewData <- function(partition, data){
  labs <- levels(partition)
  lower <- as.numeric( sub("\\((.+),.*", "\\1", labs) )
  data_range <- range(data)
  eps <- 1.0e-4
  lower[1] <- data_range[1] - eps
  breaks <- c(lower, data_range[2] + eps)
  new_partition <- cut(data, breaks)

  return(new_partition)
}


#getSampleAndPriorByOnePartition <- function(trndt, tstdt) {
# trn_partition <- factor(rep(0, nrow(trndt)))
# fit <- getFitByMLDM(trndt[, pme_count], trn_partition)

##cat('priors:', format(fit$par, digits=2, nsmall=3), "\n")
#
# tst_partition <- factor(rep(0, nrow(tstdt)))
# prior <- fit$par[tst_partition]

# outdt <- tstdt[, list(trid, pme_count)]
# outdt[, `:=`(partition = tst_partition,
#              sample    = getSampleByDM(fit$par, pme_count, tst_partition),
#              prior     = prior)]

##cat("training set's partition:")
##print(table(trn_partition))
##cat("testing set's prior:")
##print(table(prior))

# return(outdt)
#}


getSampleAndPriorByTSSPeak <- function(trndt, tstdt) {
  trn_partition <- factor(trndt[, tss_pk])
  fit <- getFitByMLDM(trndt[, pme_count], trn_partition)

  tst_partition <- factor(tstdt[, tss_pk])
  prior <- fit$par[tst_partition]

  outdt <- tstdt[, list(trid, pme_count)]
  outdt[, `:=`(partition = tst_partition,
               sample    = getSampleByDM(fit$par, pme_count, tst_partition),
               prior     = prior,
               loglikelihood = fit$value )]

  cat("training set's partition:")
  print(table(trn_partition))
  cat("testing set's prior:")
  print(table(prior))

  return(outdt)
}


getSampleAndPriorByLM <- function(trndt, tstdt, nbin=NULL) {
	if (is.null(nbin)) nbin=3
  frm <- formula( paste0('log10_count~',
         'tss_pk  + tss_pk:log10_tss_sig   + no_tss_pk:log10_tss_sig   +',
         'body_pk + body_pk:log10_body_sig + no_body_pk:log10_body_sig +',
         'tes_pk  + tes_pk:log10_tes_sig   + no_tes_pk:log10_tes_sig   +',
         'log10_eff_len + log10_GC_ov_mean') )

  trn_lm <- lm(frm, trndt)
  trn_prd <- as.vector(predict(trn_lm, trndt))
  trn_partition <- cut(trn_prd, nbin)

	fit <- getFitByMLDM(trndt[, pme_count], trn_partition)

	tst_prd <- as.vector(predict(trn_lm, tstdt))
	tst_partition <- createPartitionForNewData(trn_partition, tst_prd)

	prior <- fit$par[tst_partition]

  outdt <- tstdt[, list(trid, pme_count)]
  outdt[, `:=`(partition = tst_partition,
               sample    = getSampleByDM(fit$par, pme_count, tst_partition),
               prior     = prior)]

  cat("training set's partition:")
  print(table(trn_partition))
  cat("testing set's prior:")
  print(table(prior))

  return(outdt)
}


getSampleAndPriorByPeakLM <- function(trndt, tstdt, nbin=NULL, lm_on_wpk=NULL) {
	if (is.null(nbin)) nbin=2
  if (is.null(lm_on_wpk)) lm_on_wpk=T
  lm_pk_type <- ifelse(lm_on_wpk, 1, 0)

	slim_trndt <- trndt[, list(trid, pme_count, tss_pk)]
	slim_tstdt <- tstdt[, list(trid, pme_count, tss_pk)]

  frm <- formula( paste0('log10_count~',
         'log10_tss_sig +',
         'body_pk + body_pk:log10_body_sig + no_body_pk:log10_body_sig +',
         'tes_pk  + tes_pk:log10_tes_sig   + no_tes_pk:log10_tes_sig   +',
         'log10_eff_len + log10_GC_ov_mean') )

	pk_trn_dt <- subset(trndt, tss_pk == lm_pk_type)
	pk_trn_lm <- lm(frm, pk_trn_dt)
	pk_trn_prd <- as.vector(predict(pk_trn_lm, pk_trn_dt))
	pk_trn_partition <- cut(pk_trn_prd, nbin)
	pk_trn_dt[, partition:=pk_trn_partition]

  slim_trndt <- merge(slim_trndt, pk_trn_dt[, list(trid, partition)], by='trid',
                      all=T)
  slim_trndt <- slim_trndt[trndt[,trid]]

	slim_trndt[, partition:=factor(ifelse(tss_pk==(1-lm_pk_type), 0, partition))]

	fit <- getFitByMLDM(slim_trndt[, pme_count], slim_trndt[, partition])

	pk_tst_dt <- subset(tstdt, tss_pk == lm_pk_type)
	pk_tst_prd <- as.vector(predict(pk_trn_lm, pk_tst_dt))
	pk_tst_partition <- createPartitionForNewData(pk_trn_partition,
	                                                   pk_tst_prd)
	pk_tst_dt[, partition:=pk_tst_partition]

  slim_tstdt <- merge(slim_tstdt, pk_tst_dt[, list(trid, partition)], by='trid',
                      all=T)
  slim_tstdt <- slim_tstdt[tstdt[,trid]]
	slim_tstdt[, partition:=factor(ifelse(tss_pk==(1-lm_pk_type), 0, partition))]
	prior <- fit$par[slim_tstdt[, partition]]

  outdt <- slim_tstdt[, list(trid, pme_count, partition)]
  outdt[, `:=`(sample = getSampleByDM(fit$par, pme_count, partition),
               prior  = prior)]

  cat("training set's partition:")
  print(table(slim_trndt[,partition]))
  cat("testing set:")
  print(table(prior))
  cat("\n")

	return(outdt)
}


getSampleAndPriorByPeakLogitNoPeak <- function(trndt, tstdt) {
	slim_trndt <- trndt[, list(trid, pme_count, tss_pk)]
	slim_tstdt <- tstdt[, list(trid, pme_count, tss_pk)]

  setkey(slim_trndt, trid)
  setkey(slim_tstdt, trid)

	nopk_trn_dt <- subset(trndt, tss_pk == 0)
	nopk_tst_dt <- subset(tstdt, tss_pk == 0)
	nopk_trn_dt[, has_count:=ifelse(pme_count > 0, 1, 0)]
	nopk_tst_dt[, has_count:=ifelse(pme_count > 0, 1, 0)]

  frm <- formula( paste0('has_count~',
         'log10_tss_sig +',
         'body_pk + body_pk:log10_body_sig + no_body_pk:log10_body_sig +',
         'tes_pk  + tes_pk:log10_tes_sig   + no_tes_pk:log10_tes_sig   +',
         'log10_eff_len + log10_GC_ov_mean') )
	prt_levels <- c('no pk, no cnt', 'no pk, has cnt', 'w/ pk')

	nopk_trn_glm <- glm(frm, data=nopk_trn_dt, family='binomial')
	trn_prob <- predict(nopk_trn_glm, nopk_trn_dt, type='response')
	nopk_trn_dt[, `:=`(logit_prob=trn_prob,
	                   partition=ifelse(trn_prob > 0.5, 'no pk, has cnt',
										                  'no pk, no cnt') )]

  slim_trndt <- merge(slim_trndt,
                      nopk_trn_dt[, list(trid, partition, logit_prob)],
                      by='trid', all=T)
  slim_trndt <- slim_trndt[trndt[, trid]]

	slim_trndt[, partition:=factor(ifelse(tss_pk==1, 'w/ pk', partition),
	                               levels=prt_levels)]


	fit <- getFitByMLDM(slim_trndt[, pme_count], slim_trndt[, partition])

	tst_prob <- predict(nopk_trn_glm, nopk_tst_dt, type='response')
	nopk_tst_dt[, `:=`(logit_prob=tst_prob,
	                   partition=ifelse(tst_prob > 0.5, 'no pk, has cnt',
										                  'no pk, no cnt'))]

  slim_tstdt <- merge(slim_tstdt,
                      nopk_tst_dt[, list(trid, partition, logit_prob)],
                      by='trid', all=T)
  slim_tstdt <- slim_tstdt[tstdt[, trid]]
	slim_tstdt[, partition:=factor(ifelse(tss_pk==1, 'w/ pk', partition),
	                               levels=prt_levels)]

	prior <- fit$par[slim_tstdt[, partition]]

  outdt <- slim_tstdt[, list(trid, pme_count, partition)]
  outdt[, `:=`(sample = getSampleByDM(fit$par, pme_count, partition),
               prior  = prior)]

  cat("training set's partition:")
  print(table(slim_trndt[,partition]))
  cat("testing set:")
  print(table(prior))
  cat("\n")

	return(outdt)
}


getSampleAndPriorByLM3 <- function(trndt, tstdt) {
	return(getSampleAndPriorByLM(trndt, tstdt, nbin=3))
}

getSampleAndPriorByLM4 <- function(trndt, tstdt) {
	return(getSampleAndPriorByLM(trndt, tstdt, nbin=4))
}

getSampleAndPriorByLM5 <- function(trndt, tstdt) {
	return(getSampleAndPriorByLM(trndt, tstdt, nbin=5))
}

getSampleAndPriorByLM6 <- function(trndt, tstdt) {
	return(getSampleAndPriorByLM(trndt, tstdt, nbin=6))
}


getSampleAndPriorByNoPeakLM2Peak <- function(trndt, tstdt) {
	return(getSampleAndPriorByPeakLM(trndt, tstdt, nbin=2, lm_on_wpk=T))
}

getSampleAndPriorByNoPeakLM3Peak <- function(trndt, tstdt) {
	return(getSampleAndPriorByPeakLM(trndt, tstdt, nbin=3, lm_on_wpk=T))
}

getSampleAndPriorByNoPeakLM4Peak <- function(trndt, tstdt) {
	return(getSampleAndPriorByPeakLM(trndt, tstdt, nbin=4, lm_on_wpk=T))
}

getSampleAndPriorByNoPeakLM5Peak <- function(trndt, tstdt) {
	return(getSampleAndPriorByPeakLM(trndt, tstdt, nbin=5, lm_on_wpk=T))
}


getSampleAndPriorByPeakLM2NoPeak <- function(trndt, tstdt) {
	return(getSampleAndPriorByPeakLM(trndt, tstdt, nbin=2, lm_on_wpk=F))
}

getSampleAndPriorByPeakLM3NoPeak <- function(trndt, tstdt) {
	return(getSampleAndPriorByPeakLM(trndt, tstdt, nbin=3, lm_on_wpk=F))
}

getSampleAndPriorByPeakLM4NoPeak <- function(trndt, tstdt) {
	return(getSampleAndPriorByPeakLM(trndt, tstdt, nbin=4, lm_on_wpk=F))
}

getSampleAndPriorByPeakLM5NoPeak <- function(trndt, tstdt) {
	return(getSampleAndPriorByPeakLM(trndt, tstdt, nbin=5, lm_on_wpk=F))
}


main()
#selTrainingTr()
#prepTSSPeakFeatures()
#genPriorByTSSPeak()
#prepPeakSignalGCLenFeatures()
#genPriorByPeakSignalGCLen()
#prepTSSSignalFeatures()
#genPriorByCombinedTSSSignals()
#prepMultiTargetsFeatures()
