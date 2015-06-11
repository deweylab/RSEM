#
#  pliu 20150608
# 
#  module for processing RNA-seq data
#

main <- function() {
  name2func <- list(
    'selTrainingTr'               = selTrainingTr,
    'prepTSSPeakFeatures'         = prepTSSPeakFeatures,
    'prepPeakSignalGCLenFeatures' = prepPeakSignalGCLenFeatures,
    'genPriorByTSSPeak'           = genPriorByTSSPeak,
    'genPriorByPeakSignalGCLen'   = genPriorByPeakSignalGCLen
  )

  argv <- commandArgs(trailingOnly=T)
  name2func[[argv[1]]](argv[2:length(argv)])
}


genPriorByPeakSignalGCLen <- function(argv=NA) {
  libloc           <- argv[1]
  fall_tr_features <- argv[2]
  partition_model  <- argv[3]
  fout             <- argv[4]

# libloc <- '/ua/pliu/dev/RSEM/pRSEM/Rlib/'
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

  checkInstallCRAN('data.table', libloc)
  suppressMessages(library(data.table, lib.loc=c(.libPaths(), libloc)))

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

# libloc            <- '/ua/pliu/dev/RSEM/pRSEM/Rlib/'
# fall_tr_crd       <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.all_tr_crd'
# ftraining_tr_crd  <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.training_tr_crd' 
# fout              <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.all_tr_features' 
# fisoforms_results <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.isoforms.results' 
# flanking_width    <- 500 
# partition_model   <- 'pk_lgtnopk'
# fchipseq_peaks    <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/idr_targetRep0_vs_controlRep0.regionPeak.gz'
# fchipseq_target_signals <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/targetRep0.tagAlign.gz' 
# fall_tr_gc        <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.all_tr_gc'

  checkInstallCRAN('data.table', libloc)
  checkInstallBioc('GenomicRanges', libloc)
  suppressMessages(library(data.table,    lib.loc=c(.libPaths(), libloc)))
  suppressMessages(library(GenomicRanges, lib.loc=c(.libPaths(), libloc)))

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
                                colClasses=c('character', 'integer', 'integer',
                                             rep('NULL', 7))))
  setnames(pkdt, 1:3, c('chrom', 'start', 'end'))

  has_tss_pk_trids  <- getRegionPeakOLTrID(tssdt,  pkdt)
  has_body_pk_trids <- getRegionPeakOLTrID(bodydt, pkdt)
  has_tes_pk_trids  <- getRegionPeakOLTrID(tesdt,  pkdt)

  sigdt <- data.table(read.table(gzfile(fchipseq_target_signals), header=F, 
                                 sep="\t", colClasses=c('character', 'integer',
                                 'integer', rep('NULL', 3))))
  setnames(sigdt, 1:3, c('chrom', 'start', 'end'))
  tss_sigdt  <- countRegionSignal(tssdt,  sigdt, 'tss')
  body_sigdt <- countRegionSignal(bodydt, sigdt, 'body')
  tes_sigdt  <- countRegionSignal(tesdt,  sigdt, 'tes')

  trdt <- merge(trdt, tss_sigdt,  by='trid', all.x=T)
  trdt <- merge(trdt, body_sigdt, by='trid', all.x=T)
  trdt <- merge(trdt, tes_sigdt,  by='trid', all.x=T)

  training_trdt <- fread(ftraining_tr_crd, header=T, sep="\t", select='trid')
  trdt[, `:=`( tss_pk  = ifelse(trid %in% has_tss_pk_trids,  1, 0),
               body_pk = ifelse(trid %in% has_body_pk_trids, 1, 0),
               tes_pk  = ifelse(trid %in% has_tes_pk_trids,  1, 0),
               is_training = ifelse(trid %in% training_trdt[, trid], 1, 0))]

  setkey(trdt, trid)
  trdt <- trdt[all_trdt[, trid]] ## keep the order of original trids
  write.table(trdt, fout, quote=F, sep="\t", col.names=T, row.names=F)
}


countRegionSignal <- function(regiondt, sigdt, prefix=''){
  regiongrs <- makeGRangesFromDataFrame(regiondt[, list(chrom, start, end, 
                                                        trid)], 
                                        keep.extra.columns=T, ignore.strand=T)

  peakgrs <- makeGRangesFromDataFrame(sigdt[, list(chrom, start, end)], 
                                      ignore.strand=T)
  nol <- countOverlaps(GNCList(regiongrs), GNCList(peakgrs), type='any', 
                       ignore.strand=T)
  outdt <- data.table(trid = mcols(regiongrs)$trid)
  outdt[, eval(paste0(prefix, '_sig')) := nol/width(regiongrs)]
  return(outdt)
}


genPriorByTSSPeak <- function(argv=NA){
  libloc           <- argv[1]
  fall_tr_features <- argv[2]
  fall_tr_prior    <- argv[3]

 #libloc           <- '/ua/pliu/dev/RSEM/pRSEM/Rlib/' 
 #fall_tr_features <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.all_tr_features'
 #fall_tr_prior    <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.all_tr_prior'

  selcols <- c('trid', 'tss_pk', 'pme_count', 'is_training')
  all_trdt <- fread(fall_tr_features, header=T, sep="\t", select=selcols)
  training_trdt <- subset(all_trdt, is_training == 1)
  outdt <- getSampleAndPriorByTSSPeak(training_trdt, all_trdt)
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

# libloc            <- '/ua/pliu/dev/RSEM/pRSEM/Rlib/'
# fall_tr_crd       <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.all_tr_crd'
# ftraining_tr_crd  <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.training_tr_crd'
# fout              <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.all_tr_features' 
# fisoforms_results <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.isoforms.results' 
# flanking_width <- 500 
# fchipseq_peaks <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/idr_targetRep0_vs_controlRep0.regionPeak.gz'

  checkInstallCRAN('data.table', libloc)
  checkInstallBioc('GenomicRanges', libloc)
  suppressMessages(library(data.table,    lib.loc=c(.libPaths(), libloc)))
  suppressMessages(library(GenomicRanges, lib.loc=c(.libPaths(), libloc)))

  rsemdt <- fread(fisoforms_results, header=T, sep="\t", 
                  select=c('transcript_id', 'posterior_mean_count'))
  setnames(rsemdt, 1:2, c('trid', 'pme_count'))

  intrdt <- fread(fall_tr_crd, header=T, sep="\t")
  trdt <- merge(intrdt, rsemdt, by='trid', all.x=T)
  trdt[, tss := ifelse(strand=='+', start, end)]
  tssdt <- trdt[, list(chrom, tss, trid)]
  tssdt[, `:=`( start = tss - flanking_width,
                end   = tss + flanking_width)]

  pkdt <- data.table(read.table(gzfile(fchipseq_peaks), header=F, sep="\t", 
                                colClasses=c('character', 'integer', 'integer',
                                             rep('NULL', 7))))
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
  olgrs <- subsetByOverlaps(GNCList(regiongrs), GNCList(peakgrs), type='any',
                            ignore.strand=T)
  has_peak_trids <- unique(mcols(olgrs)$trid)
  return(has_peak_trids)
}


selTrainingTr <- function(argv) {
  libloc  <- argv[1]
  fin     <- argv[2]
  min_mpp <- as.numeric(argv[3])
  flanking_width <- as.numeric(argv[4])
  fout    <- argv[5]

 #libloc  <- '/ua/pliu/dev/RSEM/pRSEM/Rlib/'
 #fin     <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.alltrcrd'
 #min_mpp <- 0.8
 #flanking_width <- 500
 #fout    <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.training_tr'

  checkInstallCRAN('data.table', libloc)
  checkInstallBioc('GenomicRanges',  libloc)
  suppressMessages(library(data.table,    lib.loc=c(.libPaths(), libloc)))
  suppressMessages(library(GenomicRanges, lib.loc=c(.libPaths(), libloc)))

  alldt <- fread(fin, header=T, sep="\t")
  alldt[, tss := ifelse(strand == '+', start, end)]
  highmppdt <- subset(alldt, (! is.na(tss_mpp )) & ( tss_mpp  >= min_mpp ) & 
                             (! is.na(body_mpp)) & ( body_mpp >= min_mpp ) & 
                             (! is.na(tes_mpp )) & ( tes_mpp  >= min_mpp ) )

  ## select tr that don't overlap with other tr
  ol_trid <- getTrTrOLTrID(highmppdt, alldt)

  ## select tr that don't have other tr's TSS within its [TSS-width, TSS+width]
  seltrdt <- subset(highmppdt, ! trid %in%  ol_trid)

  seltr_tss_region_dt <- copy(seltrdt)
  seltr_tss_region_dt[, `:=`( start = tss - flanking_width,
                              end   = tss + flanking_width )]

  alltr_tss_dt <- copy(alldt)  
  alltr_tss_dt[, `:=`(start = tss, end=tss)]

  tss_region_ol_trid <- getTrTrOLTrID(seltr_tss_region_dt, alltr_tss_dt)

  outdt <- subset(seltrdt, ! trid %in% tss_region_ol_trid)
  write.table(outdt, fout, quote=F, sep="\t", col.names=T, row.names=F)
}


getTrTrOLTrID <- function(querydt, subjectdt) {
  querygrs <- makeGRangesFromDataFrame(
                  querydt[, list(chrom, strand, start, end, trid)], 
                  keep.extra.columns=T, ignore.strand=F)
  subjectgrs <- makeGRangesFromDataFrame(
                  subjectdt[, list(chrom, strand, start, end, trid)], 
                  keep.extra.columns=T, ignore.strand=F)

  ## select tr not overlap with any other tr
  ol <- findOverlaps(GNCList(querygrs), GNCList(subjectgrs), type='any', 
                     ignore.strand=T)

  oldt <- data.table(query=queryHits(ol), subject=subjectHits(ol))
  oldt[, `:=`(query_trid   = mcols(querygrs)$trid[query],
              subject_trid = mcols(subjectgrs)$trid[subject] )]
  ol_trid <- subset(oldt, query_trid != subject_trid)[, query_trid] 
  return(unique(ol_trid))
}


checkInstallCRAN <- function(pkg, lib) {
  if ( ! pkg %in% rownames(installed.packages(lib.loc=c(.libPaths(), lib)))){
    cat("\ninstall R package", pkg, 'to', lib, "\n\n")
    install.packages(pkgs=pkg, lib=lib, quiet=T)
  }
}

checkInstallBioc <- function(pkg, lib) {
  if ( ! pkg %in% rownames(installed.packages(lib.loc=c(.libPaths(), lib)))){
    cat("\ninstall R package", pkg, 'to', lib, "\n\n")
    source("http://bioconductor.org/biocLite.R")
    biocLite(pkgs=pkg, lib=lib, quiet=T)
  }
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

 #cat('priors:', format(fit$par, digits=2, nsmall=3), "\n")
  
  tst_partition <- factor(tstdt[, tss_pk])
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
#prepTSSPeakFeatures()
#genPriorByTSSPeak()
#prepPeakSignalGCLenFeatures()
#genPriorByPeakSignalGCLen()