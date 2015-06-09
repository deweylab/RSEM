#
#  pliu 20150608
# 
#  module for processing RNA-seq data
#

main <- function() {
  name2func <- list(
    'selTrainingTr' = selTrainingTr
  )

  argv <- commandArgs(trailingOnly=T)
  name2func[[argv[1]]](argv[2:length(argv)])
}


selTrainingTr <- function(argv) {
  libloc  <- argv[1]
  fin     <- argv[2]
  min_mpp <- argv[3]
  flanking_width <- argv[4]
  fout    <- argv[5]

 #libloc  <- '/ua/pliu/dev/RSEM/pRSEM/Rlib/'
 #fin     <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.alltrcrd'
 #min_mpp <- 0.8
 #flanking_width <- 500
 #fout    <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.training_tr'

  checkInstallCRAN('data.table', libloc)
  checkInstallBioc('GenomicRanges',  libloc)
  suppressPackageStartupMessages(library(data.table, 
                                         lib.loc=c(.libPaths(), libloc)))
  suppressPackageStartupMessages(library(GenomicRanges, 
                                         lib.loc=c(.libPaths(), libloc)))

  alldt <- fread(fin, header=T, sep="\t")
  alldt[, tss := ifelse(strand == '+', start, end)]
  highmppdt <- subset(alldt, (! is.na(tss_mpp )) & ( tss_mpp  >= min_mpp ) & 
                             (! is.na(body_mpp)) & ( body_mpp >= min_mpp ) & 
                             (! is.na(tes_mpp )) & ( tes_mpp  >= min_mpp ) )

  ## select tr that don't overlap with other tr
  ol_trid <- getOLTrID(highmppdt, alldt)

  ## select tr that don't have other tr's TSS within its [TSS-width, TSS+width]
  seltrdt <- subset(highmppdt, ! trid %in%  ol_trid)

  seltr_tss_region_dt <- copy(seltrdt)
  seltr_tss_region_dt[, `:=`( start = tss - flanking_width,
                              end   = tss + flanking_width )]

  alltr_tss_dt <- copy(alldt)  
  alltr_tss_dt[, `:=`(start = tss, end=tss)]

  tss_region_ol_trid <- getOLTrID(seltr_tss_region_dt, alltr_tss_dt)

  outdt <- subset(seltrdt, ! trid %in% tss_region_ol_trid)
  write.table(outdt, fout, quote=F, sep="\t", col.names=T, row.names=F)
}


getOLTrID <- function(querydt, subjectdt) {
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


main()
#selTrainingTr()
