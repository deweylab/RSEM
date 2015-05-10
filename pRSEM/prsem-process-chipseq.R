#!/usr/bin/env Rscript

#
#  pliu 20150509
# 
#  module for processing ChIP-seq data
#


main <- function() {
  name2func <- list(
    'guessFqEncoding' = guessFqEncoding
  )

  argv <- commandArgs(trailingOnly=T)
  name2func[[argv[1]]](argv[2:length(argv)])
  
 #guessFqEncoding('dum')
}


guessFqEncoding <- function(argv){

  if ( ! 'data.table' %in% rownames(installed.packages() ) ) {
    install.packages('data.table')
  }
  if ( ! 'ShortRead' %in% rownames(installed.packages()) ) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("ShortRead")
  }

  library(data.table)
  library(ShortRead)

  nthr      <- strtoi(argv[1])
  s_infiles <- argv[2]
  fout      <- argv[3]

 #nthr <- 16
 #s_infiles <- '/tier2/deweylab/scratch/pliu/dev/MelPol2IggmusRep1.fastq.gz,/tier2/deweylab/scratch/pliu/dev/MelPol2IggmusRep2.fastq.gz,/tier2/deweylab/scratch/pliu/dev/MelInputIggmusRep1.fastq.gz'
 #fout <- '/tier2/deweylab/scratch/pliu/dev/rsem_expr/test.temp/test_prsem.fq_encod'

  files <- strsplit(s_infiles, ',', fixed=T)[[1]]

 #lapply(files, guessFqEncodingByFile)

  ## cannot use mclapply here, ShortRead::qa function bound with BiocParallel
  ## have to use bplapply and define single core when call qa
  ## this new feature only appear in the new version of ShortRead
  register(MulticoreParam(workers=nthr))
  outdt <- rbindlist(bplapply(files, guessFqEncodingByFile))

  write.table(outdt, fout, sep="\t", quote=F, col.names=T, row.names=F)
  cat('File written:', fout, "\n")
}


## exam and guess Fastq's quality score's format
guessFqEncodingByFile <- function(fq) {
  qual <- qa(fq, BPPARAM=registered()$SerialParam)
  bq <- qual[['baseQuality']]
  score <- subset(bq, bq$count > 0, select=c('score'))$score

  encod <- 'unknown'
  if ( any( strsplit(intToUtf8(33:58), '')[[1]] %in% score ) ) {
    encod <- '--phred33-quals'
  } else if ( any( strsplit(intToUtf8(59:64), '')[[1]] %in% score ) ) {
   #encod <- 'solexa ver. <1.3'
    encod <- '--solexa-quals'
  } else if ( any( strsplit(intToUtf8(75:104), '')[[1]] %in% score ) ) {
    encod <- '--phred64-quals'
  }

  return(list(file=fq, encoding=encod))
}


system.time( main() )
