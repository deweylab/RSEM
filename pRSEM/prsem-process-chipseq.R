#!/usr/bin/env Rscript

#
#  pliu 20150509
# 
#  module for processing ChIP-seq data
#

source('prsem-util.R')

main <- function() {
  name2func <- list(
    'guessFqEncoding' = guessFqEncoding,
    'alignReads'      = alignReads
  )

  argv <- commandArgs(trailingOnly=T)
  name2func[[argv[1]]](argv[2:length(argv)])
}


alignReads <- function(argv) {
  library(data.table)

  nthr      <- strtoi(argv[1])
  s_infiles <- argv[2]
  fencod    <- argv[3]
  refname   <- argv[4]
  imdname   <- argv[5]
  bowtie    <- argv[6]
  samtools  <- argv[7]
  bedtools  <- argv[8]

  ffqs <- strsplit(s_infiles, ',', fixed=T)[[1]]
  encoddt <- fread(fencod, header=T, sep="\t")
  nthr_bowtie <- ifelse(nthr > 4, nthr-4, 1) 
  bowtie_ref <- paste0(refname, '_prsem')

  for ( ffq in ffqs ) {
    id <- Util$getFileNameSansExt(ffq)
    fout <- paste0(imdname, '_prsem.', id, '.tagAlign.gz')
    encod <- subset(encoddt, file == ffq)[, encoding]
    cmd_cat <- ifelse(Util$checkIfGzipByExt(ffq), 'zcat', 'cat')
    cmd <- paste0(cmd_cat, ' ', ffq , ' | ', 
                  bowtie, ' -q -v 2 -a --best --strata -m 1 ', encod, 
                          ' -S -p ', nthr_bowtie, ' ', bowtie_ref, ' - |',
                  samtools, ' view -S -b -F 1548 - | ',
                  bedtools, ' bamtobed -i stdin | ',
                  quote(`awk 'BEGIN{FS=\"\t";OFS="\t"}{$4="N"; print $0}' |`),
                  'gzip -c > ', fout)

    cat("\n", cmd, "\n")
    system(cmd)
  }
}


guessFqEncoding <- function(argv){
  Util$checkInstallCRAN('data.table')
  Util$checkInstallBioc('ShortRead')

  library(data.table)
  library(ShortRead)

  nthr      <- strtoi(argv[1])
  s_infiles <- argv[2]
  fout      <- argv[3]

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
