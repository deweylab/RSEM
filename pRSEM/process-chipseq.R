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
}


guessFqEncoding <- function(argv){
  .libPaths(c(argv[4], .libPaths()))
  suppressMessages(library(data.table))
  suppressMessages(library(ShortRead))

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
 #cat('File written:', fout, "\n")
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


main()
