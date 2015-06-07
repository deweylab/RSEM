# modified 06/07/2015 pliu
#
# modified 3-29-10: Qunhua Li
# add 2 columns in the output of "-overlapped-peaks.txt": local.idr and IDR

# 01-20-2010 Qunhua Li
#
# This program performs consistency analysis for a pair of peak calling outputs
# It takes narrowPeak or broadPeak formats.
# 
# usage: Rscript batch-consistency-analysis2.r peakfile1 peakfile2 half.width outfile.prefix overlap.ratio  is.broadpeak sig.value
#
# peakfile1 and peakfile2 : the output from peak callers in narrowPeak or broadPeak format
# half.width: -1 if using the reported peak width, 
#             a numerical value to truncate the peaks to
# outfile.prefix: prefix of output file
# overlap.ratio: a value between 0 and 1. It controls how much overlaps two peaks need to have to be called as calling the same region. It is the ratio of overlap / short peak of the two. When setting at 0, it means as long as overlapped width >=1bp, two peaks are deemed as calling the same region.
# is.broadpeak: a logical value. If broadpeak is used, set as T; if narrowpeak is used, set as F
# sig.value: type of significant values, "q.value", "p.value" or "signal.value" (default, i.e. fold of enrichment)

args <- commandArgs(trailingOnly=T)

# consistency between peakfile1 and peakfile2
#input1.dir <- args[1]
#input2.dir <- args[2] # directories of the two input files
peakfile1 <- args[1]
peakfile2 <- args[2]

if(as.numeric(args[3])==-1){ # enter -1 when using the reported length 
  half.width <- NULL
}else{
  half.width <- as.numeric(args[3])
}


output.prefix <- args[4]
overlap.ratio <- args[5]

if(args[6] == "T"){
  is.broadpeak <- T
}else{
  is.broadpeak <- F
}

sig.value <- args[7]


#dir1 <- "~/ENCODE/anshul/data/"
#dir2 <- dir1
#peakfile1 <- "../data/SPP.YaleRep1Gm12878Cfos.VS.Gm12878Input.PointPeak.narrowPeak"
#peakfile2 <- "../data/SPP.YaleRep3Gm12878Cfos.VS.Gm12878Input.PointPeak.narrowPeak"
#half.width <- NULL
#overlap.ratio <- 0.1
#sig.value <- "signal.value"

## pliu
idr_code_path <- args[8]
source(paste0(idr_code_path, "/functions-all-clayton-12-13.r"))
#source("functions-all-clayton-12-13.r")
##

# read the length of the chromosomes, which will be used to concatenate chr's
## pliu 20150607
#chr.file <- "genome_table.txt"
chr.file <- args[9]
##

chr.size <- read.table(chr.file)


sink(paste(output.prefix, "-Rout.txt", sep=""))

############# process the data
cat("is.broadpeak", is.broadpeak, "\n")
# process data, summit: the representation of the location of summit
rep1 <- process.narrowpeak(paste(peakfile1, sep=""), chr.size, half.width=half.width, summit="offset", broadpeak=is.broadpeak)
rep2 <- process.narrowpeak(paste(peakfile2, sep=""), chr.size, half.width=half.width, summit="offset", broadpeak=is.broadpeak)

cat(paste("read", peakfile1, ": ", nrow(rep1$data.ori), "peaks\n", nrow(rep1$data.cleaned), "peaks are left after cleaning\n", peakfile2, ": ", nrow(rep2$data.ori), "peaks\n", nrow(rep2$data.cleaned), " peaks are left after cleaning"))

if(args[3]==-1){
  cat(paste("half.width=", "reported", "\n"))
}else{
  cat(paste("half.width=", half.width, "\n"))
}  
cat(paste("significant measure=", sig.value, "\n"))

# compute correspondence profile (URI)
uri.output <- compute.pair.uri(rep1$data.cleaned, rep2$data.cleaned, sig.value1=sig.value, sig.value2=sig.value, overlap.ratio=overlap.ratio)

#uri.output <- compute.pair.uri(rep1$data.cleaned, rep2$data.cleaned)

cat(paste("URI is done\n"))

# save output
save(uri.output, file=paste(output.prefix, "-uri.sav", sep=""))
cat(paste("URI is saved at: ", output.prefix, "-uri.sav \n", sep=""))


# EM procedure for inference
em.output <- fit.em(uri.output$data12.enrich, fix.rho2=T)

#em.output <- fit.2copula.em(uri.output$data12.enrich, fix.rho2=T, "gaussian")

cat(paste("EM is done\n\n"))


save(em.output, file=paste(output.prefix, "-em.sav", sep=""))
cat(paste("EM is saved at: ", output.prefix, "-em.sav \n", sep=""))


# write em output into a file

cat(paste("EM estimation for the following files\n", peakfile1, "\n", peakfile2, "\n", sep=""))

print(em.output$em.fit$para)

# add on 3-29-10
# output both local idr and IDR
idr.local <- 1-em.output$em.fit$e.z
IDR <- c()
o <- order(idr.local)
IDR[o] <- cumsum(idr.local[o])/c(1:length(o))


write.out.data <- data.frame(chr1=em.output$data.pruned$sample1[, "chr"],
                    start1=em.output$data.pruned$sample1[, "start.ori"],
                    stop1=em.output$data.pruned$sample1[, "stop.ori"],
                    sig.value1=em.output$data.pruned$sample1[, "sig.value"],   
                    chr2=em.output$data.pruned$sample2[, "chr"],
                    start2=em.output$data.pruned$sample2[, "start.ori"],
                    stop2=em.output$data.pruned$sample2[, "stop.ori"],
                    sig.value2=em.output$data.pruned$sample2[, "sig.value"],
                    idr.local=1-em.output$em.fit$e.z, IDR=IDR)

write.table(write.out.data, file=paste(output.prefix, "-overlapped-peaks.txt", sep=""))
cat(paste("Write overlapped peaks and local idr to: ", output.prefix, "-overlapped-peaks.txt\n", sep=""))

# number of peaks passing IDR range (0.01-0.25)
IDR.cutoff <- seq(0.01, 0.25, by=0.01)
idr.o <- order(write.out.data$idr.local)
idr.ordered <- write.out.data$idr.local[idr.o]
IDR.sum <- cumsum(idr.ordered)/c(1:length(idr.ordered))

IDR.count <- c()
n.cutoff <- length(IDR.cutoff)
for(i in 1:n.cutoff){
  IDR.count[i] <- sum(IDR.sum <= IDR.cutoff[i])
}


# write the number of peaks passing various IDR range into a file
idr.cut <- data.frame(peakfile1, peakfile2, IDR.cutoff=IDR.cutoff, IDR.count=IDR.count)
write.table(idr.cut, file=paste(output.prefix,"-npeaks-aboveIDR.txt", sep=""), append=T, quote=F, row.names=F, col.names=F)
cat(paste("Write number of peaks above IDR cutoff [0.01, 0.25]: ","npeaks-aboveIDR.txt\n", sep=""))

mar.mean <- get.mar.mean(em.output$em.fit)

cat(paste("Marginal mean of two components:\n"))
print(mar.mean)

sink()


