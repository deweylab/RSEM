# 1-20-10 Qunhua Li
#
# This program first plots correspondence curve and IDR threshold plot
# (i.e. number of selected peaks vs IDR) for each pair of sample
#
# It then performs consistency analysis on merged data
# It takes the parameters estimated from pairwise consistency analysis, and
# use the same parameters to determine threshold on the merged data
#
# usage: 
# Rscript batch-consistency-plot-merged.r [npairs] [output.dir] [input.file.prefix 1, 2, 3 ...] [half.width] [idr.level] [significant measure] [pooled.filename] [write out option] [overlap.ratio] [is.broadpeak]
# [npairs]: integer, number of consistency analyses
#          (e.g. if 2 replicates, npairs=1, if 3 replicates, npairs=3
# [output.dir]: output directory for plot
# [input.file.prefix 1, 2, 3]: prefix for the output from batch-consistency-analysis2. They are the input files for merged analysis see below for examples (i.e. saved.file.prefix). It can be multiple files
#
# The parameters below are for processing merged data
# [half.width]: -1 if using reported interval
# [idr.level]: threshold for idr
# [significant measure]: choose from "p.value", "q.value" or "signal.value"
# [pooled.filename]: peak caller output in narrowpeak or broadpeak format
# [write out option]: logical, T: write out selected peaks in merged data, F: not write out 
# [overlap.ratio]: minimum overlap for two peaks to be called as calling the
#  same region. A numerical value between 0 and 1. If 0, minimum overlap
#  is >=1bp.
# [is.broadpeak]:  a logical value. If broadpeak is used, set as T;
# if narrowpeak is used, set as F

args <- commandArgs(trailingOnly=T)

npair <- args[1] # number of curves to plot on the same figure
output.file.prefix <- args[2] # file name for plot, generated from script at the outer level

df.txt <- 10

## examples for debugging
#npair <- 3
#output.file.prefix <- "~/ENCODE/anshul/results/gm12878-cfos-YALE-combined-threshold/consistency-plot" 
#combofile <- "~/ENCODE/anshul/data/SPP.YaleGm12878Cfos.VS.Gm12878Input.PointPeak.narrowPeak"
#saved.file.prefix <- list()
#saved.file.prefix[[1]] <- "~/ENCODE/anshul/results/gm12878-cfos-YALE-combined-threshold/SPP.YaleRep1Gm12878Cfos.VS.SPP.YaleRep2Gm12878Cfos"
#saved.file.prefix[[2]] <- "~/ENCODE/anshul/results/gm12878-cfos-YALE-combined-threshold/SPP.YaleRep1Gm12878Cfos.VS.SPP.YaleRep3Gm12878Cfos"
#saved.file.prefix[[3]] <- "~/ENCODE/anshul/results/gm12878-cfos-YALE-combined-threshold/SPP.YaleRep2Gm12878Cfos.VS.SPP.YaleRep3Gm12878Cfos"

#npair <- 1
#output.file.prefix <- "~/ENCODE/anshul/results/gm12878-pol2-YALE-combined-threshold/consistency-plot" 
#combofile <- "~/ENCODE/anshul/data/SPP.YaleGm12878Pol2.VS.Gm12878Input.PointPeak.narrowPeak"
#saved.file.prefix <- "~/ENCODE/anshul/results/gm12878-pol2-YALE-combined-threshold/SPP.YaleRep1Gm12878Pol2.VS.SPP.YaleRep2Gm12878Pol2"


#ori.sig.value <- "signal.value"
# nominal.sig.value <- "q.value"
# idr.level <- 0.05
# half.width <- NULL
###################

# the df for plotting the smooth spline on the consistency curve
#if(length(args)-3> npair){ # if df is specified
#  df.txt <- as.numeric(args[length(args)]) # df for plotting, default is 10
#}else{ 
#  df.txt <- 10
#}

ntemp <- as.numeric(npair)

###### this is needed for pooled data

cat(as.numeric(args[3+ntemp]))
if(as.numeric(args[3+ntemp])==-1){ # enter -1 when using the reported length 
  half.width <- NULL
}else{
  half.width <- as.numeric(args[3+ntemp])
}


idr.level <- as.numeric(args[4+ntemp])  # this is the consistency FDR, e.g. 0.05
# a string: "signal.value", "p.value" or "q.value", for specifying which
# significant value to use for thresholding the merged data
ori.sig.value <- args[5+ntemp] 


# pooled data file
combofile <- args[6+ntemp]
is.write.out <- as.logical(args[7+ntemp])
overlap.ratio <- as.numeric(args[8+ntemp]) # the minimum amount of overlap to be called as an overlap

is.broadpeak <- args[9+ntemp]

saved.file.prefix <- list() # identifier of filenames that contain the em and URI results


source("functions-all-clayton-12-13.r")

uri.list <- list()
uri.list.match <- list()
ez.list <- list()
legend.txt <- c()
#fdr.map <- c()
sig.map <- list() 
em.output.list <- list()
uri.output.list <- list()

for(i in 1:npair){
  saved.file.prefix[i] <- args[2+i]
 
  load(paste(saved.file.prefix[i], "-uri.sav", sep=""))
  load(paste(saved.file.prefix[i], "-em.sav", sep=""))

  uri.output.list[[i]] <- uri.output
  em.output.list[[i]] <- em.output

  ez.list[[i]] <- get.ez.tt.all(em.output, uri.output.list[[i]]$data12.enrich$merge1,
                                uri.output.list[[i]]$data12.enrich$merge2, idr.level=idr.level) # reverse =T for error rate

  # URI for all peaks
  uri.list[[i]] <- uri.output$uri.n
  # URI for matched peaks
  uri.match <- get.uri.matched(em.output$data.pruned, df=df.txt)
  uri.list.match[[i]] <- uri.match$uri.n

  file.name <- unlist(strsplit(as.character(saved.file.prefix[i]), "/"))
  
  legend.txt[i] <- paste(i, "=", file.name[length(file.name)])
  sig.map[[i]] <- cbind(idr.level, ez.list[[i]]$map.uv)


  # map idr computed from consistency back to the original significant measure
  # 
  # if(is.null(nominal.sig.value)){
  #  sig.map <- cbind(idr.level, ez.list[[i]]$map.uv)
  #} else {

    # for SPP, need find the significant value based on FDR
  #  temp.map <- map.sig.value(uri.output$data12.enrich, ez.list[[i]]$map.uv, nominal.value=nominal.sig.value)
  #  sig.map <- cbind(idr.level, ez.list[[i]]$map.uv)
    # this is the corresponding FDR mapped from the significant value
    # you don't need this in general
  #  fdr.map <- cbind(idr.level, temp.map)
  #}

}

plot.uri.file <- paste(output.file.prefix, "-plot.ps", sep="")

cat("plot consistency plots\n")
############# plot and report output
# plot correspondence curve for each pair,
# plot number of selected peaks vs IDR 

# plot all into 1 file
postscript(paste(output.file.prefix, "-plot.ps", sep=""))
par(mfcol=c(2,3), mar=c(5,6,4,2)+0.1)
plot.uri.group(uri.list, NULL, file.name=NULL, c(1:npair), title.txt="all peaks")
plot.uri.group(uri.list.match, NULL, file.name=NULL, c(1:npair), title.txt="matched peaks")
plot.ez.group(ez.list, plot.dir=NULL, file.name=NULL, legend.txt=c(1:npair), y.lim=c(0, 0.6))
plot(0, 1, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n") # legends
legend(0, 1, legend.txt, cex=0.6)

dev.off()

############### consistency cutoff on the replicates #############

cat("read pooled sample \n")
##################################################
########## now this part is for combined dataset
##################################################

chr.file <- "genome_table.txt"

chr.size <- read.table(chr.file)

# read combined data
combined.ori <- process.narrowpeak(paste(combofile, sep=""), chr.size, half.width=half.width, summit="offset", broadpeak=is.broadpeak)$data.cleaned

#combined <- combined.ori[, c("chr", "start", "stop", ori.sig.value, "signal.value", "p.value", "q.value", "start.chr", "stop.chr")]
#colnames(combined) <- c("chr", "start", "stop", "sig.value", "signal.value", "p.value", "q.value",  "start.chr", "stop.chr")

combined <- combined.ori[, c( ori.sig.value, "start", "stop","signal.value", "p.value", "q.value", "chr", "start.ori", "stop.ori")]
colnames(combined) <- c("sig.value", "start", "stop",  "signal.value", "p.value", "q.value", "chr", "start.ori", "stop.ori")
combined$frac.ratio <- NA

########
#  map by the matched structure
########
cat("Selecting peaks using parameters from consistency analysis\n")
sig.select.method2 <- pass.structure(uri.output.list, em.output.list, combined, idr.level=idr.level, sig.value.impute=0, chr.size)

if(is.write.out){
  write.table(sig.select.method2$combined.selected, file=paste(output.file.prefix, "-combined.selection.txt", sep=""), quote=F, row.names=F)
}

save(sig.select.method2, file=paste(output.file.prefix, "-select.sav", sep=""))


# output for ez
sink(paste(output.file.prefix, "-Rout.txt", sep=""))
cat("IDR Map for specified sig.value", "\n")
print(sig.map)

# cat("IDR Map for", nominal.sig.value, "\n")
# print(fdr.map)


# output for merged dataset
cat("Merged dataset has ", nrow(combined), "p\n")

cat("Apply parameters estimated from consistency analysis to merged data: select by ", ori.sig.value, "\n")
print(sig.select.method2$npeak.stat)
cat("Range of significant values on the selected pooled data", "\n")
print(sig.select.method2$sig.combined)

sink()

