# 1-20-10 Qunhua Li
#
# This program first plots correspondence curve and IDR threshold plot
# (i.e. number of selected peaks vs IDR) for each pair of sample
#
# usage: 
# Rscript batch-consistency-plot-merged.r [npairs] [output.dir] [input.file.prefix 1, 2, 3 ...]
# [npairs]: integer, number of consistency analyses
#          (e.g. if 2 replicates, npairs=1, if 3 replicates, npairs=3
# [output.dir]: output directory for plot
# [input.file.prefix 1, 2, 3]: prefix for the output from batch-consistency-analysis2. They are the input files for merged analysis see below for examples (i.e. saved.file.prefix). It can be multiple files
#

args <- commandArgs(trailingOnly=T)
npair <- args[1] # number of curves to plot on the same figure
output.file.prefix <- args[2] # file name for plot, generated from script at the outer level
df.txt <- 10
ntemp <- as.numeric(npair)
saved.file.prefix <- list() # identifier of filenames that contain the em and URI results
source("functions-all-clayton-12-13.r")

uri.list <- list()
uri.list.match <- list()
ez.list <- list()
legend.txt <- c()
em.output.list <- list()
uri.output.list <- list()

for(i in 1:npair){
  saved.file.prefix[i] <- args[2+i]
 
  load(paste(saved.file.prefix[i], "-uri.sav", sep=""))
  load(paste(saved.file.prefix[i], "-em.sav", sep=""))

  uri.output.list[[i]] <- uri.output
  em.output.list[[i]] <- em.output

  ez.list[[i]] <- get.ez.tt.all(em.output, uri.output.list[[i]]$data12.enrich$merge1,
                                uri.output.list[[i]]$data12.enrich$merge2) # reverse =T for error rate

  # URI for all peaks
  uri.list[[i]] <- uri.output$uri.n
  # URI for matched peaks
  uri.match <- get.uri.matched(em.output$data.pruned, df=df.txt)
  uri.list.match[[i]] <- uri.match$uri.n

  file.name <- unlist(strsplit(as.character(saved.file.prefix[i]), "/"))
  
  legend.txt[i] <- paste(i, "=", file.name[length(file.name)])

}

plot.uri.file <- paste(output.file.prefix, "-plot.ps", sep="")

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
