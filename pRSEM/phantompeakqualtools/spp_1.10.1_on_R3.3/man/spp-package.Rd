\name{spp-package}
\alias{spp-package}
\alias{spp}
\docType{package}
\title{
ChIP-seq (Solexa) Processing Pipeline
}
\description{
A set of routines for reading short sequence alignments, calculating tag
density, estimates of statistically significant enrichment/depletion
along the chromosome, identifying point binding positions (peaks), and
characterizing saturation properties related to sequencing depth.
}
\details{
\tabular{ll}{
Package: \tab spp\cr
Type: \tab Package\cr
Version: \tab 1.8\cr
Date: \tab 2008-11-14\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
See example below for typical processing sequence.y
}
\author{Peter Kharchenko <peter.kharchenko@post.harvard.edu>}
\references{
Kharchenko P., Tolstorukov M., Park P. "Design and analysis of ChIP-seq
experiments for DNA-binding proteins." Nature Biotech. doi:10.1038/nbt.1508
}

\examples{

  # load the library
  library(spp);

  ## The following section shows how to initialize a cluster of 8 nodes for parallel processing
  ## To enable parallel processing, uncomment the next three lines, and comment out "cluster<-NULL";
  ## see "snow" package manual for details.
  #library(snow)
  #cluster <- makeCluster(2);
  #invisible(clusterCall(cluster,source,"routines.r"));
  cluster <- NULL;


  
  # read in tag alignments
  chip.data <- read.eland.tags("chip.eland.alignment");
  input.data <- read.eland.tags("input.eland.alignment");
  
  # get binding info from cross-correlation profile
  # srange gives the possible range for the size of the protected region;
  # srange should be higher than tag length; making the upper boundary too high will increase calculation time
  #
  # bin - bin tags within the specified number of basepairs to speed up calculation;
  # increasing bin size decreases the accuracy of the determined parameters
  binding.characteristics <- get.binding.characteristics(chip.data,srange=c(50,500),bin=5,cluster=cluster);


  # plot cross-correlation profile
  pdf(file="example.crosscorrelation.pdf",width=5,height=5)
  par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8);
  plot(binding.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation");
  abline(v=binding.characteristics$peak$x,lty=2,col=2)
  dev.off();
  
  # select informative tags based on the binding characteristics
  chip.data <- select.informative.tags(chip.data,binding.characteristics);
  input.data <- select.informative.tags(input.data,binding.characteristics);

  # restrict or remove positions with anomalous number of tags relative
  # to the local density
  chip.data <- remove.local.tag.anomalies(chip.data);
  input.data <- remove.local.tag.anomalies(input.data);


  # output smoothed tag density (subtracting re-scaled input) into a WIG file
  # note that the tags are shifted by half of the peak separation distance
  smoothed.density <- get.smoothed.tag.density(chip.data,control.tags=input.data,bandwidth=200,step=100,tag.shift=round(binding.characteristics$peak$x/2));
  writewig(smoothed.density,"example.density.wig","Example smoothed, background-subtracted tag density");
  rm(smoothed.density);

  # output conservative enrichment estimates
  # alpha specifies significance level at which confidence intervals will be estimated
  enrichment.estimates <- get.conservative.fold.enrichment.profile(chip.data,input.data,fws=2*binding.characteristics$whs,step=100,alpha=0.01);
  writewig(enrichment.estimates,"example.enrichment.estimates.wig","Example conservative fold-enrichment/depletion estimates shown on log2 scale");
  rm(enrichment.estimates);


  # binding detection parameters
  # desired FDR. Alternatively, an E-value can be supplied to the method calls below instead of the fdr parameter
  fdr <- 1e-2; 
  # the binding.characteristics contains the optimized half-size for binding detection window
  detection.window.halfsize <- binding.characteristics$whs;
  
  # determine binding positions using wtd method
  bp <- find.binding.positions(signal.data=chip.data,control.data=input.data,fdr=fdr,method=tag.wtd,whs=detection.window.halfsize,cluster=cluster)

  # alternatively determined binding positions using lwcc method (note: this takes longer than wtd)
  # bp <- find.binding.positions(signal.data=chip.data,control.data=input.data,fdr=fdr,method=tag.lwcc,whs=detection.window.halfsize,cluster=cluster)
  
  print(paste("detected",sum(unlist(lapply(bp$npl,function(d) length(d$x)))),"peaks"));
  
  # output detected binding positions
  output.binding.results(bp,"example.binding.positions.txt");


  # ------------------------------------------------------------------------------------------- 
  # the set of commands in the following section illustrates methods for saturation analysis
  # these are separated from the previous section, since they are highly CPU intensive
  # ------------------------------------------------------------------------------------------- 

  # determine MSER
  # note: this will take approximately 10-15x the amount of time the initial binding detection did
  # The saturation criteria here is 0.99 consistency in the set of binding positions when adding 1e5 tags.
  # To ensure convergence the number of subsampled chains (n.chains) should be higher (80)
  mser <- get.mser(chip.data,input.data,step.size=1e5,test.agreement=0.99,n.chains=8,cluster=cluster,fdr=fdr,method=tag.wtd,whs=detection.window.halfsize)
  
  print(paste("MSER at a current depth is",mser));
  
  # note: an MSER value of 1 or very near one implies that the set of detected binding positions satisfies saturation criteria without
  # additional selection by fold-enrichment ratios. In other words, the dataset has reached saturation in a traditional sense (absolute saturation).

  # interpolate MSER dependency on tag count
  # note: this requires considerably more calculations than the previous steps (~ 3x more than the first MSER calculation)
  # Here we interpolate MSER dependency to determine a point at which MSER of 2 is reached
  # The interpolation will be based on the difference in MSER at the current depth, and a depth at 5e5 fewer tags (n.steps=6);
  # evaluation of the intermediate points is omitted here to speed up the calculation (excluded.steps parameter)
  # A total of 7 chains is used here to speed up calculation, whereas a higher number of chains (50) would give good convergence
  msers <- get.mser.interpolation(chip.data,input.data,step.size=1e5,test.agreement=0.99, target.fold.enrichment=2, n.chains=7,n.steps=6,excluded.steps=c(2:4),cluster=cluster,fdr=fdr,method=tag.wtd,whs=detection.window.halfsize)

  print(paste("predicted sequencing depth =",round(unlist(lapply(msers,function(x) x$prediction))/1e6,5)," million tags"))
  

  # note: the interpolation will return NA prediction if the dataset has reached absolute saturation at the current depth.
  # note: use return.chains=T to also calculated random chains (returned under msers$chains field) - these can be passed back as
  #       "get.mser.interpolation( ..., chains=msers$chains)" to calculate predictions for another target.fold.enrichment value
  #        without having to recalculate the random chain predictions.

  ## stop cluster if it was initialized
  #stopCluster(cluster);  



}
