\name{find.binding.positions}
\alias{find.binding.positions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Determine significant point protein binding positions (peaks) }
\description{
  Given the signal and optional control (input) data, determine location of the
  statistically significant point binding positions. If the control data
  is not provided, the statistical significance can be assessed based on
  tag randomization. The method also provides options for masking
  regions exhibiting strong signals within the control data.
}
\usage{
find.binding.positions(signal.data, e.value = NULL, fdr = NULL, masked.data = NULL, control.data = NULL, min.dist = 200, window.size = 4e+07, cluster = NULL, debug = T, n.randomizations = 3, shuffle.window = 1, min.thr = 0, topN = NULL, tag.count.whs = 100, enrichment.z = 2, method = tag.wtd, tec.filter = T, tec.window.size = 10000, tec.masking.window.size=tec.window.size, tec.z = 5, tec.poisson.z=5,tec.poisson.ratio=5, n.control.samples = 1, enrichment.background.scales = c(1, 5, 10), background.density.scaling = F, use.randomized.controls = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  ~~ tag data ~~
  \item{signal.data}{ signal tag vector list }
  \item{control.data}{ optional control (input) tag vector list }

  ~~ position stringency criteria ~~
  \item{e.value}{ E-value defining the desired statistical significance
    of binding positions.  }
  \item{fdr}{ FDR defining statistical significance of binding positions  }
  \item{topN}{ instead of determining statistical significance
    thresholds, return the specified number of highest-scoring
    positions}

  ~~ other params ~~
  \item{whs}{ window half-sized that should be used for binding
    detection (e.g. determined from cross-correlation profiles)}
  \item{masked.data}{ optional set of coordinates that should be masked
    (e.g. known non-unique regions) }
  \item{min.dist}{ minimal distance that must separate detected binding
    positions. In case multiple binding positions are detected within
    such distance, the position with the highest score is returned. }
  \item{window.size}{ size of the window used to segment the chromosome
    during calculations to reduce memory usage. }
  \item{cluster}{ optional \code{snow} cluster to parallelize the
    processing on }
  \item{min.thr}{ minimal score requirement for a peak }
  \item{background.density.scaling}{ If TRUE, regions of significant tag
    enrichment will be masked out when calculating size ratio of the
    signal to control datasets (to estimate ratio of the background tag
    density). If FALSE, the dataset ratio will be equal to the ratio of
    the number of tags in each dataset.}

  ~~ randomized controls ~~
  \item{n.randomizations}{ number of tag randomziations that should be
    performed (when the control data is not provided) }
  \item{use.randomized.controls}{ Use randomized tag control, even if
    \code{control.data} is supplied. }
  \item{shuffle.window}{ during tag randomizations, tags will be split
    into groups of \code{shuffle.window} and will be maintained
    together throughout the randomization. }

  ~~ fold-enrichment confidence intervals
  \item{tag.count.whs}{ half-size of a window used to assess fold
    enrichment of a binding position}
  \item{enrichment.z}{ Z-score used to define the significance level of
    the fold-enrichment confidence intervals }
    \item{enrichment.background.scales}{ In estimating the peak
      fold-enrichment confidence intervals, the background tag density is
      estimated based on windows with half-sizes of
  \code{2*tag.count.whs*enrichment.background.scales}. }
  \item{method}{ either \code{tag.wtd} for WTD method, or
    \code{tag.lwcc} for MTC method}
  \item{mle.filter}{ If turned on, will exclude predicted positions
    whose MLE enrichment ratio (for any of the background scales) is
    below a specified min.mle.threshold }
  \item{min.mle.threshold}{ MLE enrichment ratio threshold that each
    predicted position must exceed if mle.filter is turned on. }

  ~~ masking regions of significant control enrichment ~~
  \item{tec.filter}{ Whether to mask out the regions exhibiting
  significant enrichment in the control data in doing other
  calculations. The regions are identified using Poisson statistics
  within sliding windows, either relative to the scaled signal (tec.z), or
  relative to randomly-distributed expectation (tec.poisson.z).}
  \item{tec.window.size}{ size of the window used to determine
    significantly enrichent control regions }
  \item{tec.masking.window.size}{ size of the window used to mask 
    the area around significantly enrichent control regions }
  \item{tec.z}{ Z-score defining statistical stringency by which a given
    window is determined to be significantly higher in the input than in
    the signal, and masked if that is the case.}
  \item{tec.poisson.z}{ Z-score defining statistical stringency by which a given
    window is determined to be significantly higher than the
    tec.poisson.ratio above the expected uniform input background. }
  \item{tec.poisson.ratio}{ Fold ratio by which input must exceed the
    level expected from the uniform distribution. }


  
  
}
\value{
  \item{npl}{A per-chromosome list containing data frames describing
    determined binding positions. Column description:
    \itemize{
      \item{x} {position}
      \item{y}{ score }
      \item{evalue}{ E-value }
      \item{fdr}{ FDR. For peaks higher than the maximum control peak,
                  the highest dataset FDR is reported }
      \item{enr}{ lower bound of the fold-enrichment ratio confidence
                  interval. This is the estimate determined using scale of
                  1. Estimates corresponding to higher scales are returned in other enr columns
                  with scale appearing in the name.}
      \item{enr.mle}{ enrichment ratio maximum likely estimate }
    }
  }
  \item{thr}{ info on the chosen statistical threshold of the peak scores}
}

\examples{
  # find binding positions using WTD method, 200bp half-window size,
control data, 1% FDR
  bp <-
find.binding.positions(signal.data=chip.data,control.data=input.data,fdr=0.01,method=tag.wtd,whs=200);

  # find binding positions using MTC method, using 5 tag randomizations,
  #  keeping pairs of tag positions together (shuffle.window=2)
  bp <- find.binding.positions(signal.data=chip.data,control.data=input.data,fdr=0.01,method=tag.lwcc,whs=200,use.randomized.controls=T,n.randomizations=5,shuffle.window=2)

  # print out the number of determined positions  
  print(paste("detected",sum(unlist(lapply(bp$npl,function(d) length(d$x)))),"peaks"));
  

}
