\name{get.binding.characteristics}
\alias{get.binding.characteristics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate characteristics of observed DNA-binding signal from
  cross-correlation profiles }
\description{
  The methods calculates strand cross-correlation profile to determine binding
  peak separation distance and approximate window size that should be used
  for binding detection. If quality scores were given for the tags,
  which quality bins improve the cross-correlation pattern.
}
\usage{
get.binding.characteristics(data, srange = c(50, 500), bin = 5, cluster = NULL, debug = F, min.tag.count = 1000, acceptance.z.score = 3, remove.tag.anomalies = T, anomalies.z = 5,accept.all.tags=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Tag/quality data: output of \code{read.eland.tags} or similar function }
  \item{srange}{ A range within which the binding peak separation is
    expected to fall. Should be larger than probe size to avoid artifacts. }
  \item{bin}{ Resolution (in basepairs) at which cross-corrrelation
    should be calculated. bin=1 is ideal, but takes longer to calculate. }
  \item{cluster}{ optional snow cluster for parallel processing }
  \item{debug}{ whether to print debug messages }
  \item{min.tag.count}{ minimal number of tags on the chromosome to be
    considered in the cross-correlation calculations }
  \item{acceptance.z.score}{ A Z-score used to determine if a given tag
    quality bin provides significant improvement to the strand cross-correlation }
  \item{remove.tag.anomalies}{ Whether to remove singular tag count peaks prior to
  calculation. This is recommended, since such positions may distort the
  cross-correlation profile and increase the necessary computational time. }
  \item{anomalies.z}{ Z-score for determining if the number of tags at a
  given position is significantly higher about background, and should be
  considered an anomaly.}
  \item{accept.all.tags}{ Whether tag alignment quality calculations
    should be skipped and all available tags should be accepted in the
    downstream analysis.}
}
\value{
  \item{cross.correlation }{ Cross-correlation profile as an $x/$y data.frame}
  \item{peak }{Position ($x) and height ($y) of automatically detected
    cross-correlation peak.}
  \item{whs}{ Optimized window half-size for binding detection (based
    on the width of the cross-correlation peak) }
  \item{quality.bin.acceptance}{ A list structure, describing the
    effect of inclusion of different tag quality bins on
    cross-correlation, and a resolution on which bins should be
    considered.
    \itemize{
      \item{informative.bins}{ A boolean vector indicating whether the
        inclusion of tags from the tag quality bin specified in the name
        attribute significantly increases cross-correlation profile near
        the peak.}
      \item{quality.cc}{ A list giving the cross-correlation profile
        after the inclusion of the tags from different quality bins 
      }
    }
  }
}
