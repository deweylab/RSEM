\name{add.broad.peak.regions}
\alias{add.broad.peak.regions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate chromosome-wide profiles of smoothed tag density }
\description{
  Looks for broader regions of enrichment associated with the determined
  peak positions, adds them to the $npl data as $rs, $re columns. 
}
\usage{
add.broad.peak.regions(signal.tags, control.tags, binding.postions,window.size=500,z.thr=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signal.tags}{ signal chromosome tag coordinate vectors (e.g. output
    of \code{\link{select.informative.tags}} }
  \item{control.tags}{ optionall control (input) tags }
  \item{binding.positions}{ output of find.binding.positions call }
  \item{window.size}{ window size to be used in calculating enrichment }
  \item{z.thr}{ Z-score corresponding to the Poisson ratio threshold
    used to flag significantly enriched windows}
}
\value{
  A structure identical to binding.postions with two additional columns
  added (rs and re) corresponding to start and end of the associated
  significantly enriched region. If no region was associated with a
  particular peak, NAs values are reported.
}
