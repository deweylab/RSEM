\name{get.mser}
\alias{get.mser}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate minimal saturated enrichment fold ratio }
\description{
  Determine if the dataset has reached absolute saturation, or otherwise
  find minimal fold enrichment ratio above which the detection of peaks
  has stabilized enough to meet the saturation criteria.
}
\usage{
get.mser(signal.data, control.data, n.chains = 5, step.size = 1e+05, chains = NULL, cluster = NULL, test.agreement = 0.99, return.chains = F, enrichment.background.scales = c(1), n.steps = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signal.data}{ signal tag vector list }
  \item{control.data}{ control tag vector list }
  \item{n.chains}{ number of dataset subsamples to use }
  \item{step.size}{ subsampling step describing the saturation
    criteria. The criteria requires the set of detected binding sites to
    be stable (as described by the \code{test.agreement} param) when the
    number of tags in the dataset is reduced by \code{step.size}. The
    value can either be an integer above one, in which case it specifies a fixed
    number of tags, or a real value below one, in which case it
    specifies the fraction of tags that should be removed (e.g. 0.1 will
    remove 10% of tags).
  }
  \item{test.agreement}{ Fraction of the detected peaks that should
    agree between the full and subsampled datasets. }
  \item{chains}{ optional parameter, giving pre-calculated chains }
  \item{cluster}{ optional \code{snow} cluster to parallelize processing }

  \item{return.chains}{ whether subsampled dataset results should be returned as
  well }
  \item{enrichment.background.scales}{ one or multiple window scales at
  which the background tag density should be assessed. See
  \code{enrichment.background.scales} in
  \code{\link{find.binding.positions}}. If multiple scales are provided,
  multiple MSER estimates will be returned.}
  \item{\dots}{ additional parameters should be the same as those passed
  to the \code{\link{find.binding.positions}}}
}
\value{
  A single, or multple (if multiple \code{enrichment.background.scales} were
  provided) MSER value. A value of 1 or very close to it implies that
  the dataset has reached absolute saturation based on the given criteria.
}
