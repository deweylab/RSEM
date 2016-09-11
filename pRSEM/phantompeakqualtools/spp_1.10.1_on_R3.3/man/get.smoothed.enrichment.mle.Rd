\name{get.smoothed.enrichment.mle}
\alias{get.smoothed.enrichment.mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate chromosome-wide profiles of smoothed enrichment estimate }
\description{
  Given signal and control tag positions, the method calculates log2
  signal to control enrichment esimates (maximum likelihood) for each
  chromosome, based on the smoothed tag density profile (see \link{get.smoothed.tag.density}).
}
\usage{
get.smoothed.enrichment.mle(signal.tags, control.tags, bandwidth = 150,tag.shift = 146/2, step = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signal.tags}{ signal chromosome tag coordinate vectors (e.g. output
    of \code{\link{select.informative.tags}} }
  \item{control.tags}{ control (input) tags }
  \item{pseudocount}{ pseudocount value to be added to tag density -
    defaults to 1 }
  other parameters (such as bandwidth, step.size and tag.shift) are
  passed to \link{get.smoothed.tag.density} - see appropriate reference
  for details.
}
\value{
  A list of elements corresponding to chromosomes, with each element
  being an $x/$y data.frame giving the position and associated
  log2 signal/control enrichment estimate. 
}
\seealso{ \code{\link{writewig}} }
\examples{
  # get smoothed enrichment estimate profile using 500bp bandwidth at
  # 50bp steps
  smoothed.M <- get.smoothed.enrichment.mle(chip.data,bandwidth=500,step=50);
  writewig(smoothed.M,"example.smoothedM.wig","Example smoothed log2 intensity ratio estimate");
}