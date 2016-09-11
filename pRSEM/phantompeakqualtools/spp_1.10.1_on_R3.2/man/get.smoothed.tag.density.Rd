\name{get.smoothed.tag.density}
\alias{get.smoothed.tag.density}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate chromosome-wide profiles of smoothed tag density }
\description{
  Given tag positions, the method calculates for each chromosome a tag
  density profile, smoothed by the Gaussian kernel. If the optional
  control tags are provided, the difference between ChIP and control tag
  density is returned.
}
\usage{
get.smoothed.tag.density(signal.tags, control.tags = NULL, bandwidth = 150, bg.weight = NULL, tag.shift = 146/2, step = round(bandwidth/3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signal.tags}{ signal chromosome tag coordinate vectors (e.g. output
    of \code{\link{select.informative.tags}} }
  \item{control.tags}{ optional control (input) tags }
  \item{bandwidth}{ standard deviation of the Gaussian kernel }
  \item{bg.weight}{ optional weight by which the background density
    should be multipled for scaling. If not supplied, the weight is
    calculated based on the ratio of the reduced ChIP to input dataset sizes. }
  \item{tag.shift}{ Distance by which the positive and negative strand
    tags should be shifted towards eachother. This
    normally corresponds to the half of the cross-correlation peak
    position (e.g. \code{get.binding.characteristics()}$peak$x/2) }
  \item{step}{ The distance between the regularly spaced points for
    which the values should be calculated. }
  \item{background.density.scaling}{ If TRUE, regions of significant tag
    enrichment will be masked out when calculating size ratio of the
    signal to control datasets (to estimate ratio of the background tag
    density). If FALSE, the dataset ratio will be equal to the ratio of
    the number of tags in each dataset.}
}
\value{
  A list of elements corresponding to chromosomes, with each element
  being an $x/$y data.frame giving the position and associated tag
  density. Use \code{\link{writewig}} to output the structure to a WIG
  file.
}
\seealso{ \code{\link{writewig}} }
\examples{
  smoothed.density <- get.smoothed.tag.density(chip.data,control.tags=input.data,bandwidth=200,step=100,tag.shift=round(binding.characteristics$peak$x/2));
  writewig(smoothed.density,"example.density.wig","Example smoothed, background-subtracted tag density");
}