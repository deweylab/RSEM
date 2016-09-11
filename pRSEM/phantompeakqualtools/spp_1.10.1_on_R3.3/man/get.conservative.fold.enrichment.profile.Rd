\name{get.conservative.fold.enrichment.profile}
\alias{get.conservative.fold.enrichment.profile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate minimal fold enrichment/depletion along the chromosomes }
\description{
  The method provides a statistical assessment of enrichment/depletion
  along the chromosomes. To assess tag density enrichment/depletion, a
  sliding window of a specified size (\code{fws}) is used to calculate
  the density of the foreground tags (\code{ftl}). Multiple, typically
  larger windows are used to estimate background tag (\code{btl}) density around the
  same location. The densities are compared as ratios of two Poisson
  processes to estimate lower bound of foreground enrichment, or upper
  bound of foreground depletion. If multiple window sizes were used to
  estimate the background tag density, the most conservative one is
  chosen for each point. 
}
\usage{
get.conservative.fold.enrichment.profile(ftl, btl, fws, bwsl = c(1, 5, 25, 50) * fws, step = 50, tag.shift = 146/2, alpha = 0.05, use.most.informative.scale = F, quick.calculation = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ftl}{ foreground tag vector list }
  \item{btl}{ background tag vector list }
  \item{fws}{ foreground window size }
  \item{bwsl}{ background window scales. The size(s) of background windows
    will be \code{fws*bwsl}. }
  \item{step}{ spacing between positions at which the
    enrichment/depletion is evaluated }
  \item{tag.shift}{ number of basepairs by which positive and negative
    tag coordinates should be shifted towards eachother (half of binding
  peak separation distance)}
  \item{alpha}{ desired level of statistical significance }
  \item{use.most.informative.scale}{ for each position, instead of
    evaluating enrichment ratio bounds for all background window scales,
    choose the one with the highest observed density to speed up the calculations}
  \item{quick.calculation}{ Use square root transformation method
    instead of a Bayesian method. This speeds up the caclulation
    considerably and is turned on by default. }
  \item{background.density.scaling}{ If TRUE, regions of significant tag
    enrichment will be masked out when calculating size ratio of the
    signal to control datasets (to estimate ratio of the background tag
    density). If FALSE, the dataset ratio will be equal to the ratio of
    the number of tags in each dataset.}
}
\value{
  A list of elements corresponding to chromosomes, with each element
  being an $x/$y data.frame giving the position and the log2
  conservative estimate of enrichment/depletion fold ratios around that
  position. 
  Use \code{\link{writewig}} to output the structure to a WIG
  file.
}
\references{ R.M.Price, D.G. Bonett "Estimating the ratio fo two Poisson
  rates", Comp. Stat & Data Anal. 32(2000) 345}
\seealso{ \code{\link{get.smoothed.tag.density}} }
\examples{
  enrichment.estimates <- get.conservative.fold.enrichment.profile(chip.data,input.data,fws=2*binding.characteristics$whs,step=100,alpha=0.01);
  writewig(enrichment.estimates,"example.enrichment.estimates.wig","Example conservative fold-enrichment/depletion estimates shown on log2 scale");
}
