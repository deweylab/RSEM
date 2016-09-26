\name{remove.local.tag.anomalies}
\alias{remove.local.tag.anomalies}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Restrict or remove positions with too many tags relative to
  local background. }
\description{
  In Solexa ChIP-seq experiments some anomalous positions contain
  extremely high number of tags at the exact coordinates. The function
  scans the chromosomes, determining local tag density based on a
  provided \code{window.size}, doing two types of corrections:
  1. removing all tags from positions that exceed local density by
  \code{eliminate.fold}; 2. reducing the tag count at positions
  exceeding \code{cap.fold} to the maximal allowed count. The
  statistical significance of counts exceeding either of these two
  threshold densities is calculated based on Poisson model, with
  confidence interval determined by the \code{z.threshold} Z-score parameter.
}
\usage{
remove.local.tag.anomalies(tags, window.size = 200, eliminate.fold = 10, cap.fold = 4, z.threshold = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tags}{ Chromosome-list of tag vectors }
  \item{window.size}{ Size of the window used to assess local
    density. Increasing the window size considerably beyond the size of
    the binding features will result in flattened profiles, with bound
    positions exhibiting a difference of just 1 tag beyond the background. }
  \item{eliminate.fold}{ Threshold definining fold-over background
    density above which the position is considered anomalous and removed
  completely.}
  \item{cap.fold}{ Threshold fold-over background density above which
    the position is capped to the maximum statistically likely given
    local tag density }
  \item{z.threshold}{ Z-score used to assess significance of a given
    position exceeding either of the two density thresholds. }
}
\value{
  A modified chromosome-wise tag vector list.
}
\references{ ~put references to the literature/web site here ~ }

\note{ ~~further notes~~
  Increasing window.size to very large values will result in flat
  profiles similar to those described by Zhang et al. "Model-based
  Analysis of ChIP-Seq (MACS)." Genome Biol. 2008 Sep 17;9(9):R137.
}
