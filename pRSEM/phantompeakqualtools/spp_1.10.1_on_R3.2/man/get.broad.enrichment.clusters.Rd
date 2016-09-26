\name{get.broad.enrichment.clusters}
\alias{get.broad.enrichment.clusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Determine broad clusters of enrichment }
\description{
  Scan chromosomes with a pre-defined window size, comparing scaled ChIP
  and input tag coutns to see if their ratio exceeds that expected from
  a Poisson process (normalized for dataset size).
}
\usage{
get.broad.enrichment.clusters(chip.tags, input.tags, window.size=1e3,z.thr=3,tag.shift=146/2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chip.tags}{ foreground tag vector list }
  \item{input.tags}{ background tag vector list }
  \item{window.size}{ window size to be used for tag counting }
  \item{z.thr}{ Z-score to be used as a significance threshold }
  \item{tag.shift}{ number of base pairs by which positive and negative
    tag coordinates should be shifted towards eachother (half of binding
  peak separation distance)}
}
\value{
  A list of elements corresponding to chromosomes, with each element
  being an $s/$e/$rv data.frame giving the starting, ending positions and the log2
  enrichment estimate for that region.
}
