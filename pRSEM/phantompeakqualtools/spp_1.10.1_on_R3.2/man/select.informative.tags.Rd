\name{select.informative.tags}
\alias{select.informative.tags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Choose informative tags }
\description{
  For datasets with tag alignment quality information (e.g. number of
  mismatches for Eland alignments),
  \code{\link{get.binding.characteristics}} determines whether inclusion
  of tags from each specific quality bin improves the cross-correlation
  profile. The present function is then used to actually select these
  informative tags, discarding all other information, including quality
  scores that are not used in further processing.
}
\usage{
select.informative.tags(data, binding.characteristics)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Full alignment data (a list with $tags and $quality elements) }
  \item{binding.characteristics}{ result of a
    \code{\link{get.binding.characteristics}}  call. If NULL value is
    supplied,all tags will be accepted. }
}
\value{
  A chromosome-wise tag list. Each element of the list corresponds to a
  chromosome and is a numeric vector of 5' tag coordinates, with sign
  designating DNA strand.
  This form of tag data is used for most of the other processing.
}
