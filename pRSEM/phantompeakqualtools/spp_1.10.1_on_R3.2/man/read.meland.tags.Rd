\name{read.meland.tags}
\alias{read.meland.tags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read modified BED tag alignment file that contains variable match
  length information }
\description{
  Reads in an extended BED tag alignment file. An example line given  below:
  \code{49 . U1 . 1 . . 23 chr2 -234567}
  The line above specifies a 23-bp portion of the tag tag with id 49 was
  aligned with 1 mismatch to the negative strand of chr2 at position 234567.
}
\usage{
read.meland.tags(filename, read.tag.names = F, fix.chromosome.names = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ name of the extended BED file }
  \item{read.tag.names}{ whether to read in tag names }
  \item{fix.chromosome.names}{ whether to remove ".fa" from the sequence
  name ends. }
}
\value{
  \item{tags }{ A vector of 5' tag coordinates, with negative values
    corresponding to tags mapped to the negative strand. }
  \item{quality }{ Quality expressed as a float x.y, where x is
    tag.length - aligned.tag.portion.length, and y is the number of
    mismatches (must be less than 10). }
  \item{names }{ Tag names, if \code{read.tag.names} was set }
}