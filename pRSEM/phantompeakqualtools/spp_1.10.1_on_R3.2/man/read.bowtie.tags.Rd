\name{read.bowtie.tags}
\alias{read.bowtie.tags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read bowtie text alignment output file }
\description{
  Reads in bowtie alignment results in text format
}
\usage{
read.bowtie.tags(filename, read.tag.names = F, fix.chromosome.names = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ bowtie text output file }
  \item{read.tag.names}{ Whether the tag names should be read in }
  \item{fix.chromosome.names}{ Whether to remove ".fa" from the end of
    the sequence names }
}
\value{
  \item{tags }{ A vector of 5' tag coordinates, with negative values
    corresponding to tags mapped to the negative strand. }
  \item{quality }{ Number of mismatches }
  \item{names }{ Tag names, if \code{read.tag.names} was set }
}