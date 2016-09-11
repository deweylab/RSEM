\name{read.bam.tags}
\alias{read.bam.tags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read BAM alignment file }
\description{
  Reads in aligned reads from BAM file. Note: no split (non-unique)
  alignemnts should be reported in the BAM file.
}
\usage{
read.bam.tags(filename, read.tag.names = F, fix.chromosome.names = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ BAM file }
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