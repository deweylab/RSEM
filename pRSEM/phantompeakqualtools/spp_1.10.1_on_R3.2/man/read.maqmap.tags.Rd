\name{read.maqmap.tags}
\alias{read.maqmap.tags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read MAQ text alignment output file }
\description{
  Reads in MAQ alignment results in text format (that results from "maq mapview" command.)
}
\usage{
read.maqmap.tags(filename, read.tag.names = F, fix.chromosome.names = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ MAQ text output file }
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