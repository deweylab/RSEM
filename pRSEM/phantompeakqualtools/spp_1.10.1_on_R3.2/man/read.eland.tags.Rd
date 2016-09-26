\name{read.eland.tags}
\alias{read.eland.tags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read eland output file }
\description{
  Reads in ELAND output file, returning 5'-end tag coordinates and
  number of mismatches associated with each mapped tag.
}
\usage{
read.eland.tags(filename, read.tag.names = F, fix.chromosome.names = T, max.eland.tag.length = -1,extended=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ ELAND output file }
  \item{read.tag.names}{ Whether the tag names should be read in }
  \item{fix.chromosome.names}{ Whether to remove ".fa" from the end of
    the sequence names }
  \item{max.eland.tag.length}{ Specifies max length of the tag sequence
    considered by ELAND. This needs to be specified if the tags are
    longer than the sequences considred by ELAND during alignment. }
  \item{extended}{ Whether the file is written out in "extended" format
    provided in GA pipeline 1.0. }
  \item{multi}{ Whether the file is written in "multi" format, showing multiple alignments of the reads }
}
\value{
  \item{tags }{ A vector of 5' tag coordinates, with negative values
    corresponding to tags mapped to the negative strand. }
  \item{quality }{ Number of mismatches }
  \item{names }{ Tag names, if \code{read.tag.names} was set }
}
