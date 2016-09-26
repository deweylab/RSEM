\name{output.binding.results}
\alias{output.binding.results}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write out determined binding peaks into a text file table }
\description{
  Writes out determined binding positions into a text file. The file
  will contain a table with each row corresponding to a detected
  position, with the following columns: 
  \item{chr}{ chromosome or target sequence }
  \item{pos}{ position of detected binding site on the chromosome/sequence}
  \item{score}{a score reflecting magnitude of the binding}
  \item{Evalue}{E-value corresponding to the peak magnitude}
  \item{FDR}{FDR corresponding to the peak magnitude}
  \item{enrichment.lb}{lower bound of the fold-enrichment ratio}
  \item{enrichment.mle}{maximum likelihood estimate of the fold-enrichment ratio}
}
\usage{
output.binding.results(results, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results}{ output of the \code{\link{find.binding.positions}} }
  \item{filename}{ file name }
}
