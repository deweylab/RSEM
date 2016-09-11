\name{write.narrowpeak.binding}
\alias{write.narrowpeak.binding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write out determined binding peaks using narrowPeak format }
\description{
  Writes out determined binding positions into a narrowPeak file.
  The region will correspond to associated broad enrichment region, if
  such were added using add.broad.peak.regions method. Otherwise the
  region size will be determined using margin (which defaults to the
  window half size that was used to determine binding positions)
}
\usage{
write.narrowpeak.binding(results, filename,margin=results$whs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results}{ output of the \code{\link{find.binding.positions}} }
  \item{filename}{ file name }
  \item{margin}{ explicit value of the margin to be used if the borad
    region information is absent (defaults to peak detection window half-size}
}
