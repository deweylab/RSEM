\name{writewig}
\alias{writewig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to save a list of chromosome-wise x/y data frames
  into a WIG file format. }
\description{
  Takes a list that contains an $x and $y data.frame for a number of
  chromosomes and writes it out to a WIG BED style format.
}
\usage{
writewig(dat, fname, feature, threshold = 5, zip = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ Chromosome coordinate-value data. \code{dat} is a list,
    each member of a list is a data frame with $x and $y columns
    containing chromosome positions and associated values. The names of
    the list elements correspond to the chromosomes. }
  \item{fname}{ Filename to which the output should be written }
  \item{feature}{ Data description to be incorporated into the WIG header }
  \item{threshold}{ Optional threshold to be saved in the WIG file}
  \item{zip}{ Wheter to invoke a zip program to compress the file }
}

\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

data <- list("chr1"=data.frame(x=c(100,130,200),y=c(1.2,4.0,2.3)));
writewig(data,"filename");

}
