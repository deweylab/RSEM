\name{beta.mom}
\alias{beta.mom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit the beta distribution by method of moments
}
\description{
Fit the beta distribution by method of moments
}
\usage{
beta.mom(qs.in)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{qs.in}{
A vector contains the numbers that are assumed to follow a beta distribution 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{alpha.hat }{Return the estimation of alpha}
  \item{beta.hat}{Return the estimation of beta}
}
\references{
}
\author{
Ning Leng
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
DenNHist, DenNHistTable
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
tmp=rbeta(5,5,100)
param=beta.mom(tmp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ beta }
