\name{crit_fun}
\alias{crit_fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the adjusted FDR threshold 
}
\description{
Calculate the adjusted FDR threshold using the posterior probabilities at a target FDR
}
\usage{
crit_fun(PPEE, thre)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PPEE}{
The posterior probabilities of being EE.
}
  \item{thre}{
The target FDR.
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The adjusted FDR threshold of target FDR.
}
\references{
}
\author{
Ning Leng
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ FDR }
