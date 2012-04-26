\name{GetMultiPP}
\alias{GetMultiPP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate the Posterior Probability of each transcript.
}
\description{
Generate the Posterior Probability of being each pattern of each transcript based on the EBMultiTest output.
}
\usage{
GetMultiPP(EBout)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EBout}{
The output of EBMultiTest function.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{PP}{The poster probabilities of being each pattern.}
\item{MAP}{The most likely pattern each gene belongs to}
\item{Patterns}{The Patterns}
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
\keyword{ Posterior Probability }
