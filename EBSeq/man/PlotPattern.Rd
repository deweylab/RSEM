\name{PlotPattern}
\alias{PlotPattern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize the patterns
}
\description{
visualize the patterns
}
\usage{
PlotPattern(PosParti)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PosParti}{
The output of GetPatterns function.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A heatmap to visualize the patterns of interest.
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
Conditions=c("C1","C1","C2","C2","C3","C3")
PosParti=GetPatterns(Conditions)
PlotPattern(PosParti)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
