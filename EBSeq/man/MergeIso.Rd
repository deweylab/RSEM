\name{MergeIso}
\alias{MergeIso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Isoforms of gene simulation result
}
\description{
}
\usage{
MergeIso(IsoSIMout, Num, Path = "./")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IsoSIMout}{
The output of IsoSimu with OnlyData="F".
}
  \item{Num}{
How many times the simulation ran.

}
  \item{Path}{
	  The path to store the plots.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
3 plots will be generated.
1 FPR vs TPR of each method
2 FDR vs TPR of each method
2 Top counts vs FDR of each method

A table will be generated which contains the FDR and TPR of each method.
Each method will be ran on all the data and within group.
(Using p-value=.05 or Posterior Probability=.95).


}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ning Leng
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
MergeGene
}
\examples{
IsoGenerate=IsoSimu(DVDconstant=NULL, DVDqt1=.97, DVDqt2=.98, Conditions=as.factor(rep(c(1,2),each=5)), NumofSample=10, NumofIso=NULL, DEIsoProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, OnlyData="F" )

IsoTable=MergeIso(IsoGenerate,1,"./")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
