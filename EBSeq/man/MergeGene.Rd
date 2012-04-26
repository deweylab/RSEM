\name{MergeGene}
\alias{MergeGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots of gene simulation result
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MergeGene(GeneSIMout, Num, Path = "./")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GeneSIMout}{
The output of GeneSimu with OnlyData="F".
}
  \item{Num}{
How many times the simulation ran.
}
  \item{Path}{
The path to store the plots
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
IsoMerge
}
\examples{
GeneGenerate=GeneSimu(DVDconstant=4, DVDqt1=NULL, DVDqt2=NULL, Conditions=rep(c(1,2),each=5), NumofSample=10, NumofGene=10000, DEGeneProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, Meanconstant=NULL, OnlyData="F")

GeneTable=MergeGene(GeneGenerate,1,"./")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
