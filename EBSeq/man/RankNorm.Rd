\name{RankNorm}
\alias{RankNorm}
\title{
Rank Normalization
}
\description{
}
\usage{
RankNorm(Data)
}
\arguments{

  \item{Data}{
The data matrix with transcripts in rows and lanes in columns.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function will return a matrix contains the normalization factor for each lane and each transcript.
% ...
}
\references{
}
\author{
Ning Leng
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
GeneGenerate=GeneSimu(DVDconstant=4, DVDqt1=NULL, DVDqt2=NULL, Conditions=rep(c(1,2),each=5), NumofSample=10, NumofGene=10000, DEGeneProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, Meanconstant=NULL, OnlyData="Y")
GeneData=GeneGenerate$data

Sizes=RankNorm(GeneData)
# Run EBSeq
EBres=EBTest(Data=GeneData, NgVector=rep(1,10^4), Vect5End=rep(1,10^4), Vect3End=rep(1,10^4), Conditions=as.factor(rep(c(1,2),each=5)), sizeFactors=Sizes,maxround=5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
