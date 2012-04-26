\name{QuantileNorm}
\alias{QuantileNorm}
\title{
Quantile Normalization
}
\description{
The Quantile normalization
}
\usage{
QuantileNorm(Data, Quantile)
}
\arguments{

  \item{Data}{
The data matrix with transcripts in rows and lanes in columns.
}
\item{Quantile}{
The quantile the user wishs to use. Should be a number between 0 and 1.	
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Use a quantile point to normalize the data.
}
\value{
The function will return a vector contains the normalization factor for each lane.
% ...
}
\references{}
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

Sizes=QuantileNorm(GeneData)
# Run EBSeq
EBres=EBTest(Data=GeneData, NgVector=rep(1,10^4), Vect5End=rep(1,10^4), Vect3End=rep(1,10^4), Conditions=as.factor(rep(c(1,2),each=5)), sizeFactors=Sizes,maxround=5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
