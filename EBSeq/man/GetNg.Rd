\name{GetNg}
\alias{GetNg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate the Ng vector
}
\description{
Generate the Ng vector for the isoforms
}
\usage{
GetNg(IsoformName, GeneName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IsoformName}{
A vector contains the isoform names
}
  \item{GeneName}{
The gene names of the isoforms in IsoformNames. (Should be in the same order)
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{GeneNg}{
  The number of isoforms each gene contains
  }
  \item{GeneNgTrun}{
  The truncated Ng of each gene. (The genes contain more than 3 isoforms are with Ng 3.)
  }
   \item{IsoformNg}{
  The Ng of each isoform 
  }
    \item{IsoformNgTrun}{
   The truncated Ng of each isoform. 
  }


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
IsoGenerate=IsoSimu(DVDconstant=NULL, DVDqt1=.97, DVDqt2=.98, Conditions=as.factor(rep(c(1,2),each=5)), NumofSample=10, NumofIso=c(1000,2000,3000), DEIsoProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, OnlyData=T )

IsoMat=do.call(rbind,IsoGenerate$data)
IsoNames=rownames(IsoMat)

Ngvector=GetNg(IsoNames, IsosGeneNames)

IsoEBres=EBTest(Data=IsoMat, NgVector=IsoNgTrun, Conditions=as.factor(rep(c(1,2),each=5)),sizeFactors=rep(1,10), maxround=5)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Ng }
