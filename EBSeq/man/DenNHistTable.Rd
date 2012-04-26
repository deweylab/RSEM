\name{DenNHistTable}
\alias{DenNHistTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
9 Density plots to compare the empirical q's and the simulated q's from the fitted beta distribution.
}
\description{
Check the beta fit of 9 different groups on isoform level data. 
}
\usage{
DenNHistTable(QList, Alpha, Beta,  AList = "F")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{QList}{
The estimated q's from the output of NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi_NoNormVar
. Input should be a list of different groups of transcripts. The number of lists here should be 9.

}

  \item{Alpha}{
The fitted parameter alpha from the output of NBBetaEB.bias.uniqueP_PoolVarSpeed
Up_MDFPoi_NoNormVar. Input should be a number if AList is not defined.
}
  \item{Beta}{
The fitted parameter beta from the output of NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi_NoNormVar. 
Input could be one single number or a vector of several numbers. The length of the input should be 9.

}

  \item{AList}{
Whether a list of alpha's are used

}

  }
\details{

}
\value{
A plot contains 9 dub-plots will be generated. 
The empirical estimation of q's will be represented as blue histogram and the density of 
the fitted beta distribution will be represented as the green line.
The main title of the plot will be "GroupName name"  
}
\references{
DenNHist, beta.mom, QQP, NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi_NoNormVar
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
# Generate Data
IsoGenerate=IsoSimu(DVDconstant=NULL, DVDqt1=.97, DVDqt2=.98, Conditions=as.factor(rep(c(1,2),each=5)), NumofSample=10, NumofIso=NULL, DEIsoProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, OnlyData="Y" )
IsoList=IsoGenerate$data

# Get Ng Vector, 5End Vector and 3End Vector
ngv=c(1,2,3,2,3,2,3,2,3)
b3v=c(1,0,0,1,1,0,0,1,1)
b5v=c(1,0,0,0,0,1,1,1,1)
NgV=unlist(sapply(1:9,function(i)rep(ngv[i],dim(IsoList[[i]])[1])))
Bias3V=unlist(sapply(1:9,function(i)rep(b3v[i],dim(IsoList[[i]])[1])))
Bias5V=unlist(sapply(1:9,function(i)rep(b5v[i],dim(IsoList[[i]])[1])))

#Run EBSeq
IsoData=do.call(rbind,IsoList)
IsoEBres=NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi_NoNormVar(Data=IsoData, NgVector=NgV, Vect5End=Bias5V, Vect3End=Bias3V, Conditions=as.factor(rep(c(1,2),each=5)),sizeFactors=rep(1,10), maxround=5)

DenNHistTable(QList=IsoEBres$QList1, Alpha=IsoEBres[[1]][5,], Beta=IsoEBres[[2]][5,], AList="F")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
