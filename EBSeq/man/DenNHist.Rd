\name{DenNHist}
\alias{DenNHist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Density plot to compare th empirical q's and the simulated q's from the fitted beta distribution.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DenNHist(QList, Alpha, Beta, name, AList = "F", GroupName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{QList}{
The estimated q's from the output of NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi_NoNormVar
. Input could be a vector or a list of different groups of transcripts. The number of lists here should be the same as the length of Beta.

}
  \item{Alpha}{
The fitted parameter alpha from the output of NBBetaEB.bias.uniqueP_PoolVarSpeed
Up_MDFPoi_NoNormVar. Input should be a number if AList is not defined.
}
  \item{Beta}{
The fitted parameter beta from the output of NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi_NoNormVar. 
Input could be one single number or a vector of several numbers. The length of the input should be the same as the number of lists of QList.

}
  \item{name}{
The name of the plots
}
  \item{AList}{
Whether a list of alpha's are used

}
  \item{GroupName}{
The names of each sub plot. The l
ength of the input should be the same as the number of lists of QList.

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plots will be generated. Each plot represents a sub-list of the QList.
The empirical estimation of q's will be represented as blue histogram and the density of
the fitted beta distribution will be represented as the green line.
The main title of the plot will be "GroupName name".

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
beta.mom, DenNHistTable, QQP, NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi_NoNormVar
}
\examples{
### Simulate Gene Level Data
GeneGenerate=GeneSimu(DVDconstant=4, DVDqt1=NULL, DVDqt2=NULL, Conditions=rep(c(1,2),each=5), NumofSample=10, NumofGene=10000, DEGeneProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, Meanconstant=NULL, OnlyData="Y")
GeneData=GeneGenerate$data

# Run EBSeq
EBres=NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi_NoNormVar(Data=GeneData, NgVector=rep(1,10^4), Vect5End=rep(1,10^4), Vect3End=rep(1,10^4), Conditions=as.factor(rep(c(1,2),each=5)), sizeFactors=rep(1,10),maxround=5)
# Plot
DenNHist(QList=EBres$QList1, Alpha=EBres[[1]][5,1], Beta=EBres[[2]][5,1], name="Gene", AList="F", GroupName="")

### Simulate Isoform Level Data
IsoGenerate=IsoSimu(DVDconstant=NULL, DVDqt1=.97, DVDqt2=.98, Conditions=as.factor(rep(c(1,2),each=5)), NumofSample=10, NumofIso=NULL, DEIsoProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, OnlyData="Y" )
IsoList=IsoGenerate$data

# Get Vectors and Run EBSeq
ngv=c(1,2,3,2,3,2,3,2,3)
b3v=c(1,0,0,1,1,0,0,1,1)
b5v=c(1,0,0,0,0,1,1,1,1)
NgV=unlist(sapply(1:9,function(i)rep(ngv[i],dim(IsoList[[i]])[1])))
Bias3V=unlist(sapply(1:9,function(i)rep(b3v[i],dim(IsoList[[i]])[1])))
Bias5V=unlist(sapply(1:9,function(i)rep(b5v[i],dim(IsoList[[i]])[1])))

IsoData=do.call(rbind,IsoList)
IsoEBres=NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi_NoNormVar(Data=IsoData, NgVector=NgV, Vect5End=Bias5V, Vect3End=Bias3V, Conditions=as.factor(rep(c(1,2),each=5)),sizeFactors=rep(1,10), maxround=5)

# Plot
par(mfrow=c(3,3))
DenNHist(QList=IsoEBres$QList1, Alpha=IsoEBres[[1]][5,], Beta=IsoEBres[[2]][5,], name="Isoform", AList="F", GroupName=paste("group",c(1:9),sep=""))







}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ beta }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
