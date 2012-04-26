\name{CheckNg}
\alias{CheckNg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CheckNgStructure
}
\description{
Provide the best polynomial fit of log variance and log mean in each Ng group. 
}
\usage{
CheckNg(NewMean, NewVar,nterm, xlim, ylim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NewMean}{
A list contains  
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
4 plots will be returned. The first 3 are the data and the best fit line from 
Ng=1, Ng=2 and Ng=3 group.
The 4th plot is the scatetr plot of all the data.  
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

### Simulate Isoform Level Data
IsoGenerate=IsoSimu(DVDconstant=NULL, DVDqt1=.97, DVDqt2=.98, Conditions=as.factor(rep(c(1,2),each=5)), NumofSample=10, NumofIso=NULL, DEIsoProp=.1, Phiconstant=NULL, Phi.qt1=.4, Phi.qt2=.6, OnlyData="Y" )
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
CheckNg(IsoEBres$C1Mean, IsoEBres$C1EstVar,5, c(-1,5),c(-1,7))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Ng }
