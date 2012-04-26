\name{QQP}
\alias{QQP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The QQ Plot of empirical q's and simulated q's from fitted beta distribution
}
\description{

}
\usage{
QQP(QList, AlphaResult, BetaResult, name, AList="F", GroupName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{QList}{
The estimated q's from the output of NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi. Input could be a vector or a list of different groups of transcripts. The number of lists here should be the same as the length of BetaResult.

}
  \item{AlphaResult}{
The fitted parameter alpha from the output of NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi. Input should be a number if AList is not defined.
}
  \item{BetaResult}{
The fitted parameter beta from the output of NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi. Input could be one single number or a vector of several numbers. The length of the input should be the same as the number of lists of QList.
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
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
 NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi , DenNHist
}
\examples{
GeneData=GeneSimu(DVDconstant=4, DVDqt1=NULL, DVDqt2=NULL, Conditions=rep(c(1,2),each=5), NumofSample=10, NumofGene=10000, DEGeneProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, Meanconstant=NULL, OnlyData="Y")

EBres=NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi(Data=GeneData, NgVector=rep(1,10^4), Vect5End=rep(1,10^4), Vect3End=rep(1,10^4), Conditions=as.factor(rep(c(1,2),each=5)), maxround=5)

QQP(QList=EBres$QList1, AlphaResult=EBres[[1]][5,1], BetaResult=EBres[[2]][5,1], name="Gene", AList="F", GroupName=NULL)

## The function is currently defined as
function(QList,AlphaResult,BetaResult,name,AList="F",GroupName){
		    for (i in 1:length(BetaResult)){
				tmpSize=length(QList[[i]][QList[[i]]<1])
			if (AList=="F") rdpts=rbeta(tmpSize,AlphaResult,BetaResult[i])
				else rdpts=rbeta(tmpSize,AlphaResult[i],BetaResult[i])
	qqplot(QList[[i]][QList[[i]]<1], rdpts, xlab="estimated q's", "simulated q's from fitted beta distribution",main=paste(names(name,GroupName[i],sep=" "),xlim=c(0,1),ylim=c(0,1))
			}
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
