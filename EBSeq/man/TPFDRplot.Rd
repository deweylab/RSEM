\name{TPFDRplot}
\alias{TPFDRplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the number of top counts vs FDR for each method in simulation data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TPFDRplot(DESeqP, EBZ, TrueDE, main, FDR = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DESeqP}{
%%     ~~Describe \code{DESeqP} here~~
}
  \item{EBZ}{
%%     ~~Describe \code{EBZ} here~~
}
  \item{TrueDE}{
%%     ~~Describe \code{TrueDE} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{FDR}{
%%     ~~Describe \code{FDR} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(DESeqP, EBZ, TrueDE, main, FDR=NULL){
	Seq=seq(0.001,0.5,by=0.001)
	DETPR=rep(0,length(Seq))
	EBTPR=rep(0,length(Seq))
	DEFDR=rep(0,length(Seq))
	EBFDR=rep(0,length(Seq))
	DETPNum=rep(0,length(Seq))
    EBTPNum=rep(0,length(Seq))
    DEFDNum=rep(0,length(Seq))
    EBFDNum=rep(0,length(Seq))
	for (i in 1:length(Seq)){
		DESeqOnes=names(DESeqP)[DESeqP<=Seq[i]]
		if (length(FDR)==0) EBOnes=names(EBZ)[EBZ>=crit.fun(1-EBZ, Seq[i])]
		else if (FDR=="H") EBOnes=names(EBZ)[EBZ>=(1-Seq[i])]
			else EBOnes=names(EBZ)[EBZ>=FDR[i]]

		DETPNum[i]=sum(DESeqOnes\%in\%TrueDE)
		EBTPNum[i]=sum(EBOnes\%in\%TrueDE)
		DEFDNum[i]=sum(!DESeqOnes\%in\%TrueDE)
		EBFDNum[i]=sum(!EBOnes\%in\%TrueDE)
		
		DETPR[i]=DETPNum[i]/length(TrueDE)
		EBTPR[i]=EBTPNum[i]/length(TrueDE)
		DEFDR[i]=DEFDNum[i]/length(TrueDE)
		EBFDR[i]=EBFDNum[i]/length(TrueDE)
	}
	plot(Seq,DETPR,ylim=c(0,1),xlim=c(0,.5),type="l",col="red", main=paste(main, "TPR"),xlab="controled FDR level", ylab="TPR",lwd=2)
	lines(Seq,EBTPR,col="blue",lwd=2)
	legend("bottomright",lwd=2, col=c("red","blue"), c("DESeq","EBSeq"))

	plot(Seq,DEFDR,ylim=c(0,1),xlim=c(0,.5),type="l",col="red", main=paste(main, "FDR"),xlab="controled FDR level", ylab="TPR",lwd=2)
	lines(Seq,EBFDR,col="blue",lwd=2)
	legend("topleft", lwd=2, col=c("red","blue"), c("DESeq","EBSeq"))


	output=cbind( DETPR,EBTPR, DEFDR,EBFDR,DETPNum,EBTPNum,DEFDNum,EBFDNum)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
