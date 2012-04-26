\name{TopCts}
\alias{TopCts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get FDR of Top Counts 
}
\description{

}
\usage{
TopCts(pvalue, PP = NULL, TrueNames, TopNum)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalue}{
A matrix contains the p values (Posterior probabilities) for each transcript and each method.
Rows are for different methods and columns are for different transcripts.
}
  \item{PP}{
The length of PP vector should be the same as the number of columns in pvalue matrix.
The value in PP either 0 or 1.
If the ith value of PP is 0, it means the ith method (the ith row of pvalue) provided p-values.
If the ith value of PP is 1, it means the ith method (the ith row of pvalue) provided posterior probabilities. 
}
  \item{TrueNames}{
The names of the transcripts who defined to be DE.
}
  \item{TopNum}{
The number of top counts we are interested in.
For example, if TopNum=1000, we'll calculate the FDR's of each method if we pick the top 1, 2, ... 1000 genes.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A metrix contains the FDR's.
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
function(pvalue, PP=NULL, TrueNames, TopNum){
	NumOfMethods=ncol(pvalue)
	puse=pvalue
	if(1\%in\%PP)puse[,PP==1]=1-pvalue[,PP==1]
	#puse.list=data.frame(puse)
	FD=matrix(rep(0,NumOfMethods*TopNum),ncol=NumOfMethods)
#	Rank=apply(puse,2,rank)
#	for(i in 1:TopNum)
#		FD[i,]=sapply(1:NumOfMethods, function(j)sum(!rownames(Rank)[Rank[,j]<=i]\%in\%TrueNames))	
#	FD=sapply(1:TopNum, function(i)sapply(1:NumOfMethods, function(j)sum(!rownames(Rank)[Rank[,j]<=i]\%in\%TrueNames)))
	for (s in 1:NumOfMethods){
		tmp=puse[,s]
		names(tmp)=rownames(puse)
		sorttmp=sort(tmp)
		for( c in 2:TopNum)
			FD[c, s]=FD[(c-1),s]+as.numeric(!names(sorttmp)[c]\%in\%TrueNames)
	}
	FD
	#matplot(TopNum,FD,type="l",ylim=c(0,1),xlab="Top DE selected", ylab="FDR")
	#legend("rightbottom",col=1:TopNum, lty=1:TopNum, names)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
