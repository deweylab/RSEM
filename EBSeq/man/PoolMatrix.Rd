\name{PoolMatrix}
\alias{PoolMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate the expression matrix from the output of GetData
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PoolMatrix(Data, reads, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
The output from GetData function.
}
  \item{reads}{
The total number of reads in each lane. Could be obtained from the RSEM outputs. 
}
  \item{type}{
If type="S", the outputs will be the a matrix which transcript names in row and sample names in column.
If type="G", the first column will be the group information. 

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{PoolM}{The matrix of nu values}
\item{PoolValue}{The matrix of expression values}
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
GetData
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(Data,reads,type)
{
poolnames=names(Data)
poolM=NULL
for (po in 1:8)
	poolM=cbind(poolM,Data[[po]][,1])
rownames(poolM)=rownames(Data[[1]])
colnames(poolM)=poolnames

#poolValue=poolM*reads
poolValue=poolM
for (col in 1:8)
	poolValue[,col]=poolM[,col]*reads[col]
poolValue=round(poolValue)
if (type=="G")
	{
		poolM=cbind(Data[[1]][,2],poolM)
		poolValue=cbind(Data[[1]][,2],poolValue)
		colnames(poolM)=c("Groups",poolnames)
		colnames(poolValue)=c("Groups",poolnames)
	}
poolOutput=list(poolM=poolM,poolValue=poolValue)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
