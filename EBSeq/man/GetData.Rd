\name{GetData}
\alias{GetData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read in RSEM output of Gould data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GetData(path, Name1, Name2, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
The path of RSEM outputs
}
  \item{Name1}{
The output names of the files from Condition 1
}
  \item{Name2}{
The output names of the files from Condition 2
}
  \item{type}{
If type="G", read in the gene level output
If type="I", read in the isoform level output
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output is the "nu values" from RSEM.
To generate a expression matrix, the user need to run the PoolMatrix function.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(path,Name1,Name2,type)
{
Data=vector("list",8)
Filenames=NULL
Tablenames=NULL
for (name in 1:4)
	{
		if (type=="I")
			Filenames=c(Filenames,paste(path,Name1,name,"_isoform_nus.tab",sep=""))  
		if (type=="G")	
			Filenames=c(Filenames,paste(path,Name1,name,"_gene_nus.tab",sep=""))  
		Tablenames=c(Tablenames,paste(Name1,name,sep=""))
	}
for (name in 1:4)
	{
		if (type=="I")
			Filenames=c(Filenames,paste(path,Name2,name,"_isoform_nus.tab",sep=""))
		if (type=="G")
			Filenames=c(Filenames,paste(path,Name2,name,"_gene_nus.tab",sep=""))
		Tablenames=c(Tablenames,paste(Name2,name,sep=""))
	}


names(Data)=Tablenames
for (file in 1:8)
	{
		temp=read.table(Filenames[file],header=T)
		temp2=as.matrix(temp[-1])
		rownames(temp2)=as.vector(as.matrix(temp[1]))
		Data[[file]]=temp2
	}
	Data
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
