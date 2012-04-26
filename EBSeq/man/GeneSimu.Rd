\name{GeneSimu}
\alias{GeneSimu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gene Level Simulation 
}
\description{
Simulate Gene level expression data from a Negative Binomial assumption. (Without outliers)
}
\usage{
GeneSimu(DVDconstant = NULL, DVDqt1 = NULL, DVDqt2 = NULL, Conditions, NumofSample, NumofGene = NULL, DEGeneProp, Phiconstant = NULL, Phi.qt1 = NULL, Phi.qt2 = NULL, Meanconstant = NULL,  NormFactor=NULL, OnlyData = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DVDconstant}{
Whether want to use constant fold change value for all the DE genes. 
If set DVDconstant=4, all the DE genes will have fold change of 4 across two condtions. 
If this parameter is specified, the settings of DVDqt1 and DVDqt2 will be ignored. 
}
  \item{DVDqt1, DVDqt2}{
If DVDconstant is not specified, the user could use a range of empirical DVD's f
rom Gould' data. The lower and upper bound ( quantile) could be specified.
The suggested value is c(.96, .97). DVD for each gene will be randomly choosed within the range.

}
  \item{Conditions}{
A vector of charecters to show each sample's condition. 
(Only the two-condition case is supported now)
}
  \item{NumofSample}{
Number of samples to generte.
}
  \item{NumofGene}{
Number of genes to generate.
}
  \item{DEGeneProp}{
The proportion of genes to be generated as DE. The value should be in [0, 1].
Besides, the same proportion of genes will be generated as EE genes with outlier. 
The genes will be generated as EE at first, then the count of one of the samples 
(randomly selected) will be setted as its original counts multiplied by one of (4, 6, 8, 10).
}
  \item{Phiconstant}{
Whether set the disperse parameter phi to be a constant. If this parameter is specified, the settings of Phi.qt1 and Phi.qt2 will be ignored.
}
  \item{Phi.qt1, Phi.qt2}{
If Phiconstant is not specified, the user could use a range of empirical phi's from Gould' data. The lower and upper bound ( quantile) could be specified.
The suggested value is c(.25, .75). phi for each gene will be randomly choosed w
ithin the range.

}
  \item{Meanconstant}{
Wether set the mean of each gene to be a constant.
}
  \item{OnlyData}{
Wether the user only want the generated data matrix. If OnlyData=T, the function will return the simulated matrix
and the name of the DE genes.
Otherwise the funtion will run DESeq, EBSeq, edgeR, baySeq and BBSeq and provide the results of each method. 
}

\item{NormFactor}{
If NormFactor is NULL, each lane will be set to be with the same library size. Otherwise NormFactor should be a 
vector with length NumofSample. 
}

}
\details{
For each gene, we assumed that the expression follows a Negative Binomial distribution with mean mu_g and variance mu_g * (1 + mu_g * phi_g). 
For DE genes, we assumed that in one condition the genes are with mean mu_g * DVD.
mu, phi and DVD could be specified by the parameter settings.
}
\value{
\item{data}{
A matrix of expression values will be generated. The rows of the matrix refer to the genes and the columns of the matrix are the samples. The genes are named "G_1", "G_2", ... The first part of the genes will be the DE ones. (The number depends on the DEGeneProp parameter.)
}
\item{TrueDE}{The names of the genes who are defined to be DE.}
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
IsoSimu, IsoSimuAt, GeneSimuAt
}
\examples{
GeneGenerate=GeneSimu(DVDconstant=4, DVDqt1=NULL, DVDqt2=NULL, Conditions=rep(c(1,2),each=5), NumofSample=10, NumofGene=10000, DEGeneProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, Meanconstant=NULL, NormFactor=rep(c(.9,1.2),5),OnlyData=T)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ simulation }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
