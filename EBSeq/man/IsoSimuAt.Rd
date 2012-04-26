\name{IsoSimuAt}
\alias{IsoSimuAt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Isoform level simulation with outliers
}
\description{
Simulate isoform level expression data from a Negative Binomial assumption. (With outliers)
}
\usage{
IsoSimuAt(DVDconstant = NULL, DVDqt1 = NULL, DVDqt2 = NULL, Conditions, NumofSample, NumofIso = NULL, DEIsoProp, Phiconstant = NULL, Phi.qt1 = NULL, Phi.qt2 = NULL, NormFactor = NULL, OnlyData = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DVDconstant}{
Whether want to use constant fold change value for all the DE genes.
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
Number of samples the user want to generate.
}
  \item{NumofIso}{
Input should be a vector with length 3. All values should be non-negative.
The ith value represents how many isoforms the user want to generate for isoform group i. 
}
  \item{DEIsoProp}{
The proportion of isoforms to be generated as DE. The value should be in [0, 1].
Besides, the same proportion of isoforms will be generated as EE isoforms with outlier.
The genes will be generated as EE at first, then the count of one of the samples
(randomly selected) will be setted as its original counts multiplied by one of (4, 6, 8, 10).

}
  \item{Phiconstant}{
Whether set the disperse parameter phi to be a constant. If this parameter is specified, the settings of Phi.qt1 and Phi.qt2 will be ignored.
Input should be a vector with length 3. The ith value indicates the overdisperse parameter of isoform group i.
}
  \item{Phi.qt1, Phi.qt2}{
If Phiconstant is not specified, the user could use a range of empirical phi's from each group of Gould' data. The lower and upper bound ( quantile) could be specified.
The suggested value is c(.25, .75). phi for each gene will be randomly choosed w
ithin the range.

}
  \item{NormFactor}{
Wether set the mean of each isoform to be a constant.
}
  \item{OnlyData}{
Wether the user only want the generated data matrix. If OnlyData = T, the function will return the simulated matrix
and the name of the DE genes.
Otherwise the funtion will run DESeq, EBSeq, edgeR, baySeq and BBSeq and provide the results of each method.
Currently only OnlyData=T is supported
}
}
\details{
For each isoform, we assumed that the expression follows a Negative Binomial distribution with mean mu_gi and variance mu_gi * (1 + mu_gi * phi_gi).
For DE genes, we assumed that in one condition the genes are with mean mu_gi * DVD.
mu, phi and DVD could be specified by the parameter settings.

}
\value{
\item{data}{
A list of expression values will be generated. 
Each list represents a group of isoforms.
Group1: Ng1
Group2: Ng2
Group3: Ng3
The rows refer to the isoforms and the columns are the samples. 
The isoforms are named "I_GroupNumber_IsoformNumber". The first part of the isoforms of each group will be the DE ones. (The number depends on the DEIsoProp parameter.)
}
\item{TrueDE}{The names of the isoforms who are defined to be DE.}

\item{Outliers}{The names of the genes who are defined to be outliers at each level of (4,6,8,10).}
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
GeneSimu, IsoSimu
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
IsoGenerate=IsoSimuAt(DVDconstant=NULL, DVDqt1=.97, DVDqt2=.98, Conditions=as.factor(rep(c(1,2),each=5)), NumofSample=10, NumofIso=NULL, DEIsoProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, OnlyData=T )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
